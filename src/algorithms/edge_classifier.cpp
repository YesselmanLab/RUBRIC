/**
 * @file edge_classifier.cpp
 * @brief LW edge classification.
 */
#include <pairfinder/algorithms/edge_classifier.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>

namespace pairfinder::algorithms::classification {

namespace {

using core::BaseTyping;
using core::Residue;
using core::ReferenceFrame;
using geometry::Vector3d;
using hbond::HBond;

// ---- tuning constants ---------------------------------
constexpr double kPi = 3.14159265358979323846;
constexpr double kN7GuardMargin = 0.2;
constexpr double kVotingContactDistance = 4.0;
constexpr double kPyrWsTiebreakerPlaneDeg = 45.0;
constexpr double kPyrWsTiebreakerMarginA = 0.15;
constexpr double kBothBoundaryAmbigMargin = 0.15;
constexpr double kPyrWsGap = 0.20;
constexpr double kCentroidAmbiguityGap = 0.15;
constexpr double kRingBondMarginThreshold = 0.5;
constexpr double kCosThreshold = 0.7;

const std::unordered_set<std::string> kBaseHbondAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::unordered_set<std::string> kBaseRingAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6",
    "C2", "C4", "C5", "C6", "C8"};
const std::unordered_set<std::string> kDnaBaseCodes = {"DA", "DT", "DC", "DG", "DU", "T"};

int face_priority(char f) { return f == 'W' ? 0 : (f == 'H' ? 1 : 2); }

// ---- per-base tables (resolved to single-letter key A/G/C/U/T/I) -----------

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

bool is_known_key(const std::string& k) {
    return k == "A" || k == "G" || k == "C" || k == "U" || k == "T" || k == "I";
}

// Resolve a base_type to a known single-letter key, or "" if none.
std::string resolve_key(const std::string& base, const BaseTyping& typing) {
    const std::string u = upper(base);
    if (is_known_key(u)) return u;
    const std::string p = typing.normalize(u);
    if (is_known_key(p)) return p;
    return "";
}

using FaceList = std::vector<std::pair<char, std::vector<std::string>>>;  // W,H,S order

FaceList face_atoms_table(const std::string& k) {
    if (k == "G") return {{'W', {"N1", "N2", "O6"}}, {'H', {"N7", "C8", "O6"}}, {'S', {"O2'", "N3"}}};
    if (k == "A") return {{'W', {"N1", "N6"}}, {'H', {"N7", "C8", "N6"}}, {'S', {"O2'", "N3"}}};
    if (k == "C") return {{'W', {"N3", "N4", "O2"}}, {'H', {"C5", "C6", "N4"}}, {'S', {"O2'", "O2"}}};
    if (k == "U" || k == "T") return {{'W', {"N3", "O4", "O2"}}, {'H', {"C5", "C6", "O4"}}, {'S', {"O2'", "O2"}}};
    if (k == "I") return {{'W', {"N1", "O6"}}, {'H', {"N7", "C8", "O6"}}, {'S', {"O2'", "N3"}}};
    return {};
}

// Face atoms for a base: fallback to parent, then to "A".
FaceList face_atoms_for(const std::string& base, const BaseTyping& typing) {
    std::string k = resolve_key(base, typing);
    if (k.empty()) k = "A";
    return face_atoms_table(k);
}

// Unambiguous face atoms for a resolved key (empty map if unknown).
std::unordered_map<std::string, char> unambiguous_face(const std::string& base,
                                                       const BaseTyping& typing) {
    const std::string k = resolve_key(base, typing);
    if (k == "C" || k == "U" || k == "T")
        return {{"N3", 'W'}, {"C5", 'H'}, {"C6", 'H'}, {"O2'", 'S'}};
    if (k == "G" || k == "A" || k == "I")
        return {{"N1", 'W'}, {"N7", 'H'}, {"C8", 'H'}, {"O2'", 'S'}, {"N3", 'S'}};
    return {};
}

// Boundary atoms for a resolved key.
std::unordered_map<std::string, std::vector<char>> boundary_atoms(const std::string& base,
                                                                  const BaseTyping& typing) {
    const std::string k = resolve_key(base, typing);
    if (k == "G") return {{"N2", {'S', 'W'}}, {"O6", {'H', 'W'}}};
    if (k == "A") return {{"N6", {'H', 'W'}}};
    if (k == "C") return {{"O2", {'S', 'W'}}, {"N4", {'H', 'W'}}};
    if (k == "U" || k == "T") return {{"O2", {'S', 'W'}}, {"O4", {'H', 'W'}}};
    if (k == "I") return {{"O6", {'H', 'W'}}};
    return {};
}

// Ring bonds (atom1, atom2, face). Purine vs pyrimidine chosen by caller.
const std::vector<std::tuple<std::string, std::string, char>>& ring_bonds(bool purine) {
    static const std::vector<std::tuple<std::string, std::string, char>> pur = {
        {"N1", "C2", 'W'}, {"C6", "N1", 'W'}, {"C5", "C6", 'H'}, {"N7", "C8", 'H'},
        {"C5", "N7", 'H'}, {"C4", "C5", 'H'}, {"N3", "C4", 'S'}, {"C4", "N9", 'S'},
        {"C2", "N3", 'S'}};
    static const std::vector<std::tuple<std::string, std::string, char>> pyr = {
        {"N3", "C4", 'W'}, {"C2", "N3", 'W'}, {"C4", "C5", 'H'}, {"C5", "C6", 'H'},
        {"N1", "C2", 'S'}, {"C6", "N1", 'S'}};
    return purine ? pur : pyr;
}

const std::array<const char*, 9> kPurineRingTie = {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"};
const std::array<const char*, 6> kPyrRingTie = {"N1", "C2", "N3", "C4", "C5", "C6"};
const std::array<const char*, 9> kPurineRingCenter = {"C4", "N3", "C2", "N1", "C6", "C5", "N7", "C8", "N9"};
const std::array<const char*, 6> kPyrRingCenter = {"C2", "N1", "C6", "C5", "C4", "N3"};

// ---- small geometry / lookup helpers ---------------------------------------

const Vector3d* atom_coords(const Residue& r, const std::string& name) {
    const core::Atom* a = r.get_atom(name);
    return a ? &a->coords : nullptr;
}

bool is_base_hbond_atom(const std::string& a) { return kBaseHbondAtoms.count(a) > 0; }

// Glycosidic-N coordinate (typing-driven), or nullptr.
const Vector3d* gly_coords(const Residue& r, const BaseTyping& typing) {
    const std::string name = typing.glycosidic_n_name(r.base_type());
    if (name.empty()) return nullptr;
    return atom_coords(r, name);
}

Vector3d centroid(const std::vector<Vector3d>& pts) {
    Vector3d s;
    for (const auto& p : pts) s = s + p;
    return s / static_cast<double>(pts.size());
}

// planar projection of v perpendicular to unit-ish z
Vector3d planar(const Vector3d& v, const Vector3d& z) { return v - z * v.dot(z); }

double dotget(const std::unordered_map<char, double>& d, char f, double def) {
    const auto it = d.find(f);
    return it == d.end() ? def : it->second;
}

// ---- H-bond predicate helpers ----------------------------------------------

bool n7_in_hbonds(const std::string& res_id, const std::vector<HBond>& hbonds) {
    for (const auto& hb : hbonds) {
        if (!is_base_hbond_atom(hb.donor_atom) || !is_base_hbond_atom(hb.acceptor_atom)) continue;
        if (hb.donor_res_id == res_id && hb.donor_atom == "N7") return true;
        if (hb.acceptor_res_id == res_id && hb.acceptor_atom == "N7") return true;
    }
    return false;
}

bool n3_is_only_base_hbond(const std::string& res_id, const std::vector<HBond>& hbonds) {
    bool has_n3 = false;
    int other = 0;
    for (const auto& hb : hbonds) {
        for (auto [rid, atom] : {std::pair{std::cref(hb.donor_res_id), std::cref(hb.donor_atom)},
                                 std::pair{std::cref(hb.acceptor_res_id), std::cref(hb.acceptor_atom)}}) {
            if (rid.get() != res_id || !is_base_hbond_atom(atom.get())) continue;
            if (atom.get() == "N3") has_n3 = true;
            else ++other;
        }
    }
    return has_n3 && other == 0;
}

bool has_boundary_atom_in_hbonds(const Residue& res, const std::vector<HBond>& hbonds,
                                 const BaseTyping& typing) {
    const auto boundary = boundary_atoms(res.base_type(), typing);
    if (boundary.empty()) return false;
    for (const auto& hb : hbonds) {
        if (!is_base_hbond_atom(hb.donor_atom) || !is_base_hbond_atom(hb.acceptor_atom)) continue;
        if (hb.donor_res_id == res.res_id() && boundary.count(hb.donor_atom)) return true;
        if (hb.acceptor_res_id == res.res_id() && boundary.count(hb.acceptor_atom)) return true;
    }
    return false;
}

bool pyr_O2_in_hbonds(const std::string& res_id, const std::vector<HBond>& hbonds) {
    for (const auto& hb : hbonds) {
        if (!is_base_hbond_atom(hb.donor_atom) || !is_base_hbond_atom(hb.acceptor_atom)) continue;
        if (hb.donor_res_id == res_id && hb.donor_atom == "O2") return true;
        if (hb.acceptor_res_id == res_id && hb.acceptor_atom == "O2") return true;
    }
    return false;
}

// Face from H-bonds: if all unambiguous atoms agree, else 0.
char unambiguous_face_from_hbonds(const Residue& res, const std::vector<HBond>& hbonds,
                                  const BaseTyping& typing) {
    if (hbonds.empty()) return 0;
    const auto unamb = unambiguous_face(res.base_type(), typing);
    if (unamb.empty()) return 0;
    std::set<char> seen;
    for (const auto& hb : hbonds) {
        if (!is_base_hbond_atom(hb.donor_atom) || !is_base_hbond_atom(hb.acceptor_atom)) continue;
        if (hb.donor_res_id == res.res_id()) {
            auto it = unamb.find(hb.donor_atom);
            if (it != unamb.end()) seen.insert(it->second);
        }
        if (hb.acceptor_res_id == res.res_id()) {
            auto it = unamb.find(hb.acceptor_atom);
            if (it != unamb.end()) seen.insert(it->second);
        }
    }
    return seen.size() == 1 ? *seen.begin() : 0;
}

// Boundary competing faces: sorted unique competing faces (size>=2) or empty.
std::vector<char> boundary_competing_faces(const Residue& res, const std::vector<HBond>& hbonds,
                                           const BaseTyping& typing) {
    const auto boundary = boundary_atoms(res.base_type(), typing);
    if (boundary.empty()) return {};
    std::set<char> competing;
    for (const auto& hb : hbonds) {
        if (!is_base_hbond_atom(hb.donor_atom) || !is_base_hbond_atom(hb.acceptor_atom)) continue;
        if (hb.donor_res_id == res.res_id()) {
            auto it = boundary.find(hb.donor_atom);
            if (it != boundary.end()) competing.insert(it->second.begin(), it->second.end());
        }
        if (hb.acceptor_res_id == res.res_id()) {
            auto it = boundary.find(hb.acceptor_atom);
            if (it != boundary.end()) competing.insert(it->second.begin(), it->second.end());
        }
    }
    if (competing.size() >= 2) return {competing.begin(), competing.end()};  // set => sorted
    return {};
}

bool is_rna_pyrimidine(const std::string& base, const BaseTyping& typing) {
    if (kDnaBaseCodes.count(upper(base))) return false;
    return !typing.is_purine(base);
}

// ---- Layer 1: centroid face detection --------------------------------------

struct DetectResult {
    char face = 'W';
    std::unordered_map<char, double> dots;
};

DetectResult detect_face_with_dots(const ReferenceFrame& frame, const Residue& res,
                                   const Vector3d& partner_gly, const BaseTyping& typing) {
    DetectResult out;
    const Vector3d* gly = gly_coords(res, typing);
    if (gly == nullptr) return out;  // ("W", {})
    const Vector3d anchor = *gly;
    const Vector3d z = frame.z_axis();
    const Vector3d v_planar = planar(partner_gly - anchor, z);
    const double v_norm = v_planar.norm();
    if (v_norm < 1e-6) return out;
    const Vector3d v_unit = v_planar / v_norm;

    char best_face = 'W';
    double best_dot = -std::numeric_limits<double>::infinity();
    for (const auto& [face, atom_names] : face_atoms_for(res.base_type(), typing)) {
        std::vector<Vector3d> coords;
        for (const auto& n : atom_names)
            if (const Vector3d* c = atom_coords(res, n)) coords.push_back(*c);
        if (coords.empty()) continue;
        const Vector3d fv = centroid(coords) - anchor;
        const Vector3d fv_planar = planar(fv, z);
        const double fv_norm = fv_planar.norm();
        if (fv_norm < 1e-6) continue;
        const double dot = v_unit.dot(fv_planar / fv_norm);
        out.dots[face] = dot;
        if (dot > best_dot) {
            best_dot = dot;
            best_face = face;
        }
    }
    // Pyrimidine W/S bias.
    if (!typing.is_purine(res.base_type()) && best_face == 'S') {
        if (out.dots.count('W') && out.dots.count('S') &&
            out.dots['S'] - out.dots['W'] < kPyrWsGap)
            best_face = 'W';
    }
    out.face = best_face;
    return out;
}

// Face by ring-bond proximity -> (face or 0, margin)
std::pair<char, double> face_by_ring_bond_proximity(const Residue& res, const Vector3d& partner,
                                                    const BaseTyping& typing) {
    const bool purine = typing.is_purine(res.base_type());
    std::unordered_map<char, double> best_per_face;
    for (const auto& [a1, a2, face] : ring_bonds(purine)) {
        const Vector3d* c1 = atom_coords(res, a1);
        const Vector3d* c2 = atom_coords(res, a2);
        if (!c1 || !c2) continue;
        const Vector3d mid = (*c1 + *c2) / 2.0;
        const double dist = (partner - mid).norm();
        auto it = best_per_face.find(face);
        if (it == best_per_face.end() || dist < it->second) best_per_face[face] = dist;
    }
    if (best_per_face.empty()) return {0, std::numeric_limits<double>::infinity()};
    std::vector<std::pair<char, double>> ranked(best_per_face.begin(), best_per_face.end());
    std::sort(ranked.begin(), ranked.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    const double margin =
        ranked.size() >= 2 ? ranked[1].second - ranked[0].second : std::numeric_limits<double>::infinity();
    return {ranked[0].first, margin};
}

// Layer 2: confirm the face by atoms
char confirm_face_by_atoms(const Residue& res, char face_guess, const ReferenceFrame& frame,
                           const Vector3d& partner_gly, const std::unordered_map<char, double>& dots,
                           const std::vector<HBond>& hbonds, const BaseTyping& typing) {
    const bool have_hb = !hbonds.empty();
    if (have_hb) {
        char unamb = unambiguous_face_from_hbonds(res, hbonds, typing);
        if (unamb != 0) return unamb;
    }
    const Vector3d* gly = gly_coords(res, typing);
    const Vector3d anchor = gly ? *gly : frame.origin;
    const Vector3d z = frame.z_axis();
    const Vector3d partner_planar = planar(partner_gly - anchor, z);
    const double p_norm = partner_planar.norm();
    if (p_norm < 1e-6) return face_guess;
    const Vector3d p_unit = partner_planar / p_norm;

    if (atom_coords(res, "N7")) {
        if (have_hb && n7_in_hbonds(res.res_id(), hbonds)) return 'H';
        const Vector3d n7v = planar(*atom_coords(res, "N7") - anchor, z);
        const double n7_norm = n7v.norm();
        if (n7_norm > 0.5 && p_unit.dot(n7v / n7_norm) >= kCosThreshold) {
            if (!dots.empty()) {
                const double h_dot = dotget(dots, 'H', -1.0);
                const double w_dot = dotget(dots, 'W', -1.0);
                if (h_dot >= w_dot - kN7GuardMargin) return 'H';
            } else {
                return 'H';
            }
        }
    }
    if (atom_coords(res, "N7") && have_hb && n3_is_only_base_hbond(res.res_id(), hbonds))
        return 'S';

    if (atom_coords(res, "O2'")) {
        const Vector3d o2pv = planar(*atom_coords(res, "O2'") - anchor, z);
        const double o2p_norm = o2pv.norm();
        if (o2p_norm > 0.5 && p_unit.dot(o2pv / o2p_norm) >= kCosThreshold) return 'S';
    } else if (face_guess == 'S') {
        return 'W';
    }

    bool trigger = false;
    char hint = 0;
    if (have_hb && has_boundary_atom_in_hbonds(res, hbonds, typing)) {
        trigger = true;
    } else if (dots.size() >= 2) {
        std::vector<double> vals;
        for (const auto& [f, d] : dots) vals.push_back(d);
        std::sort(vals.rbegin(), vals.rend());
        if (vals[0] - vals[1] < kCentroidAmbiguityGap) trigger = true;
    }
    if (!trigger && have_hb) {
        hint = unambiguous_face_from_hbonds(res, hbonds, typing);
        if (hint != 0 && hint != face_guess && (hint == 'W' || hint == 'S') &&
            (face_guess == 'W' || face_guess == 'S'))
            trigger = true;
    }
    if (trigger) {
        auto [ring_face, ring_margin] = face_by_ring_bond_proximity(res, partner_gly, typing);
        if (ring_face != 0) {
            if (ring_margin < kRingBondMarginThreshold && have_hb) {
                char unamb = hint != 0 ? hint : unambiguous_face_from_hbonds(res, hbonds, typing);
                if (unamb != 0 && unamb != ring_face && (unamb == 'W' || unamb == 'S') &&
                    (ring_face == 'W' || ring_face == 'S'))
                    return unamb;
            }
            return ring_face;
        }
    }
    return face_guess;
}

// Resolve an ambiguous face from the partner (cross-check)
char resolve_ambiguous_face_by_partner(const Residue& amb, const std::vector<char>& competing,
                                       const Residue& partner, char partner_face,
                                       const BaseTyping& typing) {
    const bool partner_pur = typing.is_purine(partner.base_type());
    std::vector<Vector3d> pmids;
    for (const auto& [a1, a2, face] : ring_bonds(partner_pur)) {
        if (face != partner_face) continue;
        const Vector3d* c1 = atom_coords(partner, a1);
        const Vector3d* c2 = atom_coords(partner, a2);
        if (!c1 || !c2) continue;
        pmids.push_back((*c1 + *c2) / 2.0);
    }
    if (pmids.empty()) return 0;
    const Vector3d pref = centroid(pmids);

    const bool amb_pur = typing.is_purine(amb.base_type());
    std::unordered_map<char, double> avg;
    for (char face : competing) {
        std::vector<double> distances;
        for (const auto& [a1, a2, bface] : ring_bonds(amb_pur)) {
            if (bface != face) continue;
            const Vector3d* c1 = atom_coords(amb, a1);
            const Vector3d* c2 = atom_coords(amb, a2);
            if (!c1 || !c2) continue;
            distances.push_back((pref - (*c1 + *c2) / 2.0).norm());
        }
        if (distances.size() >= 2) {
            std::sort(distances.begin(), distances.end());
            avg[face] = (distances[0] + distances[1]) / 2.0;
        } else if (!distances.empty()) {
            avg[face] = distances[0];
        }
    }
    if (avg.size() < 2) return 0;
    char best = 0;
    double bestd = std::numeric_limits<double>::infinity();
    // Ties resolve to the first key in insertion order; iterate
    // competing in its sorted order for determinism.
    for (char face : competing) {
        auto it = avg.find(face);
        if (it != avg.end() && it->second < bestd) {
            bestd = it->second;
            best = face;
        }
    }
    return best;
}

void cross_check_faces(const Residue& res1, char& face1, const Residue& res2, char& face2,
                       const std::vector<HBond>& hbonds, const BaseTyping& typing) {
    char known1 = unambiguous_face_from_hbonds(res1, hbonds, typing);
    char known2 = unambiguous_face_from_hbonds(res2, hbonds, typing);
    if (known1 == 0 && known2 != 0) {
        auto competing = boundary_competing_faces(res1, hbonds, typing);
        if (!competing.empty()) {
            char r = resolve_ambiguous_face_by_partner(res1, competing, res2, known2, typing);
            if (r != 0) face1 = r;
        }
    }
    if (known2 == 0 && known1 != 0) {
        auto competing = boundary_competing_faces(res2, hbonds, typing);
        if (!competing.empty()) {
            char r = resolve_ambiguous_face_by_partner(res2, competing, res1, known1, typing);
            if (r != 0) face2 = r;
        }
    }
}

// H-bond face signal -> (primary, alt or 0)
std::pair<char, char> hbond_face_signal(const Residue& res, const std::vector<HBond>& hbonds,
                                        char geometric_guess, const std::unordered_map<char, double>& dots,
                                        const BaseTyping& typing) {
    char unamb = unambiguous_face_from_hbonds(res, hbonds, typing);
    if (unamb != 0) return {unamb, 0};
    auto competing = boundary_competing_faces(res, hbonds, typing);
    if (competing.empty()) return {geometric_guess, 0};

    if (competing.size() == 2) {
        char primary;
        if (!dots.empty()) {
            primary = dotget(dots, competing[0], 0.0) >= dotget(dots, competing[1], 0.0)
                          ? competing[0]
                          : competing[1];
            if (is_rna_pyrimidine(res.base_type(), typing) && primary == 'S' &&
                (competing[0] == 'W' || competing[1] == 'W') && dots.count('S') && dots.count('W') &&
                dots.at('S') - dots.at('W') < kPyrWsGap)
                primary = 'W';
        } else {
            bool guess_in = geometric_guess == competing[0] || geometric_guess == competing[1];
            if (guess_in) primary = geometric_guess;
            else if (competing[0] == 'W' || competing[1] == 'W') primary = 'W';
            else primary = std::min(competing[0], competing[1]);
        }
        char alt = (competing[0] == primary) ? competing[1] : competing[0];
        return {primary, alt};
    }
    if (competing.size() >= 3) {
        bool has_w = std::find(competing.begin(), competing.end(), 'W') != competing.end();
        return {has_w ? 'W' : geometric_guess, 0};
    }
    return {geometric_guess, 0};
}

// Layer 3: refine the face by close contacts
char refine_face_by_close_contacts(const Residue& res, const Residue& partner, char current_face,
                                   const BaseTyping& typing) {
    if (current_face != 'W' && current_face != 'H') return current_face;
    const auto unamb = unambiguous_face(res.base_type(), typing);
    if (unamb.empty()) return current_face;
    std::vector<Vector3d> partner_pos;
    for (const auto& n : kBaseRingAtoms)
        if (const Vector3d* c = atom_coords(partner, n)) partner_pos.push_back(*c);
    if (partner_pos.empty()) return current_face;

    int W = 0, H = 0;
    for (const auto& [atom_name, face] : unamb) {
        if (face != 'W' && face != 'H') continue;
        const Vector3d* ap = atom_coords(res, atom_name);
        if (!ap) continue;
        double min_d = std::numeric_limits<double>::infinity();
        for (const auto& p : partner_pos) min_d = std::min(min_d, (*ap - p).norm());
        if (min_d < kVotingContactDistance) (face == 'W' ? W : H)++;
    }
    if (current_face == 'W' && H > 0 && W == 0) return 'H';
    if (current_face == 'H' && W > 0 && H == 0) return 'W';
    return current_face;
}

std::optional<Vector3d> ring_center_for_tiebreaker(const Residue& res, const BaseTyping& typing) {
    std::vector<Vector3d> coords;
    if (typing.is_purine(res.base_type())) {
        for (const char* a : kPurineRingTie)
            if (const Vector3d* c = atom_coords(res, a)) coords.push_back(*c);
    } else {
        for (const char* a : kPyrRingTie)
            if (const Vector3d* c = atom_coords(res, a)) coords.push_back(*c);
    }
    if (coords.empty()) return std::nullopt;
    return centroid(coords);
}

// Layer 4: pyrimidine W/S tiebreaker
char apply_pyr_WS_tiebreaker(const Residue& res, const Residue& partner, char current_face,
                             const std::vector<HBond>& hbonds, const BaseTyping& typing) {
    if (current_face != 'W' && current_face != 'S') return current_face;
    if (!is_rna_pyrimidine(res.base_type(), typing)) return current_face;
    if (!pyr_O2_in_hbonds(res.res_id(), hbonds)) return current_face;
    auto pc = ring_center_for_tiebreaker(partner, typing);
    const Vector3d* n3 = atom_coords(res, "N3");
    const Vector3d* n1 = atom_coords(res, "N1");
    if (!pc || !n3 || !n1) return current_face;
    const double d_n3 = (*pc - *n3).norm();
    const double d_n1 = (*pc - *n1).norm();
    if (d_n3 < d_n1 - kPyrWsTiebreakerMarginA) return 'W';
    if (d_n1 < d_n3 - kPyrWsTiebreakerMarginA) return 'S';
    return current_face;
}

// ---- orientation -----------------------------------------------------------

std::optional<Vector3d> ring_center(const Residue& res, const BaseTyping& typing) {
    std::vector<Vector3d> coords;
    if (typing.is_purine(res.base_type())) {
        for (const char* a : kPurineRingCenter)
            if (const Vector3d* c = atom_coords(res, a)) coords.push_back(*c);
    } else {
        for (const char* a : kPyrRingCenter)
            if (const Vector3d* c = atom_coords(res, a)) coords.push_back(*c);
    }
    if (coords.size() < 3) return std::nullopt;
    return centroid(coords);
}

char orientation_by_yaxis(const ReferenceFrame& f1, const ReferenceFrame& f2) {
    return f1.y_axis().dot(f2.y_axis()) > 0 ? 'c' : 't';
}

char classify_orientation(const ReferenceFrame& f1, const ReferenceFrame& f2, const Residue& res1,
                          const Residue& res2, const BaseTyping& typing) {
    const Vector3d* c1p1 = atom_coords(res1, "C1'");
    const Vector3d* c1p2 = atom_coords(res2, "C1'");
    const Vector3d* g1n = gly_coords(res1, typing);
    const Vector3d* g2n = gly_coords(res2, typing);
    auto rc1 = ring_center(res1, typing);
    auto rc2 = ring_center(res2, typing);
    if (!c1p1 || !c1p2 || !g1n || !g2n || !rc1 || !rc2) return orientation_by_yaxis(f1, f2);

    const Vector3d g1 = *g1n - *c1p1;
    const Vector3d g2 = *g2n - *c1p2;
    const Vector3d r12 = *rc2 - *rc1;
    const double r12_len = r12.norm();
    if (r12_len < 0.1) return 'c';
    const Vector3d r12u = r12 / r12_len;
    const Vector3d g1p = g1 - r12u * g1.dot(r12u);
    const Vector3d g2p = g2 - r12u * g2.dot(r12u);
    const double n1 = g1p.norm(), n2 = g2p.norm();
    if (n1 < 1e-6 || n2 < 1e-6) return orientation_by_yaxis(f1, f2);
    return (g1p / n1).dot(g2p / n2) > 0 ? 'c' : 't';
}

// Canonical LW class.
LwResult canonicalize_lw(char orientation, char face1, char face2) {
    std::string o(1, orientation);
    if (face_priority(face1) <= face_priority(face2))
        return {o + face1 + face2, false};
    return {o + face2 + face1, true};
}

// ---- classification by geometry --------------------------------------------

struct GeoResult {
    std::string lw_class;
    bool swapped = false;
    bool uncertain = false;
    std::string alt_class;
    bool alt_swapped = false;
};

GeoResult classify_by_geometry(const ReferenceFrame& f1, const ReferenceFrame& f2,
                               const Residue& res1, const Residue& res2, char orientation,
                               const std::vector<HBond>& hbonds, const BaseTyping& typing) {
    const Vector3d* g2 = gly_coords(res2, typing);
    const Vector3d* g1 = gly_coords(res1, typing);
    const Vector3d partner1 = g2 ? *g2 : f2.origin;
    const Vector3d partner2 = g1 ? *g1 : f1.origin;

    DetectResult d1 = detect_face_with_dots(f1, res1, partner1, typing);
    char face1 = confirm_face_by_atoms(res1, d1.face, f1, partner1, d1.dots, hbonds, typing);
    DetectResult d2 = detect_face_with_dots(f2, res2, partner2, typing);
    char face2 = confirm_face_by_atoms(res2, d2.face, f2, partner2, d2.dots, hbonds, typing);

    const bool have_hb = !hbonds.empty();
    if (have_hb) cross_check_faces(res1, face1, res2, face2, hbonds, typing);

    char alt_face1 = 0, alt_face2 = 0;
    if (have_hb) {
        auto [p1, a1] = hbond_face_signal(res1, hbonds, face1, d1.dots, typing);
        face1 = p1;
        alt_face1 = a1;
        auto [p2, a2] = hbond_face_signal(res2, hbonds, face2, d2.dots, typing);
        face2 = p2;
        alt_face2 = a2;
        if ((alt_face1 == 0) != (alt_face2 == 0)) {
            alt_face1 = 0;
            alt_face2 = 0;
        } else if (alt_face1 != 0 && alt_face2 != 0) {
            const double m1 = !d1.dots.empty()
                                  ? std::abs(dotget(d1.dots, face1, 0.0) - dotget(d1.dots, alt_face1, 0.0))
                                  : 1.0;
            const double m2 = !d2.dots.empty()
                                  ? std::abs(dotget(d2.dots, face2, 0.0) - dotget(d2.dots, alt_face2, 0.0))
                                  : 1.0;
            if (m1 >= kBothBoundaryAmbigMargin && m2 >= kBothBoundaryAmbigMargin) {
                alt_face1 = 0;
                alt_face2 = 0;
            }
        }
    }
    bool uncertain = alt_face1 != 0 || alt_face2 != 0;

    face1 = refine_face_by_close_contacts(res1, res2, face1, typing);
    face2 = refine_face_by_close_contacts(res2, res1, face2, typing);

    if (have_hb) {
        const Vector3d z1 = f1.z_axis(), z2 = f2.z_axis();
        const double z1n = z1.norm(), z2n = z2.norm();
        if (z1n > 1e-9 && z2n > 1e-9) {
            double cos = std::abs((z1 / z1n).dot(z2 / z2n));
            cos = std::max(-1.0, std::min(1.0, cos));
            const double plane_deg = std::acos(cos) * 180.0 / kPi;
            if (plane_deg > kPyrWsTiebreakerPlaneDeg) {
                face1 = apply_pyr_WS_tiebreaker(res1, res2, face1, hbonds, typing);
                face2 = apply_pyr_WS_tiebreaker(res2, res1, face2, hbonds, typing);
            }
        }
    }

    LwResult canon = canonicalize_lw(orientation, face1, face2);
    GeoResult out;
    out.lw_class = canon.lw_class;
    out.swapped = canon.swapped;
    out.uncertain = uncertain;

    if (uncertain) {
        char af1 = alt_face1 != 0 ? alt_face1 : face1;
        char af2 = alt_face2 != 0 ? alt_face2 : face2;
        LwResult canon_alt = canonicalize_lw(orientation, af1, af2);
        if (canon_alt.lw_class == out.lw_class) {
            // Both-flip lands on the SAME canonical class. When BOTH residues are
            // independently boundary-uncertain (we only reach here when a
            // centroid margin is below the ambiguity threshold), the genuine
            // ambiguity is a SINGLE-face flip (e.g. cWH vs cWW/cHH), which the
            // both-flip alt cannot express. Flip only the less-confident residue.
            if (alt_face1 != 0 && alt_face2 != 0) {
                const double m1 = std::abs(dotget(d1.dots, face1, 0.0) -
                                           dotget(d1.dots, alt_face1, 0.0));
                const double m2 = std::abs(dotget(d2.dots, face2, 0.0) -
                                           dotget(d2.dots, alt_face2, 0.0));
                LwResult single = (m1 <= m2)
                    ? canonicalize_lw(orientation, alt_face1, face2)
                    : canonicalize_lw(orientation, face1, alt_face2);
                if (single.lw_class != out.lw_class) {
                    out.alt_class = single.lw_class;
                    out.alt_swapped = single.swapped;
                } else {
                    out.uncertain = false;
                }
            } else {
                out.uncertain = false;
            }
        } else if (canon_alt.swapped != out.swapped &&
                   (out.lw_class.find('S') != std::string::npos ||
                    canon_alt.lw_class.find('S') != std::string::npos)) {
            // For a W/S boundary the swapped-flag mismatch is a canonical-sort
            // artifact and that path isn't reliably separable, so commit to the
            // geometric primary. A *pure W/H* boundary ambiguity (tWW vs tWH,
            // cWW vs cHH) is a genuine class difference - keep it rather than
            // silently committing.
            out.uncertain = false;
        } else {
            out.alt_class = canon_alt.lw_class;
            out.alt_swapped = canon_alt.swapped;
        }
    }
    return out;
}

}  // namespace

LwResult classify_lw_class(const ReferenceFrame& frame1, const ReferenceFrame& frame2,
                           const Residue& res1, const Residue& res2,
                           const std::vector<HBond>& hbonds, const BaseTyping& typing) {
    const char orientation = classify_orientation(frame1, frame2, res1, res2, typing);
    GeoResult geo =
        classify_by_geometry(frame1, frame2, res1, res2, orientation, hbonds, typing);
    if (geo.uncertain) return {geo.lw_class + "|" + geo.alt_class, geo.swapped};
    return {geo.lw_class, geo.swapped};
}

}  // namespace pairfinder::algorithms::classification
