/**
 * @file finder.cpp
 * @brief Slot-based H-bond finder.
 */
#include <pairfinder/algorithms/hbond/finder.hpp>

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <numeric>
#include <unordered_set>

namespace pairfinder::algorithms::hbond {

namespace {

using ResMap = std::unordered_map<std::string, const core::Residue*>;

// Base atoms used to reject intra-residue base-to-base contacts.
const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};

const std::unordered_set<std::string> kAminoAcids = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"};

bool is_intra_base_pair(const std::string& res1, const std::string& res2,
                        const std::string& atom1, const std::string& atom2) {
    if (res1 != res2) return false;
    return kBaseAtoms.count(atom1) > 0 && kBaseAtoms.count(atom2) > 0;
}

// Split "CHAIN-RESNAME-NUM" -> (chain, num). Returns false if malformed.
bool split_res_id(const std::string& res_id, std::string& chain, long& num) {
    const auto p2 = res_id.rfind('-');
    if (p2 == std::string::npos || p2 == 0) return false;
    const auto p1 = res_id.rfind('-', p2 - 1);
    if (p1 == std::string::npos) return false;
    chain = res_id.substr(0, p1);
    try {
        num = std::stol(res_id.substr(p2 + 1));
    } catch (const std::exception&) {
        return false;
    }
    return true;
}

// Peptide-bond geometry (consecutive protein N...O) - not a real H-bond.
bool is_consecutive_backbone_bond(const HBondCandidate& c, const ResMap& residues) {
    if (c.donor_atom != "N" || c.acceptor_atom != "O") return false;
    const auto d_it = residues.find(c.donor_res_id);
    const auto a_it = residues.find(c.acceptor_res_id);
    if (d_it == residues.end() || a_it == residues.end()) return false;
    if (kAminoAcids.count(d_it->second->base_type()) == 0) return false;
    if (kAminoAcids.count(a_it->second->base_type()) == 0) return false;
    std::string c1, c2;
    long n1 = 0, n2 = 0;
    if (!split_res_id(c.donor_res_id, c1, n1)) return false;
    if (!split_res_id(c.acceptor_res_id, c2, n2)) return false;
    return c1 == c2 && std::abs(n1 - n2) == 1;
}

}  // namespace

// ---------------------------------------------------------------------------
//  ResidueSlotCache
// ---------------------------------------------------------------------------

ResidueSlotCache::ResidueSlotCache(std::string res_id, std::string base_type,
                                   AtomMap atoms, const HBondChemistry& chem)
    : res_id_(std::move(res_id)),
      base_type_(std::move(base_type)),
      atoms_(std::move(atoms)),
      chem_(&chem) {}

const std::vector<HbSlot>& ResidueSlotCache::get_h_slots(const std::string& atom_name) {
    const auto it = h_slots_.find(atom_name);
    if (it != h_slots_.end()) return it->second;
    if (!base_normal_) base_normal_ = compute_base_normal(atoms_);
    auto [ins, _] = h_slots_.emplace(
        atom_name, chem_->predict_h_slots(base_type_, atom_name, atoms_, *base_normal_));
    return ins->second;
}

const std::vector<HbSlot>& ResidueSlotCache::get_lp_slots(const std::string& atom_name) {
    const auto it = lp_slots_.find(atom_name);
    if (it != lp_slots_.end()) return it->second;
    if (!base_normal_) base_normal_ = compute_base_normal(atoms_);
    auto [ins, _] = lp_slots_.emplace(
        atom_name, chem_->predict_lp_slots(base_type_, atom_name, atoms_, *base_normal_));
    return ins->second;
}

// ---------------------------------------------------------------------------
//  HBondFinder
// ---------------------------------------------------------------------------

HBondFinder::HBondFinder(const HBondChemistry& chemistry, double max_distance,
                         double min_alignment, double short_distance_threshold)
    : chem_(&chemistry),
      max_distance_(max_distance),
      min_alignment_(min_alignment),
      short_distance_threshold_(short_distance_threshold) {}

ResidueSlotCache& HBondFinder::get_cache(const core::Residue& residue) {
    const auto it = residue_cache_.find(residue.res_id());
    if (it != residue_cache_.end()) return it->second;
    AtomMap atoms;
    for (const auto& atom : residue.atoms()) atoms[atom.name] = atom.coords;
    auto [ins, _] = residue_cache_.emplace(
        residue.res_id(),
        ResidueSlotCache(residue.res_id(), residue.base_type(), std::move(atoms), *chem_));
    return ins->second;
}

void HBondFinder::add_directional_candidates(const core::Residue& donor,
                                             const core::Residue& acceptor,
                                             std::vector<HBondCandidate>& out) {
    // Which acceptor atoms can accept, looked up once rather than once per donor atom.
    const auto& a_atoms = acceptor.atoms();
    std::vector<char> can_accept(a_atoms.size());
    for (std::size_t k = 0; k < a_atoms.size(); ++k)
        can_accept[k] = chem_->acceptor_capacity(acceptor.base_type(), a_atoms[k].name).has_value();
    for (const auto& d_atom : donor.atoms()) {
        if (!chem_->donor_capacity(donor.base_type(), d_atom.name)) continue;
        for (std::size_t k = 0; k < a_atoms.size(); ++k) {
            if (!can_accept[k]) continue;
            const auto& a_atom = a_atoms[k];
            if (is_intra_base_pair(donor.res_id(), acceptor.res_id(), d_atom.name,
                                   a_atom.name))
                continue;
            const double dist = d_atom.coords.distance_to(a_atom.coords);
            if (dist > max_distance_) continue;
            out.push_back({donor.res_id(), acceptor.res_id(), d_atom.name, a_atom.name,
                           dist, d_atom.coords, a_atom.coords, -1, -1, 0.0});
        }
    }
}

std::vector<HBondCandidate> HBondFinder::find_candidates(const core::Residue& res1,
                                                         const core::Residue& res2) {
    std::vector<HBondCandidate> candidates;
    add_directional_candidates(res1, res2, candidates);
    add_directional_candidates(res2, res1, candidates);
    return candidates;
}

void HBondFinder::compute_alignments(std::vector<HBondCandidate>& candidates,
                                     const ResidueMap& residues) {
    for (auto& c : candidates) {
        const auto d_it = residues.find(c.donor_res_id);
        const auto a_it = residues.find(c.acceptor_res_id);
        if (d_it == residues.end() || a_it == residues.end()) continue;
        ResidueSlotCache& donor_cache = get_cache(*d_it->second);
        const std::vector<HbSlot>& h_slots = donor_cache.get_h_slots(c.donor_atom);
        ResidueSlotCache& acceptor_cache = get_cache(*a_it->second);
        const std::vector<HbSlot>& lp_slots = acceptor_cache.get_lp_slots(c.acceptor_atom);
        if (h_slots.empty() || lp_slots.empty()) continue;

        const Vector3d d2a = hb_normalize(c.acceptor_pos - c.donor_pos);
        const Vector3d a2d = d2a * -1.0;
        double best_h = -2.0;
        int best_h_idx = 0;
        for (int i = 0; i < static_cast<int>(h_slots.size()); ++i) {
            const double align = h_slots[i].direction.dot(d2a);
            if (align > best_h) { best_h = align; best_h_idx = i; }
        }
        double best_lp = -2.0;
        int best_lp_idx = 0;
        for (int i = 0; i < static_cast<int>(lp_slots.size()); ++i) {
            const double align = lp_slots[i].direction.dot(a2d);
            if (align > best_lp) { best_lp = align; best_lp_idx = i; }
        }
        c.h_slot_idx = best_h_idx;
        c.lp_slot_idx = best_lp_idx;
        c.alignment_score = best_h + best_lp;
    }
}

std::vector<int> HBondFinder::find_best_selection(const std::vector<HBondCandidate>& valid,
                                                  const ResidueMap& /*residues*/) {
    const int n = static_cast<int>(valid.size());
    if (n == 0) return {};

    // Conflict graph: share donor-slot or acceptor-slot key.
    std::vector<std::vector<int>> conflicts(n);
    auto shares = [&](int i, int j) {
        const auto& a = valid[i];
        const auto& b = valid[j];
        const bool donor_same = a.donor_res_id == b.donor_res_id &&
                                a.donor_atom == b.donor_atom && a.h_slot_idx == b.h_slot_idx;
        const bool acc_same = a.acceptor_res_id == b.acceptor_res_id &&
                              a.acceptor_atom == b.acceptor_atom &&
                              a.lp_slot_idx == b.lp_slot_idx;
        return donor_same || acc_same;
    };
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (shares(i, j)) { conflicts[i].push_back(j); conflicts[j].push_back(i); }

    // Greedy for large sets (>50).
    if (n > 50) {
        std::vector<int> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            if (valid[a].alignment_score != valid[b].alignment_score)
                return valid[a].alignment_score > valid[b].alignment_score;
            return valid[a].distance < valid[b].distance;
        });
        std::vector<char> excluded(n, 0);
        std::vector<int> selected;
        for (int idx : order) {
            if (excluded[idx]) continue;
            selected.push_back(idx);
            excluded[idx] = 1;
            for (int j : conflicts[idx]) excluded[j] = 1;
        }
        return selected;
    }

    // Backtracking maximum-weight (count, then min total distance) independent set.
    std::vector<int> best_selection;
    long best_count = 0;
    double best_dist = std::numeric_limits<double>::infinity();
    std::vector<int> current;
    std::vector<char> excluded(n, 0);

    std::function<void(int)> search = [&](int idx) {
        const int remaining = n - idx;
        if (static_cast<int>(current.size()) + remaining < best_count) return;
        if (idx == n) {
            double total = 0.0;
            for (int i : current) total += valid[i].distance;
            const long count = static_cast<long>(current.size());
            if (count > best_count || (count == best_count && total < best_dist)) {
                best_selection = current;
                best_count = count;
                best_dist = total;
            }
            return;
        }
        // Option 1: skip candidate idx.
        search(idx + 1);
        // Option 2: include candidate idx if not excluded.
        if (!excluded[idx]) {
            std::vector<int> newly;
            for (int j : conflicts[idx]) {
                if (!excluded[j]) { excluded[j] = 1; newly.push_back(j); }
            }
            current.push_back(idx);
            search(idx + 1);
            current.pop_back();
            for (int j : newly) excluded[j] = 0;
        }
    };
    search(0);
    return best_selection;
}

std::vector<HBond> HBondFinder::select_optimal(std::vector<HBondCandidate>& candidates,
                                               const ResidueMap& residues) {
    if (candidates.empty()) return {};
    compute_alignments(candidates, residues);

    std::vector<HBondCandidate> valid;
    for (const auto& c : candidates) {
        if (c.h_slot_idx < 0 || c.lp_slot_idx < 0) continue;
        if (is_consecutive_backbone_bond(c, residues)) continue;
        if (c.distance < short_distance_threshold_ || c.alignment_score >= min_alignment_)
            valid.push_back(c);
    }

    const std::vector<int> best = find_best_selection(valid, residues);
    std::vector<HBond> selected;
    selected.reserve(best.size());
    for (int idx : best) {
        const auto& c = valid[idx];
        selected.push_back({c.donor_res_id, c.acceptor_res_id, c.donor_atom, c.acceptor_atom,
                            c.distance, c.h_slot_idx, c.lp_slot_idx, c.alignment_score});
    }
    return selected;
}

std::vector<HBond> HBondFinder::find_between(const core::Residue& res1,
                                            const core::Residue& res2) {
    std::vector<HBondCandidate> candidates = find_candidates(res1, res2);
    ResidueMap residues;
    residues[res1.res_id()] = &res1;
    residues[res2.res_id()] = &res2;
    return select_optimal(candidates, residues);
}

std::vector<HBond> HBondFinder::find_face_voters(const core::Residue& res1,
                                                 const core::Residue& res2,
                                                 double distance_threshold) {
    std::vector<HBond> out = find_between(res1, res2);
    std::unordered_set<std::string> seen;
    auto key = [](const std::string& dr, const std::string& da, const std::string& ar,
                  const std::string& aa) { return dr + '\t' + da + '\t' + ar + '\t' + aa; };
    for (const auto& hb : out)
        seen.insert(key(hb.donor_res_id, hb.donor_atom, hb.acceptor_res_id, hb.acceptor_atom));

    std::vector<HBondCandidate> dropped = find_candidates(res1, res2);
    std::vector<HBondCandidate> rev = find_candidates(res2, res1);
    dropped.insert(dropped.end(), rev.begin(), rev.end());
    for (const auto& c : dropped) {
        if (c.distance > distance_threshold) continue;
        const std::string k = key(c.donor_res_id, c.donor_atom, c.acceptor_res_id, c.acceptor_atom);
        if (seen.count(k)) continue;
        out.push_back({c.donor_res_id, c.acceptor_res_id, c.donor_atom, c.acceptor_atom,
                       c.distance, 0, 0, c.alignment_score});
        seen.insert(k);
    }
    return out;
}

}  // namespace pairfinder::algorithms::hbond
