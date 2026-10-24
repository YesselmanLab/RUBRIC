/**
 * @file candidate_finder.cpp
 * @brief Pair-candidate finding + geometric validation.
 */
#include <pairfinder/algorithms/candidate_finder.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <optional>
#include <string>
#include <unordered_map>

#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms {

namespace {

constexpr double kPi = 3.14159265358979323846;
using geometry::Vector3d;

double clamp_unit(double x) { return x > 1.0 ? 1.0 : (x < -1.0 ? -1.0 : x); }

// Single-letter res_name from res_id: strip a leading 'D' from 2-char
// D-prefixed DNA codes; keep everything else.
std::string extract_res_name(const std::string& res_id) {
    const auto p1 = res_id.find('-');
    if (p1 == std::string::npos) return res_id;
    const auto p2 = res_id.find('-', p1 + 1);
    const std::string name =
        res_id.substr(p1 + 1, (p2 == std::string::npos ? res_id.size() : p2) - p1 - 1);
    if (name.size() == 2 && name[0] == 'D') return name.substr(1);
    return name;
}

// Same-chain adjacent residues (|i-j| == 1) cannot form a real base pair - they
// are backbone-connected neighbours (stacked), not paired. Excluded here at
// candidate generation so they never become candidates. (Identical parsing to
// selection::same_chain_adjacent, which applies the same rule later in the
// pipeline.)
bool same_chain_adjacent(const std::string& a, const std::string& b) {
    auto split = [](const std::string& r, std::string& chain, long& num) {
        const auto p2 = r.rfind('-');
        if (p2 == std::string::npos || p2 == 0) return false;
        const auto p1 = r.rfind('-', p2 - 1);
        if (p1 == std::string::npos) return false;
        chain = r.substr(0, p1);
        try { num = std::stol(r.substr(p2 + 1)); } catch (...) { return false; }
        return true;
    };
    std::string c1, c2;
    long n1 = 0, n2 = 0;
    if (!split(a, c1, n1) || !split(b, c2, n2)) return false;
    return c1 == c2 && std::abs(n1 - n2) == 1;
}

// --- ring-centroid stacking signature -------------------------------------
// A candidate is "ring-stacked" when the vector between the two base RING
// centroids is offset more ALONG the mean base normal (vertical) than in the
// base plane (lateral) - i.e. one base sits over the other. Unlike the
// dNN-based is_stacking_signature, this is class-independent and works for the
// sugar-edge classes (cWS/cSS/tWS) where real pairs have low dNN by design.
constexpr double kRingStackVertMin = 2.8;

std::optional<Vector3d> ring_centroid(const core::Residue& res) {
    static const std::array<const char*, 9> purine = {"N9", "C8", "N7", "C5", "C6",
                                                      "N1", "C2", "N3", "C4"};
    static const std::array<const char*, 6> pyrim = {"N1", "C2", "N3", "C4", "C5", "C6"};
    Vector3d sum{};
    int n = 0;
    auto add = [&](const char* nm) {
        if (const core::Atom* a = res.get_atom(nm)) { sum = sum + a->coords; ++n; }
    };
    if (core::is_purine(res.base_type()))
        for (const char* nm : purine) add(nm);
    else
        for (const char* nm : pyrim) add(nm);
    if (n < 3) return std::nullopt;
    return sum / static_cast<double>(n);
}

bool is_ring_stacked(const Vector3d& c1, const Vector3d& c2,
                     const core::ReferenceFrame& f1, const core::ReferenceFrame& f2) {
    Vector3d zave = f1.z_axis().dot(f2.z_axis()) > 0 ? f1.z_axis() + f2.z_axis()
                                                     : f2.z_axis() - f1.z_axis();
    const double zl = zave.norm();
    if (zl < 1e-10) return false;
    zave = zave / zl;
    const Vector3d vec = c2 - c1;
    const double vert = std::abs(vec.dot(zave));
    const double lat = std::sqrt(std::max(0.0, vec.dot(vec) - vert * vert));
    return vert > lat && vert > kRingStackVertMin;
}

}  // namespace

ValidationResult validate_pair(const core::ReferenceFrame& f1,
                               const core::ReferenceFrame& f2,
                               const Vector3d& n1n9_1, const Vector3d& n1n9_2,
                               const ValidationThresholds& thr) {
    ValidationResult r;
    const Vector3d dorg_vec = f1.origin - f2.origin;
    r.dorg = dorg_vec.norm();

    r.dir_x = f1.x_axis().dot(f2.x_axis());
    r.dir_y = f1.y_axis().dot(f2.y_axis());
    r.dir_z = f1.z_axis().dot(f2.z_axis());

    // Mean helix axis: add z-axes if aligned, subtract if anti.
    Vector3d zave = r.dir_z > 0 ? f1.z_axis() + f2.z_axis() : f2.z_axis() - f1.z_axis();
    const double zlen = zave.norm();
    zave = zlen > 1e-10 ? zave / zlen : f1.z_axis();
    r.d_v = std::abs(dorg_vec.dot(zave));

    r.plane_angle = std::acos(std::abs(clamp_unit(r.dir_z))) * 180.0 / kPi;
    r.dNN = (n1n9_1 - n1n9_2).norm();
    r.quality_score =
        r.dorg + thr.d_v_weight * r.d_v + r.plane_angle / thr.plane_angle_divisor;

    r.distance_check = r.dorg <= thr.max_dorg;
    r.d_v_check = r.d_v <= thr.max_d_v;
    r.plane_angle_check = r.plane_angle <= thr.max_plane_angle;
    r.dNN_check = r.dNN >= thr.min_dNN;
    r.is_valid =
        r.distance_check && r.d_v_check && r.plane_angle_check && r.dNN_check;
    return r;
}

std::vector<PairCandidate> find_candidates(
    const core::Structure& structure,
    const std::vector<std::pair<std::string, core::ReferenceFrame>>& frames,
    const core::BaseTyping& typing, double max_distance,
    const ValidationThresholds& thr) {
    // res_id -> Residue* for glycosidic-N lookup.
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    // Resolve glycosidic-N coordinate once per framed residue (nullptr-safe).
    struct Node {
        const std::string* res_id;
        const core::ReferenceFrame* frame;
        Vector3d gly;
        bool has_gly;
        Vector3d ring;      ///< base ring centroid (for the stacking signature)
        bool has_ring;
    };
    std::vector<Node> nodes;
    nodes.reserve(frames.size());
    for (const auto& [res_id, frame] : frames) {
        Node n{&res_id, &frame, {}, false, {}, false};
        const auto it = res_by_id.find(res_id);
        if (it != res_by_id.end()) {
            const std::string name = typing.glycosidic_n_name(it->second->base_type());
            if (!name.empty()) {
                if (const core::Atom* a = it->second->get_atom(name)) {
                    n.gly = a->coords;
                    n.has_gly = true;
                }
            }
            if (const auto rc = ring_centroid(*it->second)) {
                n.ring = *rc;
                n.has_ring = true;
            }
        }
        nodes.push_back(n);
    }

    // --- spatial grid (cell = max_distance) ----------------------------------
    // Replaces the O(n^2) all-pairs scan: bin node origins into cells of side
    // max_distance, so all neighbors within max_distance lie in the 27 cells
    // around a node. For each i we gather j>i in those cells, sort ascending,
    // and emit in (i asc, j asc) order.
    const int n = static_cast<int>(nodes.size());
    const double cell = max_distance > 0 ? max_distance : 1.0;
    auto cell_of = [cell](double v) { return static_cast<long>(std::floor(v / cell)); };
    // Collision-free packed cell key: offset each axis index into [0, 2^21) and
    // pack into 63 bits. Avoids the false duplicates an XOR hash would produce
    // when two of the 27 query cells alias to the same bucket.
    constexpr long kOff = 1L << 20;
    auto cell_key = [](long ix, long iy, long iz) {
        return ((ix + kOff) << 42) | ((iy + kOff) << 21) | (iz + kOff);
    };
    std::unordered_map<long, std::vector<int>> grid;
    grid.reserve(n * 2);
    for (int i = 0; i < n; ++i) {
        const auto& o = nodes[i].frame->origin;
        grid[cell_key(cell_of(o.x), cell_of(o.y), cell_of(o.z))].push_back(i);
    }

    std::vector<PairCandidate> out;
    std::vector<int> nbrs;
    for (int i = 0; i < n; ++i) {
        if (!nodes[i].has_gly) continue;  // no candidate can include a gly-less i
        const auto& oi = nodes[i].frame->origin;
        const long cx = cell_of(oi.x), cy = cell_of(oi.y), cz = cell_of(oi.z);
        nbrs.clear();
        for (long dx = -1; dx <= 1; ++dx)
            for (long dy = -1; dy <= 1; ++dy)
                for (long dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find(cell_key(cx + dx, cy + dy, cz + dz));
                    if (it == grid.end()) continue;
                    for (int j : it->second) {
                        if (j <= i || !nodes[j].has_gly) continue;
                        // Use norm() (not squared) so the boundary test is exact;
                        // sqrt is cheap here since the grid already restricted j
                        // to nearby cells.
                        if ((oi - nodes[j].frame->origin).norm() > max_distance) continue;
                        nbrs.push_back(j);
                    }
                }
        std::sort(nbrs.begin(), nbrs.end());  // (i asc, j asc) order preserved
        for (int j : nbrs) {
            // Adjacency is a candidate-generation exclusion: same-chain neighbours
            // are backbone-stacked, never a base pair.
            if (same_chain_adjacent(*nodes[i].res_id, *nodes[j].res_id)) continue;
            PairCandidate c;
            c.res_id1 = *nodes[i].res_id;
            c.res_id2 = *nodes[j].res_id;
            c.res_name1 = extract_res_name(c.res_id1);
            c.res_name2 = extract_res_name(c.res_id2);
            c.validation = validate_pair(*nodes[i].frame, *nodes[j].frame, nodes[i].gly,
                                         nodes[j].gly, thr);
            if (nodes[i].has_ring && nodes[j].has_ring)
                c.validation.ring_stacked = is_ring_stacked(
                    nodes[i].ring, nodes[j].ring, *nodes[i].frame, *nodes[j].frame);
            // Stacking is a GEOMETRY exclusion (like d_v/plane/dNN): a ring-stacked
            // pair is not a base pair, so it never becomes a valid candidate. (This
            // is class-independent, unlike selection::is_stacking_signature which is
            // per-class and must run after classification.)
            if (c.validation.ring_stacked) c.validation.is_valid = false;
            out.push_back(std::move(c));
        }
    }
    return out;
}

}  // namespace pairfinder::algorithms
