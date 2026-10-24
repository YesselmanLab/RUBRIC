/**
 * @file rna_chains.cpp
 * @brief Spatial RNA chain connectivity.
 */
#include <pairfinder/algorithms/rna_chains.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace pairfinder::algorithms {

namespace {

constexpr double kBackboneCutoff = 2.75;
using core::Residue;
using geometry::Vector3d;

const Vector3d* coord(const Residue& r, const std::string& a) {
    const core::Atom* at = r.get_atom(a);
    return at ? &at->coords : nullptr;
}

// Uniform spatial hash on a single atom type (cell size = backbone cutoff).
// query() returns residue indices whose indexed atom lies in the 27 cells
// around a point - a superset that the caller distance-filters. Cuts the O(n^2)
// all-pairs O3'<->P scan to ~O(n) (the rna_chains ribosome bottleneck).
struct AtomGrid {
    std::unordered_map<std::int64_t, std::vector<int>> cells;
    static std::int64_t key(int ix, int iy, int iz) {
        return (static_cast<std::int64_t>(ix) * 73856093LL) ^
               (static_cast<std::int64_t>(iy) * 19349663LL) ^
               (static_cast<std::int64_t>(iz) * 83492791LL);
    }
    static int cell_of(double v) { return static_cast<int>(std::floor(v / kBackboneCutoff)); }
    void insert(const Vector3d& p, int idx) {
        cells[key(cell_of(p.x), cell_of(p.y), cell_of(p.z))].push_back(idx);
    }
    template <class F>
    void for_near(const Vector3d& p, F&& f) const {
        const int cx = cell_of(p.x), cy = cell_of(p.y), cz = cell_of(p.z);
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const auto it = cells.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == cells.end()) continue;
                    for (int idx : it->second) f(idx);
                }
    }
};

}  // namespace

RNAChains RNAChains::from_structure(const core::Structure& structure) {
    // RNA residues (have C1' or P), in structure order.
    std::vector<const Residue*> residues;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues())
            if (res.get_atom("C1'") || res.get_atom("P")) residues.push_back(&res);

    const int n = static_cast<int>(residues.size());
    std::vector<const Vector3d*> o3(n), pp(n);  // per-residue O3' and P (or PA)
    AtomGrid grid_p, grid_o3;
    for (int i = 0; i < n; ++i) {
        o3[i] = coord(*residues[i], "O3'");
        const Vector3d* p = coord(*residues[i], "P");
        if (!p) p = coord(*residues[i], "PA");
        pp[i] = p;
        if (p) grid_p.insert(*p, i);
        if (o3[i]) grid_o3.insert(*o3[i], i);
    }

    // 3' extension: residue whose P is within cutoff of `from`.O3' (connected==1).
    // 5' extension: residue whose O3' is within cutoff of `from`.P  (connected==-1).
    // Both pick the lowest residue index (= first in structure order) on ties.
    std::vector<char> assigned(n, 0);
    auto best_3prime = [&](int from) -> int {
        if (!o3[from]) return -1;
        int best = -1;
        grid_p.for_near(*o3[from], [&](int r) {
            if (assigned[r] || r == best) return;
            if (pp[r] && (*pp[r] - *o3[from]).norm() < kBackboneCutoff && (best < 0 || r < best))
                best = r;
        });
        return best;
    };
    auto best_5prime = [&](int from) -> int {
        if (!pp[from]) return -1;
        int best = -1;
        grid_o3.for_near(*pp[from], [&](int r) {
            if (assigned[r]) return;
            if (o3[r] && (*pp[from] - *o3[r]).norm() < kBackboneCutoff && (best < 0 || r < best))
                best = r;
        });
        return best;
    };

    std::vector<std::vector<std::string>> chains;
    for (int start = 0; start < n; ++start) {
        if (assigned[start]) continue;
        std::vector<int> chain = {start};
        assigned[start] = 1;
        bool changed = true;
        while (changed) {
            changed = false;
            const int nxt = best_3prime(chain.back());
            if (nxt >= 0) {
                chain.push_back(nxt);
                assigned[nxt] = 1;
                changed = true;
            }
            const int prv = best_5prime(chain.front());
            if (prv >= 0) {
                chain.insert(chain.begin(), prv);
                assigned[prv] = 1;
                changed = true;
            }
        }
        std::vector<std::string> ids;
        ids.reserve(chain.size());
        for (int idx : chain) ids.push_back(residues[idx]->res_id());
        chains.push_back(std::move(ids));
    }

    RNAChains out;
    for (int ci = 0; ci < static_cast<int>(chains.size()); ++ci) {
        const auto& chain = chains[ci];
        for (int pos = 0; pos < static_cast<int>(chain.size()); ++pos) {
            out.pos_[chain[pos]] = {ci, pos};
            if (pos + 1 < static_cast<int>(chain.size())) out.next_[chain[pos]] = chain[pos + 1];
            if (pos > 0) out.prev_[chain[pos]] = chain[pos - 1];
        }
    }
    return out;
}

const std::string& RNAChains::next_residue(const std::string& res_id) const {
    const auto it = next_.find(res_id);
    return it == next_.end() ? empty_ : it->second;
}

const std::string& RNAChains::prev_residue(const std::string& res_id) const {
    const auto it = prev_.find(res_id);
    return it == prev_.end() ? empty_ : it->second;
}

std::optional<std::pair<int, int>> RNAChains::position(const std::string& res_id) const {
    const auto it = pos_.find(res_id);
    if (it == pos_.end()) return std::nullopt;
    return it->second;
}

}  // namespace pairfinder::algorithms
