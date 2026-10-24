/**
 * @file clashes.cpp
 * @brief Heavy-atom steric clash detection.
 */
#include <pairfinder/scoring/clashes.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string_view>
#include <tuple>
#include <unordered_set>

namespace pairfinder::scoring {

bool is_base_atom(std::string_view atom_name);  // defined below

namespace {

// Richardson / MolProbity heavy-atom van der Waals radii (A).
constexpr double kRadiusC = 1.70;
constexpr double kRadiusN = 1.55;
constexpr double kRadiusO = 1.52;
constexpr double kRadiusP = 1.80;
constexpr double kRadiusS = 1.80;
constexpr double kRadiusDefault = 1.20;  ///< metals / unknown

/// Largest possible min_ok (P + P). Beyond this no pair can overlap, so the
/// inner loop rejects on squared distance without ever calling sqrt.
constexpr double kMaxMinOk = kRadiusP + kRadiusP;
constexpr double kMaxMinOkSq = kMaxMinOk * kMaxMinOk;

/// Element inferred from an atom name. Returns '\0' for hydrogen/deuterium.
/// Residue context is not needed here: only nucleotide residues reach this
/// code, so there is no C-alpha / calcium ambiguity to resolve.
char element_of(std::string_view atom_name) {
    std::size_t i = 0;
    while (i < atom_name.size() && atom_name[i] >= '0' && atom_name[i] <= '9') ++i;
    if (i >= atom_name.size()) return 'C';
    const char c = static_cast<char>(std::toupper(static_cast<unsigned char>(atom_name[i])));
    if (c == 'H' || c == 'D') return '\0';  // hydrogens are ignored entirely
    return c;
}

bool is_polar(char element) { return element == 'N' || element == 'O' || element == 'S'; }

double radius_of(char element) {
    switch (element) {
        case 'C': return kRadiusC;
        case 'N': return kRadiusN;
        case 'O': return kRadiusO;
        case 'P': return kRadiusP;
        case 'S': return kRadiusS;
        default:  return kRadiusDefault;
    }
}

/// Phosphate-group atom names. Sugar atoms are identified by their prime, so
/// these are the only non-base atoms that need naming.
bool is_phosphate_name(std::string_view n) {
    return n == "P" || n == "OP1" || n == "OP2" || n == "OP3" ||
           n == "O1P" || n == "O2P" || n == "O3P";
}

/// Compact per-atom record: everything the inner loop needs, nothing it does not.
struct HeavyAtom {
    geometry::Vector3d pos;
    double radius;
    bool donor;                ///< can donate an H-bond (HBondChemistry)
    bool acceptor;             ///< can accept an H-bond
    bool polar;                ///< N / O / S
    bool base;                 ///< base atom (not sugar / phosphate)
    const std::string* name;
};

/// Donor / acceptor capacity for one atom, resolved once at collection time.
/// Keyed on the residue's base type, exactly as HBondFinder does it.
std::pair<bool, bool> hbond_capacity(const algorithms::hbond::HBondChemistry& chem,
                                     const std::string& base_type, const std::string& atom) {
    const auto d = chem.donor_capacity(base_type, atom);
    const auto a = chem.acceptor_capacity(base_type, atom);
    return {d.has_value() && *d > 0, a.has_value() && *a > 0};
}

void collect_heavy(const core::Residue& residue, std::vector<HeavyAtom>& out,
                   const algorithms::hbond::HBondChemistry& chem) {
    out.clear();
    out.reserve(residue.atoms().size());
    const std::string& base_type = residue.base_type();
    for (const auto& atom : residue.atoms()) {
        const char e = element_of(atom.name);
        if (e == '\0') continue;
        const auto [don, acc] = hbond_capacity(chem, base_type, atom.name);
        out.push_back(HeavyAtom{atom.coords, radius_of(e), don, acc, is_polar(e),
                                is_base_atom(atom.name), &atom.name});
    }
}

/// The H-bond floor applies only when the pair could actually hydrogen-bond --
/// one side donating, the other accepting. Two acceptors (two carbonyls, say)
/// have no tautomer that bonds, so they are held to the full van der Waals sum
/// and a close approach registers as the clash it is.
///
/// An atom absent from the capacity tables (the phosphodiester O3' / O5', which
/// HBondChemistry does not list either way) is treated as UNKNOWN, not as
/// incapable: 2'-OH...O5' is a real H-bond, and charging it the van der Waals
/// sum would turn a gap in the H-bond tables into a false clash. Two polar atoms
/// keep the lenient floor whenever either side's chemistry is unknown.
bool hbond_plausible(const HeavyAtom& a, const HeavyAtom& b) {
    if ((a.donor && b.acceptor) || (a.acceptor && b.donor)) return true;
    const bool a_unknown = !a.donor && !a.acceptor;
    const bool b_unknown = !b.donor && !b.acceptor;
    return (a_unknown || b_unknown) && a.polar && b.polar;
}

double min_acceptable(const HeavyAtom& a, const HeavyAtom& b) {
    return hbond_plausible(a, b) ? kHBondFloorMinOk : (a.radius + b.radius);
}

/// True when the two residues are joined by an O3'(i)-P(j) phosphodiester bond
/// (either direction). Two O(1) lookups - the fast path for the vast majority
/// of pairs, which are not covalently linked.
bool covalently_linked(const core::Residue& r1, const core::Residue& r2) {
    const auto bonded = [](const core::Atom* o3, const core::Atom* p) {
        return o3 && p && o3->coords.distance_to(p->coords) < kBondMax;
    };
    return bonded(r1.get_atom("O3'"), r2.get_atom("P")) ||
           bonded(r2.get_atom("O3'"), r1.get_atom("P"));
}

/// Cross-residue atom-index pairs within kBondDepth covalent bonds.
///
/// Only called for the rare covalently linked pair. Distance alone cannot tell
/// a covalent bond from a badly mis-modelled contact (an O3'-P bond at 1.59 A
/// and an impossible N...O at 1.60 A look identical), so a cross-residue bond
/// is admitted only for the O3'-P linkage; everything else that close stays
/// visible as the error it is.
/// Keyed by atom names ("O3'|P") rather than indices, so the same set works for
/// both the per-pair path and the grid path.
std::string name_key(const std::string& a, const std::string& b) { return a + "|" + b; }

std::unordered_set<std::string> bonded_cross_pairs(const std::vector<HeavyAtom>& a,
                                                   const std::vector<HeavyAtom>& b) {
    const std::size_t na = a.size(), nb = b.size(), n = na + nb;
    const auto pos = [&](std::size_t i) -> const geometry::Vector3d& {
        return i < na ? a[i].pos : b[i - na].pos;
    };
    const auto name = [&](std::size_t i) -> const std::string& {
        return i < na ? *a[i].name : *b[i - na].name;
    };

    // Bond graph in compressed form (edge list -> per-node offsets), built with a
    // few flat allocations rather than one vector per atom. Edges: any two atoms closer than kBondMax on one side, and across the two residues
    // only the O3'-P phosphodiester.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (pos(i).distance_to(pos(j)) >= kBondMax) continue;
            if ((i < na) != (j < na)) {
                const std::string& ni = name(i);
                const std::string& nj = name(j);
                const bool phosphodiester =
                    (ni == "O3'" && nj == "P") || (ni == "P" && nj == "O3'");
                if (!phosphodiester) continue;
            }
            edges.emplace_back(static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j));
        }
    }
    std::vector<std::uint32_t> offset(n + 1, 0), fill(n, 0), adj(2 * edges.size());
    for (const auto& [i, j] : edges) { ++offset[i + 1]; ++offset[j + 1]; }
    for (std::size_t i = 0; i < n; ++i) offset[i + 1] += offset[i];
    // Neighbours keep the order the per-atom vectors had (edges are generated in
    // (i, j) order), so the breadth-first search visits nodes in the same order.
    for (const auto& [i, j] : edges) {
        adj[offset[i] + fill[i]++] = j;
        adj[offset[j] + fill[j]++] = i;
    }

    std::unordered_set<std::string> excluded;
    std::vector<int> depth(n);
    std::vector<std::uint32_t> queue;
    queue.reserve(n);
    for (std::size_t start = 0; start < na; ++start) {
        std::fill(depth.begin(), depth.end(), -1);
        depth[start] = 0;
        queue.clear();
        queue.push_back(static_cast<std::uint32_t>(start));
        for (std::size_t head = 0; head < queue.size(); ++head) {
            const std::uint32_t node = queue[head];
            if (depth[node] >= kBondDepth) continue;
            for (std::uint32_t e = offset[node]; e < offset[node + 1]; ++e) {
                const std::uint32_t nb_idx = adj[e];
                if (depth[nb_idx] != -1) continue;
                depth[nb_idx] = depth[node] + 1;
                if (nb_idx >= na) excluded.insert(name_key(name(start), name(nb_idx)));
                queue.push_back(nb_idx);
            }
        }
    }
    return excluded;
}

/// Worst contact between two prepared atom lists. Returns (overlap, i, j).
struct BestContact {
    double overlap = -std::numeric_limits<double>::infinity();
    std::size_t i = 0, j = 0;
    bool valid = false;
};

BestContact worst_contact(const std::vector<HeavyAtom>& a, const std::vector<HeavyAtom>& b,
                          const std::unordered_set<std::string>* excluded,
                          bool base_involved_only) {
    BestContact best;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const auto& ai = a[i];
        for (std::size_t j = 0; j < b.size(); ++j) {
            const auto& bj = b[j];
            // Scoring rule: only contacts implicating a base are charged to
            // the pair. Backbone-only contacts go to pair_backbone_clash.
            if (base_involved_only && !ai.base && !bj.base) continue;
            const double dx = ai.pos.x - bj.pos.x;
            const double dy = ai.pos.y - bj.pos.y;
            const double dz = ai.pos.z - bj.pos.z;
            const double d2 = dx * dx + dy * dy + dz * dz;
            // A pair can only beat the incumbent if it is closer than
            // (min_ok - best.overlap); reject on squared distance alone, so the
            // overwhelming majority cost no sqrt and no set lookup. The bound is
            // infinite until a candidate exists, so EVERY pair is still eligible
            // to set the first best -- otherwise a residue pair whose closest
            // base contact lies beyond kMaxMinOk would return whichever pair
            // happened to be examined first rather than the true maximum.
            const double bound =
                best.valid ? min_acceptable(ai, bj) - best.overlap
                           : std::numeric_limits<double>::infinity();
            if (bound <= 0.0) continue;
            if (std::isfinite(bound) && d2 >= bound * bound) continue;
            if (excluded != nullptr && excluded->count(name_key(*ai.name, *bj.name)) != 0) {
                continue;
            }
            const double overlap = min_acceptable(ai, bj) - std::sqrt(d2);
            if (!best.valid || overlap > best.overlap) {
                best = BestContact{overlap, i, j, true};
            }
        }
    }
    return best;
}

}  // namespace

/// severity = clamp((overlap - kClashThreshold) / (kClashSeverityCap - kClashThreshold), 0, 1)
///
/// Zero at (and below) the clash threshold, rising to 1 at the cap. Rebasing on
/// the threshold means a contact that has only just become a clash costs almost
/// nothing.
/// Linear ramp from kClashThreshold to kClashSeverityCap.
double clash_severity(double overlap) {
    if (std::isnan(overlap)) return 0.0;
    constexpr double kSpan = kClashSeverityCap - kClashThreshold;
    return std::min(1.0, std::max(0.0, (overlap - kClashThreshold) / kSpan));
}

double clash_severity(std::optional<double> overlap) {
    return overlap.has_value() ? clash_severity(*overlap) : 0.0;
}

bool is_rna_residue_name(std::string_view res_name, const core::BaseTyping& typing) {
    const std::string parent = typing.normalize(res_name);
    if (parent.size() != 1) return false;
    switch (parent[0]) {
        case 'A': case 'G': case 'C': case 'U': case 'T': case 'I': return true;
        default: return false;
    }
}

bool is_rna_res_id(std::string_view res_id, const core::BaseTyping& typing) {
    const std::size_t first = res_id.find('-');
    if (first == std::string_view::npos) return false;
    const std::size_t second = res_id.find('-', first + 1);
    if (second == std::string_view::npos) return false;
    return is_rna_residue_name(res_id.substr(first + 1, second - first - 1), typing);
}

bool is_base_atom(std::string_view atom_name) {
    return atom_name.find('\'') == std::string_view::npos && !is_phosphate_name(atom_name);
}

std::optional<double> worst_overlap(const core::Residue& res1, const core::Residue& res2,
                                    const algorithms::hbond::HBondChemistry& chem) {
    std::vector<HeavyAtom> a, b;
    collect_heavy(res1, a, chem);
    collect_heavy(res2, b, chem);
    if (a.empty() || b.empty()) return std::nullopt;

    std::unordered_set<std::string> excluded;
    const bool linked = covalently_linked(res1, res2);
    if (linked) excluded = bonded_cross_pairs(a, b);

    const BestContact best = worst_contact(a, b, linked ? &excluded : nullptr, true);
    if (!best.valid || !std::isfinite(best.overlap)) return std::nullopt;
    return best.overlap;
}

std::optional<ClashContact> pair_clash(const core::Residue& res1, const core::Residue& res2,
                                       const algorithms::hbond::HBondChemistry& chem,
                                       double threshold) {
    std::vector<HeavyAtom> a, b;
    collect_heavy(res1, a, chem);
    collect_heavy(res2, b, chem);
    if (a.empty() || b.empty()) return std::nullopt;

    std::unordered_set<std::string> excluded;
    const bool linked = covalently_linked(res1, res2);
    if (linked) excluded = bonded_cross_pairs(a, b);

    const BestContact best = worst_contact(a, b, linked ? &excluded : nullptr, true);
    if (!best.valid || !std::isfinite(best.overlap) || best.overlap <= threshold) {
        return std::nullopt;
    }
    ClashContact contact;
    contact.res_id1 = res1.res_id();
    contact.atom1 = *a[best.i].name;
    contact.res_id2 = res2.res_id();
    contact.atom2 = *b[best.j].name;
    contact.distance = a[best.i].pos.distance_to(b[best.j].pos);
    contact.overlap = best.overlap;
    return contact;
}

std::optional<ClashContact> pair_backbone_clash(const core::Residue& res1,
                                                const core::Residue& res2,
                                                const algorithms::hbond::HBondChemistry& chem,
                                                double threshold) {
    std::vector<HeavyAtom> a, b;
    collect_heavy(res1, a, chem);
    collect_heavy(res2, b, chem);
    if (a.empty() || b.empty()) return std::nullopt;

    std::unordered_set<std::string> excluded;
    const bool linked = covalently_linked(res1, res2);
    if (linked) excluded = bonded_cross_pairs(a, b);

    // Unfiltered worst contact; keep it only when NO base atom is involved,
    // otherwise pair_clash already scored it.
    const BestContact best = worst_contact(a, b, linked ? &excluded : nullptr, false);
    if (!best.valid || !std::isfinite(best.overlap) || best.overlap <= threshold) {
        return std::nullopt;
    }
    if (a[best.i].base || b[best.j].base) return std::nullopt;

    ClashContact contact;
    contact.res_id1 = res1.res_id();
    contact.atom1 = *a[best.i].name;
    contact.res_id2 = res2.res_id();
    contact.atom2 = *b[best.j].name;
    contact.distance = a[best.i].pos.distance_to(b[best.j].pos);
    contact.overlap = best.overlap;
    return contact;
}

// ---------------------------------------------------------------------------
// Whole-structure detection: uniform grid, O(N) rather than O(N^2).
// ---------------------------------------------------------------------------

namespace {

struct GridAtom {
    geometry::Vector3d pos;
    double radius;
    bool donor;
    bool acceptor;
    bool polar;
    const std::string* name;
    const std::string* res_id;
    std::uint32_t residue;   ///< index into the residue list, for the same-residue test
};

struct CellKey {
    std::int64_t x, y, z;
    bool operator==(const CellKey& o) const { return x == o.x && y == o.y && z == o.z; }
};

struct CellHash {
    std::size_t operator()(const CellKey& k) const {
        // Cheap spatial hash; collisions only cost a few extra comparisons.
        return static_cast<std::size_t>(k.x * 73856093LL ^ k.y * 19349663LL ^
                                        k.z * 83492791LL);
    }
};

}  // namespace

std::vector<ClashContact> structure_clashes(
    const std::unordered_map<std::string, const core::Residue*>& residues,
    const core::BaseTyping& typing, const algorithms::hbond::HBondChemistry& chem,
    double threshold) {
    // Flatten every RNA heavy atom into one contiguous array.
    std::vector<GridAtom> atoms;
    std::vector<const core::Residue*> residue_ptrs;
    std::vector<std::uint32_t> residue_begin;   // residue k's atoms: [begin[k], begin[k + 1])
    atoms.reserve(residues.size() * 20);
    for (const auto& [res_id, residue_ptr] : residues) {
        if (!is_rna_res_id(res_id, typing)) continue;
        const core::Residue& residue = *residue_ptr;
        const auto idx = static_cast<std::uint32_t>(residue_ptrs.size());
        residue_ptrs.push_back(&residue);
        residue_begin.push_back(static_cast<std::uint32_t>(atoms.size()));
        for (const auto& atom : residue.atoms()) {
            const char e = element_of(atom.name);
            if (e == '\0') continue;
            const auto [don, acc] = hbond_capacity(chem, residue.base_type(), atom.name);
            atoms.push_back(GridAtom{atom.coords, radius_of(e), don, acc, is_polar(e),
                                     &atom.name, &residue.res_id(), idx});
        }
    }
    residue_begin.push_back(static_cast<std::uint32_t>(atoms.size()));
    if (atoms.size() < 2) return {};

    // A residue's heavy atoms as bonded_cross_pairs takes them: the same atoms in
    // the same order collect_heavy gives (both skip hydrogens only), taken from the
    // flattened array instead of being re-collected. Only position and name are read.
    const auto heavy_of = [&](std::uint32_t res, std::vector<HeavyAtom>& out) {
        out.clear();
        for (std::uint32_t k = residue_begin[res]; k < residue_begin[res + 1]; ++k)
            out.push_back(HeavyAtom{atoms[k].pos, atoms[k].radius, atoms[k].donor,
                                    atoms[k].acceptor, atoms[k].polar, false, atoms[k].name});
    };

    // Bucket atoms into cells of side kMaxMinOk so every possible contact lies
    // in one of the 27 neighbouring cells.
    const double cell = kMaxMinOk;
    std::unordered_map<CellKey, std::vector<std::uint32_t>, CellHash> grid;
    grid.reserve(atoms.size());
    const auto key_of = [cell](const geometry::Vector3d& p) {
        return CellKey{static_cast<std::int64_t>(std::floor(p.x / cell)),
                       static_cast<std::int64_t>(std::floor(p.y / cell)),
                       static_cast<std::int64_t>(std::floor(p.z / cell))};
    };
    for (std::uint32_t i = 0; i < atoms.size(); ++i) grid[key_of(atoms[i].pos)].push_back(i);

    // Linked-residue exclusions are computed lazily and memoised: a structure
    // has ~N backbone links, and each costs one small BFS. One entry per residue
    // pair holds whether it is linked and, if so, its exclusion set.
    std::unordered_map<std::uint64_t, std::unordered_set<std::string>> link_cache;
    std::unordered_map<std::uint64_t, const std::unordered_set<std::string>*> linked_cache;
    std::vector<HeavyAtom> ha, hb;

    std::vector<ClashContact> clashes;
    // Walk the grid a cell at a time: its (up to 27) neighbouring cells are looked
    // up once for all of its atoms rather than once per atom. Every atom pair is
    // still visited exactly once (j > i); only the order of discovery changes, and
    // the result is sorted below.
    std::string link_key;   // reused buffer for name_key, so a lookup allocates nothing
    std::vector<const std::vector<std::uint32_t>*> near;
    near.reserve(27);
    for (const auto& [cell_key, members] : grid) {
        near.clear();
        for (std::int64_t dx = -1; dx <= 1; ++dx)
            for (std::int64_t dy = -1; dy <= 1; ++dy)
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it =
                        grid.find(CellKey{cell_key.x + dx, cell_key.y + dy, cell_key.z + dz});
                    if (it != grid.end()) near.push_back(&it->second);
                }
        for (const std::uint32_t i : members) {
            const GridAtom& ai = atoms[i];
            for (const auto* neighbours : near) {
                for (const std::uint32_t j : *neighbours) {
                    if (j <= i) continue;                      // each pair once
                    const GridAtom& aj = atoms[j];
                    if (ai.residue == aj.residue) continue;    // same residue
                    const double ddx = ai.pos.x - aj.pos.x;
                    const double ddy = ai.pos.y - aj.pos.y;
                    const double ddz = ai.pos.z - aj.pos.z;
                    const double d2 = ddx * ddx + ddy * ddy + ddz * ddz;
                    if (d2 >= kMaxMinOkSq) continue;
                    const bool ai_unknown = !ai.donor && !ai.acceptor;
                    const bool aj_unknown = !aj.donor && !aj.acceptor;
                    const bool can_hb =
                        ((ai.donor && aj.acceptor) || (ai.acceptor && aj.donor)) ||
                        ((ai_unknown || aj_unknown) && ai.polar && aj.polar);
                    const double min_ok =
                        can_hb ? kHBondFloorMinOk : (ai.radius + aj.radius);
                    const double d = std::sqrt(d2);
                    const double overlap = min_ok - d;
                    if (overlap <= threshold) continue;

                    // Only now, for a genuine candidate, pay for the
                    // covalent-linkage check.
                    const std::uint64_t rk =
                        (static_cast<std::uint64_t>(std::min(ai.residue, aj.residue)) << 32) |
                        static_cast<std::uint64_t>(std::max(ai.residue, aj.residue));
                    const std::uint32_t lo_res = std::min(ai.residue, aj.residue);
                    const std::uint32_t hi_res = std::max(ai.residue, aj.residue);
                    auto lit = linked_cache.find(rk);
                    if (lit == linked_cache.end()) {
                        const std::unordered_set<std::string>* excl = nullptr;
                        if (covalently_linked(*residue_ptrs[lo_res], *residue_ptrs[hi_res])) {
                            heavy_of(lo_res, ha);
                            heavy_of(hi_res, hb);
                            excl = &(link_cache[rk] = bonded_cross_pairs(ha, hb));
                        }
                        lit = linked_cache.emplace(rk, excl).first;
                    }
                    if (lit->second != nullptr) {
                        // bonded_cross_pairs keys run (low residue -> high residue).
                        const std::string& lo_name = (ai.residue == lo_res) ? *ai.name
                                                                            : *aj.name;
                        const std::string& hi_name = (ai.residue == lo_res) ? *aj.name
                                                                            : *ai.name;
                        link_key.assign(lo_name).append(1, '|').append(hi_name);
                        if (lit->second->count(link_key) != 0) continue;
                    }

                    ClashContact c;
                    c.res_id1 = *ai.res_id;
                    c.atom1 = *ai.name;
                    c.res_id2 = *aj.res_id;
                    c.atom2 = *aj.name;
                    c.distance = d;
                    c.overlap = overlap;
                    clashes.push_back(std::move(c));
                }
            }
        }
    }
    // Worst first; equal overlaps (symmetric copies in an assembly) are ordered by
    // name, so the list does not depend on the order the grid happens to be walked.
    std::sort(clashes.begin(), clashes.end(), [](const ClashContact& a, const ClashContact& b) {
        if (a.overlap != b.overlap) return a.overlap > b.overlap;
        return std::tie(a.res_id1, a.atom1, a.res_id2, a.atom2) <
               std::tie(b.res_id1, b.atom1, b.res_id2, b.atom2);
    });
    return clashes;
}

}  // namespace pairfinder::scoring
