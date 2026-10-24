/**
 * @file chemistry.hpp
 * @brief Donor/acceptor capacities, base connectivity, and slot prediction.
 *
 * Per-(base,atom) donor and
 * acceptor capacities + base connectivity (hardcoded for standard nucleotides /
 * amino acids, extended from ligand_hbond_db.json for modified residues), plus
 * the geometric prediction of hydrogen (donor) and lone-pair (acceptor) slot
 * directions used to score H-bond alignment.
 */
#ifndef PAIRFINDER_ALGORITHMS_HBOND_CHEMISTRY_HPP
#define PAIRFINDER_ALGORITHMS_HBOND_CHEMISTRY_HPP

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms::hbond {

using geometry::Vector3d;

/// Normalize: returns v unchanged if |v| < 1e-10 (not zero).
Vector3d hb_normalize(const Vector3d& v);
/// Angle in degrees [0,180] between two vectors (normalizes first).
double angle_between(const Vector3d& a, const Vector3d& b);
/// Rotate v about axis (normalized) by angle_deg via Rodrigues' formula.
Vector3d rotate_vector(const Vector3d& v, const Vector3d& axis, double angle_deg);

/// Unit normal to the base plane from ring atoms C2,C4,C6,N1,N3 (first 3 present).
Vector3d compute_base_normal(const std::unordered_map<std::string, Vector3d>& atoms);

/// A donor H slot or acceptor lone-pair slot. (bond_directions/max_bonds are
/// retained for fidelity with the bifurcation model; the per-pair
/// selection path used by find_between does not mutate them.)
struct HbSlot {
    Vector3d direction;
    int max_bonds = 2;
    std::vector<Vector3d> bond_directions;
};

/// Atom-name -> coordinate view of a residue (raw PDB names).
using AtomMap = std::unordered_map<std::string, Vector3d>;

/// Per-(base,atom) H-bond chemistry tables + slot prediction.
class HBondChemistry {
public:
    /// Builds the standard tables, then extends capacities from the ligand DB
    /// at @p ligand_db_path (silently skipped if the file is unreadable).
    /// @param use_cache when true, load the precomputed ``hbond_chemistry_cache.json``
    ///   (an exact serialization of the donor/acceptor capacity tables) from the
    ///   ligand DB's directory instead of parsing the 1.5 MB DB; falls back to
    ///   parsing when it is missing. Pass false to force a fresh parse.
    explicit HBondChemistry(const std::filesystem::path& ligand_db_path, bool use_cache = true);

    std::optional<int> donor_capacity(const std::string& base, const std::string& atom) const;
    std::optional<int> acceptor_capacity(const std::string& base, const std::string& atom) const;

    /// Serialize the donor/acceptor capacity tables to JSON, so they can be
    /// reloaded without re-parsing the ligand DB.
    void save_cache(const std::filesystem::path& path) const;

    /// Predicted H slots for a donor atom (empty if not a donor / atom missing).
    std::vector<HbSlot> predict_h_slots(const std::string& base, const std::string& atom,
                                        const AtomMap& atoms,
                                        const Vector3d& base_normal) const;
    /// Predicted lone-pair slots for an acceptor atom (empty if not an acceptor).
    std::vector<HbSlot> predict_lp_slots(const std::string& base, const std::string& atom,
                                         const AtomMap& atoms,
                                         const Vector3d& base_normal) const;

    /// Bonded ring neighbors of (base, atom) - used for H-bond angle geometry.
    const std::vector<std::string>& connectivity(const std::string& base,
                                                 const std::string& atom) const;

private:
    static std::string key(const std::string& base, const std::string& atom);

    /// Populate the capacity tables from a cache file; false (tables untouched)
    /// if absent or malformed, so the caller parses the DB instead.
    bool load_cache(const std::filesystem::path& path);

    std::unordered_map<std::string, int> donor_cap_;
    std::unordered_map<std::string, int> acceptor_cap_;
    std::unordered_map<std::string, std::vector<std::string>> connectivity_;
    std::vector<std::string> empty_;  ///< returned for missing connectivity
};

}  // namespace pairfinder::algorithms::hbond

#endif  // PAIRFINDER_ALGORITHMS_HBOND_CHEMISTRY_HPP
