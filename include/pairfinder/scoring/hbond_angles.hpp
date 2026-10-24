/**
 * @file hbond_angles.hpp
 * @brief H-bond angle geometry.
 *
 * The angle at an atom is the largest (neighbor - atom - partner) angle over
 * the atom's bonded neighbors (connectivity from the chemistry DB, falling back
 * to the single-letter base type). Shared by the scorer and the CLI per-pair
 * H-bond emitter so both compute identical angles.
 */
#ifndef PAIRFINDER_SCORING_HBOND_ANGLES_HPP
#define PAIRFINDER_SCORING_HBOND_ANGLES_HPP

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::scoring {

/// Angle neighbor-atom-partner in degrees (0 if either arm is degenerate).
inline double angle_at_atom(const geometry::Vector3d& neighbor,
                            const geometry::Vector3d& atom,
                            const geometry::Vector3d& partner) {
    const geometry::Vector3d v1 = neighbor - atom, v2 = partner - atom;
    const double n1 = v1.norm(), n2 = v2.norm();
    if (n1 < 1e-10 || n2 < 1e-10) return 0.0;
    double c = v1.dot(v2) / (n1 * n2);
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c) * 180.0 / 3.14159265358979323846;
}

/// H-bond angle at atom_name: the largest neighbor-atom-partner angle over the
/// atom's bonded neighbors; std::nullopt if no connectivity/neighbor is present.
inline std::optional<double> hbond_angle(const core::Residue& res, const std::string& atom_name,
                                         const geometry::Vector3d& atom_pos,
                                         const geometry::Vector3d& partner_pos,
                                         const algorithms::hbond::HBondChemistry& chem) {
    const std::vector<std::string>* nbrs = &chem.connectivity(res.base_type(), atom_name);
    if (nbrs->empty() && !res.base_type().empty())
        nbrs = &chem.connectivity(res.base_type().substr(0, 1), atom_name);
    if (nbrs->empty()) return std::nullopt;
    std::optional<double> best;
    for (const auto& nn : *nbrs) {
        const core::Atom* na = res.get_atom(nn);
        if (!na) continue;
        const double a = angle_at_atom(na->coords, atom_pos, partner_pos);
        if (!best || a > *best) best = a;
    }
    return best;
}

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_HBOND_ANGLES_HPP
