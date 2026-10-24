/**
 * @file superposition.hpp
 * @brief Least-squares rigid superposition (Kabsch / Horn quaternion method).
 *
 * Finds the rotation R and translation t that best map a "from" point set onto
 * a "to" point set (to_i ~= R*from_i + t), minimizing RMSD. Used to fit the
 * standard base ring geometry onto observed ring atoms -> the base frame.
 *
 * Self-contained (no external linear-algebra dependency): builds Horn's 4x4
 * key matrix and extracts the top eigenvector via a symmetric Jacobi solver.
 * Gives the same optimal rotation as SVD-Kabsch.
 */
#ifndef PAIRFINDER_GEOMETRY_SUPERPOSITION_HPP
#define PAIRFINDER_GEOMETRY_SUPERPOSITION_HPP

#include <array>
#include <vector>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::geometry {

struct SuperpositionResult {
    std::array<std::array<double, 3>, 3> rotation{};  ///< row-major; to ~= R*from + t
    Vector3d translation;
    double rms = 0.0;
    bool valid = false;
};

/// Optimal rigid map from `from` onto `to` (equal sizes, >= 3 points).
SuperpositionResult superpose(const std::vector<Vector3d>& from,
                              const std::vector<Vector3d>& to);

}  // namespace pairfinder::geometry

#endif  // PAIRFINDER_GEOMETRY_SUPERPOSITION_HPP
