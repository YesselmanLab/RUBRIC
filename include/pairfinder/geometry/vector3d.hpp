/**
 * @file vector3d.hpp
 * @brief Minimal 3D vector value type for coordinate geometry.
 *
 * Header-only, constexpr-friendly.
 */
#ifndef PAIRFINDER_GEOMETRY_VECTOR3D_HPP
#define PAIRFINDER_GEOMETRY_VECTOR3D_HPP

#include <array>
#include <cmath>

namespace pairfinder::geometry {

/// A 3D Cartesian vector (Angstroms, in the conventions of the PDB).
struct Vector3d {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vector3d() = default;
    constexpr Vector3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vector3d operator+(const Vector3d& o) const { return {x + o.x, y + o.y, z + o.z}; }
    constexpr Vector3d operator-(const Vector3d& o) const { return {x - o.x, y - o.y, z - o.z}; }
    constexpr Vector3d operator*(double s) const { return {x * s, y * s, z * s}; }
    constexpr Vector3d operator/(double s) const { return {x / s, y / s, z / s}; }

    constexpr double dot(const Vector3d& o) const { return x * o.x + y * o.y + z * o.z; }

    constexpr Vector3d cross(const Vector3d& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }

    double norm() const { return std::sqrt(dot(*this)); }

    Vector3d normalized() const {
        const double n = norm();
        return n > 0.0 ? *this / n : Vector3d{};
    }

    double distance_to(const Vector3d& o) const { return (*this - o).norm(); }
};

}  // namespace pairfinder::geometry

#endif  // PAIRFINDER_GEOMETRY_VECTOR3D_HPP
