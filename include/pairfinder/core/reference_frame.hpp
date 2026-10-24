/**
 * @file reference_frame.hpp
 * @brief Base reference frame: origin + 3x3 rotation (least-squares fit).
 *
 * Serialised form: origin == "translation", rotation == "rotation_matrix"
 * (row-major), rms_fit.
 */
#ifndef PAIRFINDER_CORE_REFERENCE_FRAME_HPP
#define PAIRFINDER_CORE_REFERENCE_FRAME_HPP

#include <array>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::core {

struct ReferenceFrame {
    geometry::Vector3d origin;                          ///< translation
    std::array<std::array<double, 3>, 3> rotation{};    ///< row-major 3x3
    double rms_fit = 0.0;
    bool valid = false;

    /// Column i of the rotation matrix (frame axis i): x=col0, y=col1, z=col2.
    geometry::Vector3d axis(int i) const {
        return {rotation[0][i], rotation[1][i], rotation[2][i]};
    }
    geometry::Vector3d x_axis() const { return axis(0); }
    geometry::Vector3d y_axis() const { return axis(1); }
    geometry::Vector3d z_axis() const { return axis(2); }
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_REFERENCE_FRAME_HPP
