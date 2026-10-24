/**
 * @file rigid_body.cpp
 * @brief Rigid-body pair parameters.
 */
#include <pairfinder/algorithms/rigid_body.hpp>

#include <algorithm>
#include <array>
#include <cmath>

namespace pairfinder::algorithms {

namespace {

using geometry::Vector3d;
using Mat3 = std::array<std::array<double, 3>, 3>;

constexpr double kPi = 3.14159265358979323846;

Mat3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }

Mat3 transpose(const Mat3& m) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r[i][j] = m[j][i];
    return r;
}

Mat3 matmul(const Mat3& a, const Mat3& b) {
    Mat3 r{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double s = 0.0;
            for (int k = 0; k < 3; ++k) s += a[i][k] * b[k][j];
            r[i][j] = s;
        }
    return r;
}

Vector3d matvec(const Mat3& m, const Vector3d& v) {
    return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
}

// (axis unit vector, angle in radians) from a rotation matrix (incl. the
// theta~pi branch).
std::pair<Vector3d, double> matrix_to_axis_angle(const Mat3& R) {
    double cos_theta = std::clamp((R[0][0] + R[1][1] + R[2][2] - 1.0) / 2.0, -1.0, 1.0);
    const double theta = std::acos(cos_theta);
    if (theta < 1e-9) return {{0.0, 0.0, 1.0}, 0.0};
    if (std::abs(theta - kPi) < 1e-9) {
        Mat3 M{};
        const Mat3 I = identity();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) M[i][j] = (R[i][j] + I[i][j]) / 2.0;
        int idx = 0;
        if (M[1][1] > M[idx][idx]) idx = 1;
        if (M[2][2] > M[idx][idx]) idx = 2;
        const double denom = std::sqrt(std::max(M[idx][idx], 0.0));
        Vector3d axis{M[0][idx], M[1][idx], M[2][idx]};
        axis = axis / denom;
        axis = axis / axis.norm();
        return {axis, theta};
    }
    const double s = 2.0 * std::sin(theta);
    Vector3d axis{(R[2][1] - R[1][2]) / s, (R[0][2] - R[2][0]) / s, (R[1][0] - R[0][1]) / s};
    return {axis, theta};
}

// Rodrigues: matrix from axis-angle.
Mat3 axis_angle_to_matrix(const Vector3d& a, double angle) {
    if (angle == 0.0) return identity();
    const Mat3 K = {{{0.0, -a.z, a.y}, {a.z, 0.0, -a.x}, {-a.y, a.x, 0.0}}};
    const Mat3 KK = matmul(K, K);
    const double s = std::sin(angle), c = 1.0 - std::cos(angle);
    Mat3 r = identity();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r[i][j] += s * K[i][j] + c * KK[i][j];
    return r;
}

}  // namespace

RigidBodyParameters compute_rigid_body_parameters(const core::ReferenceFrame& frame1,
                                                  const core::ReferenceFrame& frame2) {
    const Mat3& R1 = frame1.rotation;
    // R2 = frame2.rotation @ diag(1,-1,-1): negate columns 1 and 2.
    Mat3 R2 = frame2.rotation;
    for (int i = 0; i < 3; ++i) {
        R2[i][1] = -R2[i][1];
        R2[i][2] = -R2[i][2];
    }

    const Mat3 R_rel = matmul(transpose(R1), R2);
    const auto [axis_local, theta] = matrix_to_axis_angle(R_rel);
    const Mat3 half_local = axis_angle_to_matrix(axis_local, theta / 2.0);
    const Mat3 R_mid = matmul(R1, half_local);
    const Mat3 R_mid_T = transpose(R_mid);

    const Vector3d disp = frame2.origin - frame1.origin;
    const Vector3d trans = matvec(R_mid_T, disp);
    const Vector3d axis_mid = matvec(R_mid_T, matvec(R1, axis_local));
    const double angle_deg = theta * 180.0 / kPi;

    RigidBodyParameters p;
    p.shear = trans.x;
    p.stretch = trans.y;
    p.stagger = trans.z;
    p.buckle = angle_deg * axis_mid.x;
    p.propeller = angle_deg * axis_mid.y;
    p.opening = angle_deg * axis_mid.z;
    return p;
}

}  // namespace pairfinder::algorithms
