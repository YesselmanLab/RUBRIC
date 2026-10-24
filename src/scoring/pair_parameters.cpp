/**
 * @file pair_parameters.cpp
 * @brief Six base-pair parameters.
 */
#include <pairfinder/scoring/pair_parameters.hpp>

#include <algorithm>
#include <array>
#include <cmath>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::scoring {

namespace {

using geometry::Vector3d;
using Mat3 = std::array<std::array<double, 3>, 3>;  // row-major

constexpr double kXeps = 1e-7;
constexpr double kPi = 3.14159265358979323846;

double deg2rad(double a) { return a * kPi / 180.0; }
double rad2deg(double a) { return a * 180.0 / kPi; }

Vector3d col(const Mat3& m, int i) { return {m[0][i], m[1][i], m[2][i]}; }

Vector3d mul(const Mat3& m, const Vector3d& v) {
    return {m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z,
            m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z,
            m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z};
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

double magang(const Vector3d& va, const Vector3d& vb) {
    const double la = va.norm(), lb = vb.norm();
    if (la < kXeps || lb < kXeps) return 0.0;
    double dot = (va / la).dot(vb / lb);
    dot = std::max(-1.0, std::min(1.0, dot));
    return rad2deg(std::acos(dot));
}

Mat3 arb_rotation(const Vector3d& axis, double angle_deg) {
    const double vlen = axis.norm();
    if (vlen < kXeps) return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const Vector3d a = axis / vlen;
    const double c = std::cos(deg2rad(angle_deg)), s = std::sin(deg2rad(angle_deg));
    const double dc = 1.0 - c;
    const double x = a.x, y = a.y, z = a.z;
    return {{{c + dc * x * x, x * y * dc - z * s, x * z * dc + y * s},
             {x * y * dc + z * s, c + dc * y * y, y * z * dc - x * s},
             {x * z * dc - y * s, y * z * dc + x * s, c + dc * z * z}}};
}

double vec_ang(const Vector3d& va, const Vector3d& vb, const Vector3d& vref) {
    const Vector3d n = vref / vref.norm();
    Vector3d va_p = va - n * va.dot(n);
    Vector3d vb_p = vb - n * vb.dot(n);
    const double la = va_p.norm(), lb = vb_p.norm();
    if (la < kXeps || lb < kXeps) return 0.0;
    va_p = va_p / la;
    vb_p = vb_p / lb;
    double ang = magang(va_p, vb_p);
    if (va_p.cross(vb_p).dot(n) < 0.0) ang = -ang;
    return ang;
}

Vector3d get_vector(const Vector3d& va, const Vector3d& vref, double deg_ang) {
    const Vector3d n = vref / vref.norm();
    Vector3d va_p = va;
    if (va_p.dot(n) > kXeps) va_p = va_p - n * va_p.dot(n);
    const Vector3d vo = mul(arb_rotation(n, deg_ang), va_p);
    const double l = vo.norm();
    return l < kXeps ? vo : vo / l;
}

Mat3 to_mat(const core::ReferenceFrame& f) { return f.rotation; }

}  // namespace

PairParameters compute_pair_parameters(const core::ReferenceFrame& frame1,
                                       const core::ReferenceFrame& frame2) {
    Mat3 r1 = to_mat(frame1);
    Mat3 r2 = to_mat(frame2);
    const Vector3d o1 = frame1.origin;
    const Vector3d o2 = frame2.origin;

    if (col(r1, 2).dot(col(r2, 2)) < 0.0)  // anti-parallel z: negate y,z columns of r2
        for (int row = 0; row < 3; ++row) { r2[row][1] = -r2[row][1]; r2[row][2] = -r2[row][2]; }

    const Mat3& r_first = r2;
    const Mat3& r_second = r1;
    const Vector3d o_first = o2;
    const Vector3d o_second = o1;

    const Vector3d z1 = col(r_first, 2);
    const Vector3d z2 = col(r_second, 2);
    Vector3d hinge = z1.cross(z2);
    const double bp_mag = magang(z1, z2);

    if (hinge.norm() < kXeps && (std::abs(bp_mag - 180.0) < kXeps || bp_mag < kXeps))
        hinge = col(r_first, 0) + col(r_first, 1) + col(r_second, 0) + col(r_second, 1);

    const Mat3 para2 = matmul(arb_rotation(hinge, -0.5 * bp_mag), r_second);
    const Mat3 para1 = matmul(arb_rotation(hinge, 0.5 * bp_mag), r_first);

    const Vector3d mean_z = col(para2, 2);
    const Vector3d y1_para = col(para1, 1);
    const Vector3d y2_para = col(para2, 1);
    const double opening = vec_ang(y1_para, y2_para, mean_z);

    const Vector3d mean_y = get_vector(y1_para, mean_z, 0.5 * opening);
    const Vector3d mean_x = mean_y.cross(mean_z);

    const Vector3d disp = o_second - o_first;
    PairParameters p;
    p.shear = disp.dot(mean_x);
    p.stretch = disp.dot(mean_y);
    p.stagger = disp.dot(mean_z);
    p.opening = opening;

    const double phi = deg2rad(vec_ang(hinge, mean_y, mean_z));
    p.propeller = bp_mag * std::cos(phi);
    p.buckle = bp_mag * std::sin(phi);
    return p;
}

}  // namespace pairfinder::scoring
