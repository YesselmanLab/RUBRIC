/**
 * @file superposition.cpp
 * @brief Least-squares rigid superposition (as in 3DNA).
 *
 * Method: covariance matrix -> 4x4 quaternion key matrix (Horn) -> largest
 * eigenvector via symmetric Jacobi -> rotation matrix -> translation -> RMS.
 * Full double precision (no rounding). Convention:
 *   to_i ~= R * from_i + t   (from = standard/template, to = experimental).
 */
#include <pairfinder/geometry/superposition.hpp>

#include <array>
#include <cmath>
#include <utility>

namespace pairfinder::geometry {

namespace {

using Mat3 = std::array<std::array<double, 3>, 3>;
using Mat4 = std::array<std::array<double, 4>, 4>;
using Vec4 = std::array<double, 4>;

double comp(const Vector3d& v, int i) { return i == 0 ? v.x : (i == 1 ? v.y : v.z); }

Vector3d centroid(const std::vector<Vector3d>& pts) {
    Vector3d s{0.0, 0.0, 0.0};
    for (const auto& p : pts) s = s + p;
    return s / static_cast<double>(pts.size());
}

// U[i][j] = (1/(n-1)) * sum (p1_i - c1_i)(p2_j - c2_j)   (i: from, j: to)
Mat3 covariance(const std::vector<Vector3d>& p1, const std::vector<Vector3d>& p2,
                const Vector3d& c1, const Vector3d& c2) {
    const std::size_t n = p1.size();
    Mat3 u{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k)
                sum += comp(p1[k] - c1, i) * comp(p2[k] - c2, j);
            u[i][j] = sum / static_cast<double>(n - 1);
        }
    return u;
}

// 4x4 symmetric quaternion key matrix from the 3x3 covariance (Horn).
Mat4 build_quaternion_matrix(const Mat3& u) {
    Mat4 n{};
    const double u11 = u[0][0], u22 = u[1][1], u33 = u[2][2];
    n[0][0] = u11 + u22 + u33;
    n[1][1] = u11 - u22 - u33;
    n[2][2] = -u11 + u22 - u33;
    n[3][3] = -u11 - u22 + u33;
    const double u12 = u[0][1], u21 = u[1][0], u13 = u[0][2], u31 = u[2][0],
                 u23 = u[1][2], u32 = u[2][1];
    n[0][1] = n[1][0] = u23 - u32;
    n[0][2] = n[2][0] = u31 - u13;
    n[0][3] = n[3][0] = u12 - u21;
    n[1][2] = n[2][1] = u12 + u21;
    n[1][3] = n[3][1] = u31 + u13;
    n[2][3] = n[3][2] = u23 + u32;
    return n;
}

void rotate4(Mat4& a, int i, int j, int k, int l, double s, double tau) {
    const double g = a[i][j], h = a[k][l];
    a[i][j] = g - s * (h + g * tau);
    a[k][l] = h + s * (g - h * tau);
}

// Largest-eigenvalue eigenvector of a 4x4 symmetric matrix (Jacobi, with 3DNA's
// iteration schedule); eigenvalues sorted ascending, last column returned.
Vec4 largest_eigenvector(const Mat4& nin) {
    constexpr double kEps = 1.0e-7;
    constexpr int kMaxIter = 100;
    Mat4 a = nin, v{};
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) v[i][j] = (i == j) ? 1.0 : 0.0;
    std::array<double, 4> d{}, b{}, z{};
    for (int i = 0; i < 4; ++i) {
        b[i] = d[i] = a[i][i];
        z[i] = 0.0;
    }
    for (int iter = 0; iter < kMaxIter; ++iter) {
        double sm = 0.0;
        for (int i = 0; i < 3; ++i)
            for (int j = i + 1; j < 4; ++j) sm += std::abs(a[i][j]);
        if (sm < kEps) break;
        const double tresh = (iter < 4) ? 0.2 * sm / 16.0 : 0.0;
        for (int ip = 0; ip < 3; ++ip)
            for (int iq = ip + 1; iq < 4; ++iq) {
                const double g = 100.0 * std::abs(a[ip][iq]);
                if (iter > 4 && (std::abs(d[ip]) + g) == std::abs(d[ip]) &&
                    (std::abs(d[iq]) + g) == std::abs(d[iq])) {
                    a[ip][iq] = 0.0;
                    continue;
                }
                if (std::abs(a[ip][iq]) <= tresh) continue;
                double h = d[iq] - d[ip];
                double t;
                if ((std::abs(h) + g) == std::abs(h)) {
                    t = a[ip][iq] / h;
                } else {
                    const double theta = 0.5 * h / a[ip][iq];
                    t = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
                    if (theta < 0.0) t = -t;
                }
                const double c = 1.0 / std::sqrt(1.0 + t * t);
                const double s = t * c;
                const double tau = s / (1.0 + c);
                h = t * a[ip][iq];
                z[ip] -= h;
                z[iq] += h;
                d[ip] -= h;
                d[iq] += h;
                a[ip][iq] = 0.0;
                for (int j = 0; j < ip; ++j) rotate4(a, j, ip, j, iq, s, tau);
                for (int j = ip + 1; j < iq; ++j) rotate4(a, ip, j, j, iq, s, tau);
                for (int j = iq + 1; j < 4; ++j) rotate4(a, ip, j, iq, j, s, tau);
                for (int j = 0; j < 4; ++j) rotate4(v, j, ip, j, iq, s, tau);
            }
        for (int i = 0; i < 4; ++i) {
            b[i] += z[i];
            d[i] = b[i];
            z[i] = 0.0;
        }
    }
    // selection sort ascending; largest eigenvalue ends in column 3
    for (int i = 0; i < 3; ++i) {
        int k = i;
        double p = d[i];
        for (int j = i + 1; j < 4; ++j)
            if (d[j] < p) { p = d[j]; k = j; }
        if (k != i) {
            std::swap(d[i], d[k]);
            for (int j = 0; j < 4; ++j) std::swap(v[j][i], v[j][k]);
        }
    }
    Vec4 q;
    for (int i = 0; i < 4; ++i) q[i] = v[i][3];
    return q;
}

Mat3 quaternion_to_rotation(const Vec4& q) {
    Mat4 n;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) n[i][j] = q[i] * q[j];
    Mat3 r;
    r[0][0] = n[0][0] + n[1][1] - n[2][2] - n[3][3];
    r[0][1] = 2.0 * (n[1][2] - n[0][3]);
    r[0][2] = 2.0 * (n[1][3] + n[0][2]);
    r[1][0] = 2.0 * (n[2][1] + n[0][3]);
    r[1][1] = n[0][0] - n[1][1] + n[2][2] - n[3][3];
    r[1][2] = 2.0 * (n[2][3] - n[0][1]);
    r[2][0] = 2.0 * (n[3][1] - n[0][2]);
    r[2][1] = 2.0 * (n[3][2] + n[0][1]);
    r[2][2] = n[0][0] - n[1][1] - n[2][2] + n[3][3];
    return r;
}

Vector3d matvec(const Mat3& r, const Vector3d& v) {
    return {r[0][0] * v.x + r[0][1] * v.y + r[0][2] * v.z,
            r[1][0] * v.x + r[1][1] * v.y + r[1][2] * v.z,
            r[2][0] * v.x + r[2][1] * v.y + r[2][2] * v.z};
}

}  // namespace

SuperpositionResult superpose(const std::vector<Vector3d>& from,
                              const std::vector<Vector3d>& to) {
    SuperpositionResult r;
    if (from.size() != to.size() || from.size() < 3) return r;  // valid stays false
    const std::size_t n = from.size();
    const Vector3d c1 = centroid(from), c2 = centroid(to);
    const Mat3 u = covariance(from, to, c1, c2);
    const Mat4 key = build_quaternion_matrix(u);
    const Vec4 q = largest_eigenvector(key);
    const Mat3 rot = quaternion_to_rotation(q);
    const Vector3d t = c2 - matvec(rot, c1);
    double ss = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const Vector3d diff = (matvec(rot, from[i]) + t) - to[i];
        ss += diff.x * diff.x + diff.y * diff.y + diff.z * diff.z;
    }
    r.rotation = rot;
    r.translation = t;
    r.rms = std::sqrt(ss / static_cast<double>(n));
    r.valid = true;
    return r;
}

}  // namespace pairfinder::geometry
