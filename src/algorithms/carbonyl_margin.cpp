/**
 * @file carbonyl_margin.cpp
 * @brief Pyrimidine-O2 Watson/Sugar correction.
 */
#include <pairfinder/algorithms/carbonyl_margin.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <unordered_set>
#include <utility>
#include <vector>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms::classification {

namespace {

using core::Residue;
using geometry::Vector3d;

constexpr double kPi = 3.14159265358979323846;
constexpr double kHbMax = 3.5;

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::array<const char*, 6> kPyrRing = {"N1", "C2", "N3", "C4", "C5", "C6"};

const Vector3d* coord(const Residue& r, const std::string& a) {
    const core::Atom* at = r.get_atom(a);
    return at ? &at->coords : nullptr;
}

bool is_pyrimidine(const std::string& base, const core::BaseTyping& typing) {
    const std::string p = typing.normalize(base);
    return p == "C" || p == "U" || p == "T";
}

// Best-fit plane normal of the pyrimidine ring atoms (covariance smallest
// eigenvector). Sign is irrelevant (used only for in-plane projection).
std::optional<Vector3d> plane_normal(const Residue& res) {
    std::vector<Vector3d> pts;
    for (const char* a : kPyrRing)
        if (const Vector3d* c = coord(res, a)) pts.push_back(*c);
    if (pts.size() < 3) return std::nullopt;
    Vector3d mean;
    for (const auto& p : pts) mean = mean + p;
    mean = mean / static_cast<double>(pts.size());
    // Symmetric covariance matrix.
    double c[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (const auto& p : pts) {
        const double d[3] = {p.x - mean.x, p.y - mean.y, p.z - mean.z};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) c[i][j] += d[i] * d[j];
    }
    // Cyclic Jacobi eigensolver (3x3 symmetric).
    double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = std::fabs(c[0][1]) + std::fabs(c[0][2]) + std::fabs(c[1][2]);
        if (off < 1e-14) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::fabs(c[p][q]) < 1e-18) continue;
                const double theta = (c[q][q] - c[p][p]) / (2.0 * c[p][q]);
                const double t = (theta >= 0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double cs = 1.0 / std::sqrt(t * t + 1.0);
                const double sn = t * cs;
                for (int i = 0; i < 3; ++i) {
                    const double cip = c[i][p], ciq = c[i][q];
                    c[i][p] = cs * cip - sn * ciq;
                    c[i][q] = sn * cip + cs * ciq;
                }
                for (int i = 0; i < 3; ++i) {
                    const double cpi = c[p][i], cqi = c[q][i];
                    c[p][i] = cs * cpi - sn * cqi;
                    c[q][i] = sn * cpi + cs * cqi;
                }
                for (int i = 0; i < 3; ++i) {
                    const double vip = v[i][p], viq = v[i][q];
                    v[i][p] = cs * vip - sn * viq;
                    v[i][q] = sn * vip + cs * viq;
                }
            }
        }
    }
    int smallest = 0;
    for (int i = 1; i < 3; ++i)
        if (c[i][i] < c[smallest][smallest]) smallest = i;
    Vector3d n{v[0][smallest], v[1][smallest], v[2][smallest]};
    const double nn = n.norm();
    return nn == 0 ? std::nullopt : std::optional<Vector3d>(n / nn);
}

Vector3d in_plane(const Vector3d& v, const Vector3d& n) { return v - n * v.dot(n); }

int face_pr(char f) { return f == 'W' ? 0 : (f == 'H' ? 1 : 2); }

// Canonical LW class: returns (class, swapped).
std::pair<std::string, bool> canonicalize(char o, char f1, char f2) {
    std::string s(1, o);
    if (face_pr(f1) <= face_pr(f2)) return {s + f1 + f2, false};
    return {s + f2 + f1, true};
}

}  // namespace

std::optional<double> o2_watson_angle(const Residue& res_b, const Residue& partner_res,
                                      const core::BaseTyping& typing) {
    if (!is_pyrimidine(res_b.base_type(), typing)) return std::nullopt;
    const Vector3d* o2 = coord(res_b, "O2");
    const Vector3d* c2 = coord(res_b, "C2");
    const Vector3d* n3 = coord(res_b, "N3");
    if (!o2 || !c2 || !n3) return std::nullopt;

    const Vector3d* closest = nullptr;
    double best = std::numeric_limits<double>::infinity();
    for (const auto& at : partner_res.atoms()) {
        if (!kBaseAtoms.count(at.name)) continue;
        const double d = (at.coords - *o2).norm();
        if (d <= kHbMax && d < best) { best = d; closest = &at.coords; }
    }
    if (closest == nullptr) return std::nullopt;

    const auto n = plane_normal(res_b);
    if (!n) return std::nullopt;
    Vector3d axis = in_plane(*o2 - *c2, *n);
    if (axis.norm() == 0) return std::nullopt;
    axis = axis / axis.norm();
    Vector3d w_side = in_plane(*n3 - *c2, *n);
    Vector3d w_perp = w_side - axis * w_side.dot(axis);
    if (w_perp.norm() == 0) return std::nullopt;
    w_perp = w_perp / w_perp.norm();
    Vector3d pvec = in_plane(*closest - *o2, *n);
    if (pvec.norm() == 0) return std::nullopt;
    pvec = pvec / pvec.norm();
    return std::atan2(pvec.dot(w_perp), pvec.dot(axis)) * 180.0 / kPi;
}

O2Correction correct_cis_o2(const std::string& lw_class, const std::string& res_id1,
                            const std::string& res_id2, const Residue& res1,
                            const Residue& res2,
                            const std::unordered_set<std::string>& committed,
                            const core::BaseTyping& typing) {
    O2Correction out;
    if (lw_class.size() != 3 || lw_class[0] != 'c') return out;
    const char edge1 = lw_class[1], edge2 = lw_class[2];
    if (edge1 != 'S' && edge2 != 'S') return out;
    char new1 = edge1, new2 = edge2;
    if (edge1 == 'S' && edge2 != 'S' && !committed.count(res_id1)) {
        const auto a = o2_watson_angle(res1, res2, typing);
        if (a && *a >= kDecisiveWatsonDeg) new1 = 'W';
    }
    if (edge2 == 'S' && edge1 != 'S' && !committed.count(res_id2)) {
        const auto a = o2_watson_angle(res2, res1, typing);
        if (a && *a >= kDecisiveWatsonDeg) new2 = 'W';
    }
    if (new1 == edge1 && new2 == edge2) return out;
    const auto [new_class, swapped] = canonicalize(lw_class[0], new1, new2);
    if (new_class == lw_class) return out;
    out.fired = true;
    out.new_class = new_class;
    out.swapped = swapped;
    return out;
}

}  // namespace pairfinder::algorithms::classification
