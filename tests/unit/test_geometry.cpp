/**
 * Geometry building blocks: vectors, rigid superposition, base reference frames,
 * residue-ID parsing and angular gaps. Each test uses coordinates with a known
 * exact answer.
 */
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <pairfinder/geometry/superposition.hpp>
#include <pairfinder/geometry/vector3d.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <pairfinder/scoring/torsions.hpp>
#include <pairfinder/core/resource_locator.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

using pairfinder::geometry::Vector3d;
using pairfinder::geometry::superpose;

namespace {

constexpr double kPi = 3.14159265358979323846;

/// Rotation about an arbitrary unit axis (Rodrigues), row-major.
std::array<std::array<double, 3>, 3> rotation_about(Vector3d axis, double degrees) {
    axis = axis.normalized();
    const double t = degrees * kPi / 180.0, c = std::cos(t), s = std::sin(t), v = 1.0 - c;
    const double x = axis.x, y = axis.y, z = axis.z;
    return {{{c + x * x * v, x * y * v - z * s, x * z * v + y * s},
             {y * x * v + z * s, c + y * y * v, y * z * v - x * s},
             {z * x * v - y * s, z * y * v + x * s, c + z * z * v}}};
}

Vector3d apply(const std::array<std::array<double, 3>, 3>& r, const Vector3d& p, const Vector3d& t) {
    return {r[0][0] * p.x + r[0][1] * p.y + r[0][2] * p.z + t.x,
            r[1][0] * p.x + r[1][1] * p.y + r[1][2] * p.z + t.y,
            r[2][0] * p.x + r[2][1] * p.y + r[2][2] * p.z + t.z};
}

/// A non-planar set of seven points (roughly the size of a base).
std::vector<Vector3d> sample_points() {
    return {{0.0, 0.0, 0.0}, {1.4, 0.0, 0.0}, {2.1, 1.2, 0.0}, {1.4, 2.4, 0.1},
            {0.0, 2.4, -0.1}, {-0.7, 1.2, 0.0}, {3.5, 1.3, 0.8}};
}

}  // namespace

/** Vector arithmetic: dot, cross, norm and distance give their exact values. */
TEST(vector_arithmetic) {
    const Vector3d a{1, 2, 3}, b{4, -5, 6};
    CHECK_NEAR(a.dot(b), 12.0, 1e-12);
    const Vector3d c = a.cross(b);
    CHECK_NEAR(c.x, 27.0, 1e-12);
    CHECK_NEAR(c.y, 6.0, 1e-12);
    CHECK_NEAR(c.z, -13.0, 1e-12);
    CHECK_NEAR(Vector3d(3, 4, 0).norm(), 5.0, 1e-12);
    CHECK_NEAR(a.distance_to(b), std::sqrt(9.0 + 49.0 + 9.0), 1e-12);
    CHECK_NEAR(Vector3d(0, 0, 7).normalized().z, 1.0, 1e-12);
    CHECK_NEAR(Vector3d().normalized().norm(), 0.0, 1e-12);  // zero vector stays zero
}

/** A rotated and translated copy superposes back with RMSD 0 and the same rotation. */
TEST(superpose_recovers_rigid_motion) {
    const auto from = sample_points();
    const auto rot = rotation_about({1.0, 2.0, -0.5}, 73.0);
    const Vector3d shift{4.0, -2.5, 10.0};
    std::vector<Vector3d> to;
    for (const auto& p : from) to.push_back(apply(rot, p, shift));

    const auto r = superpose(from, to);
    CHECK(r.valid);
    CHECK_NEAR(r.rms, 0.0, 1e-8);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) CHECK_NEAR(r.rotation[i][j], rot[i][j], 1e-8);
    CHECK_NEAR(r.translation.x, shift.x, 1e-8);
    CHECK_NEAR(r.translation.y, shift.y, 1e-8);
    CHECK_NEAR(r.translation.z, shift.z, 1e-8);
}

/** Moving one point by 1 A gives the RMSD expected for a single displaced point. */
TEST(superpose_reports_misfit) {
    const auto from = sample_points();
    auto to = from;
    to[6].z += 1.0;
    const auto r = superpose(from, to);
    CHECK(r.valid);
    CHECK(r.rms > 0.05);
    CHECK(r.rms < 1.0 / std::sqrt(7.0) + 1e-9);  // never worse than leaving it unfitted
}

/** Fewer than three points, or unequal sets, cannot be superposed. */
TEST(superpose_rejects_bad_input) {
    const std::vector<Vector3d> two{{0, 0, 0}, {1, 0, 0}};
    CHECK(!superpose(two, two).valid);
    const auto pts = sample_points();
    const std::vector<Vector3d> fewer(pts.begin(), pts.begin() + 4);
    CHECK(!superpose(pts, fewer).valid);
}

/** Base frames of an idealized G-C pair are valid, orthonormal and right-handed. */
TEST(base_frames_are_orthonormal) {
    const auto path = pairfinder::resources::basepair_idealized() / "cWW" / "GC.pdb";
    const auto prep = rubric_test::pipeline().prepare(path.string());
    CHECK_EQ(prep.frames.size(), std::size_t{2});
    for (const auto& [res_id, f] : prep.frames) {
        CHECK_MSG(f.valid, res_id);
        CHECK_MSG(f.rms_fit < 0.1, res_id << " rms_fit " << f.rms_fit);
        for (int i = 0; i < 3; ++i) {
            CHECK_NEAR(f.axis(i).norm(), 1.0, 1e-9);
            for (int j = i + 1; j < 3; ++j) CHECK_NEAR(f.axis(i).dot(f.axis(j)), 0.0, 1e-9);
        }
        CHECK_NEAR(f.x_axis().cross(f.y_axis()).dot(f.z_axis()), 1.0, 1e-9);
    }
    // In a cWW pair the two base normals are antiparallel.
    const double dz = prep.frames[0].second.z_axis().dot(prep.frames[1].second.z_axis());
    CHECK_MSG(dz < -0.9, "z-axis dot " << dz);
}

/** Residue IDs: plain, negative and insertion-code numbers, and an ID with no number. */
TEST(parse_residue_ids) {
    using pairfinder::scoring::parse_res_seq;
    auto a = parse_res_seq("A-G-12");
    CHECK(a.has_value());
    CHECK_EQ(a->chain, std::string("A"));
    CHECK_EQ(a->num, 12L);
    CHECK_EQ(a->icode, std::string(""));
    auto b = parse_res_seq("G-C--8");
    CHECK(b.has_value());
    CHECK_EQ(b->num, -8L);
    auto c = parse_res_seq("A-G-20A");
    CHECK(c.has_value());
    CHECK_EQ(c->num, 20L);
    CHECK_EQ(c->icode, std::string("A"));
    auto d = parse_res_seq("AA-PSU-1050");
    CHECK(d.has_value());
    CHECK_EQ(d->chain, std::string("AA"));
    CHECK_EQ(d->num, 1050L);
    CHECK(!parse_res_seq("A-G").has_value());
}

/** Angular gaps take the short way round the circle. */
TEST(signed_angle_gap_wraps) {
    using pairfinder::scoring::signed_angle_gap;
    CHECK_NEAR(signed_angle_gap(10.0, 350.0), 20.0, 1e-12);
    CHECK_NEAR(signed_angle_gap(350.0, 10.0), -20.0, 1e-12);
    CHECK_NEAR(signed_angle_gap(200.0, 150.0), 50.0, 1e-12);
    CHECK_NEAR(signed_angle_gap(90.0, 90.0), 0.0, 1e-12);
    const double g = signed_angle_gap(0.0, 180.0);
    CHECK_NEAR(std::fabs(g), 180.0, 1e-12);
}

RUBRIC_TEST_MAIN()
