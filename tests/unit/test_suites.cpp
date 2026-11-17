/**
 * Backbone suite classification (Richardson conformers): every conformer centre,
 * graded suiteness away from a centre, outliers, the delta pucker bands, and the
 * consistency of the two classification entry points.
 */
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/scoring/richardson.hpp>

#include "test_helper.hpp"

using pairfinder::scoring::RichardsonClassifier;
using Suite = std::array<double, 7>;

namespace {

struct Conformer {
    std::string name, bin;
    Suite angles{};
};

const std::filesystem::path& suites_file() {
    static const auto path = pairfinder::resources::file("richardson_suites.json");
    return path;
}

const RichardsonClassifier& classifier() {
    static const RichardsonClassifier c(suites_file());
    return c;
}

std::vector<Conformer> conformers() {
    std::ifstream in(suites_file());
    const auto j = nlohmann::json::parse(in);
    std::vector<Conformer> out;
    for (const auto& c : j["conformers"]) {
        Conformer conf;
        conf.name = c["name"].get<std::string>();
        conf.bin = c["bin"].get<std::string>();
        for (int k = 0; k < 7; ++k) conf.angles[k] = c["angles"][k].get<double>();
        out.push_back(conf);
    }
    return out;
}

double wrap(double a) {
    a = std::fmod(a, 360.0);
    return a < 0.0 ? a + 360.0 : a;
}

}  // namespace

/** The table has the 54 Richardson conformers. */
TEST(table_has_54_conformers) {
    CHECK_EQ(conformers().size(), std::size_t{54});
}

/** Each conformer's own centre is assigned to that conformer with suiteness 1. */
TEST(every_centre_is_its_own_conformer) {
    for (const auto& c : conformers()) {
        const auto r = classifier().classify(c.angles);
        CHECK_MSG(!r.is_outlier, c.name);
        CHECK_EQ(r.conformer, c.name);
        CHECK_EQ(r.bin, c.bin);
        CHECK_MSG(r.suiteness > 0.99, c.name << " suiteness " << r.suiteness);
        CHECK_NEAR(r.distance, 0.0, 1e-6);
    }
}

/** suiteness() and classify() agree on centres, near-centre and outlier suites. */
TEST(fast_and_full_paths_agree) {
    const auto confs = conformers();
    for (std::size_t i = 0; i < confs.size(); i += 5) {
        Suite s = confs[i].angles;
        for (double step : {0.0, 7.0, 15.0, 60.0}) {
            Suite t = s;
            t[4] = wrap(t[4] + step);  // beta
            CHECK_NEAR(classifier().suiteness(t), classifier().classify(t).suiteness, 1e-12);
        }
    }
}

/** Moving one torsion away from a centre lowers suiteness step by step to 0. */
TEST(suiteness_falls_with_distance) {
    const auto confs = conformers();
    const auto& a = confs.front();  // 1a, the A-form conformer
    CHECK_EQ(a.name, std::string("1a"));
    double previous = 1.01;
    for (double step : {0.0, 5.0, 10.0, 20.0, 35.0, 60.0}) {
        Suite s = a.angles;
        s[3] = wrap(s[3] + step);  // alpha
        const double v = classifier().suiteness(s);
        CHECK_MSG(v <= previous + 1e-12, "step " << step << ": " << v << " > " << previous);
        CHECK(v >= 0.0 && v <= 1.0);
        previous = v;
    }
    CHECK_MSG(previous < 0.2, "60 degrees off alpha still gives suiteness " << previous);
}

/** A suite far from every conformer is an outlier ("!!", suiteness 0) with a nearest conformer. */
TEST(far_suite_is_outlier) {
    const Suite odd{80.0, 100.0, 180.0, 180.0, 60.0, 300.0, 80.0};
    const auto r = classifier().classify(odd);
    CHECK(r.is_outlier);
    CHECK_EQ(r.conformer, std::string("!!"));
    CHECK_NEAR(r.suiteness, 0.0, 1e-12);
    CHECK(!classifier().nearest_conformer(odd).empty());
}

/** A delta between the C3'-endo and C2'-endo bands fails the sieve. */
TEST(delta_between_bands_is_outlier) {
    Suite s = conformers().front().angles;
    s[0] = 115.0;  // delta(i-1)
    const auto r = classifier().classify(s);
    CHECK(r.is_outlier);
    CHECK_NEAR(r.suiteness, 0.0, 1e-12);
}

/** Pucker from delta: C3'-endo band, C2'-endo band, and between or outside them. */
TEST(pucker_bands_from_delta) {
    CHECK_EQ(classifier().pucker(82.0), 3);
    CHECK_EQ(classifier().pucker(145.0), 2);
    CHECK_EQ(classifier().pucker(115.0), 0);
    CHECK_EQ(classifier().pucker(30.0), 0);
    CHECK_EQ(classifier().pucker(200.0), 0);
}

/** Conformer centres are looked up by name; unknown names give nothing. */
TEST(centre_lookup) {
    const auto confs = conformers();
    const auto* c = classifier().center_for(confs.front().name);
    CHECK(c != nullptr);
    for (int k = 0; k < 7; ++k) CHECK_NEAR((*c)[k], confs.front().angles[k], 1e-9);
    CHECK(classifier().center_for("not-a-conformer") == nullptr);
    for (double w : classifier().normal_widths()) CHECK(w > 0.0);
}

RUBRIC_TEST_MAIN()
