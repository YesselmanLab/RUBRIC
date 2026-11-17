/**
 * Backbone suites, their recommendations, and ribose puckers, computed by the
 * scorer on real coordinates: an 8-nucleotide fragment of 1GID (chain A
 * 195-202) and the whole 1GID structure. Checks the rules that define each
 * reported field rather than individual numbers.
 */
#include <cmath>
#include <map>
#include <sstream>
#include <string>

#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/scoring/richardson.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

using rubric_test::run_engine;

namespace {

const std::filesystem::path kFragment = rubric_test::tests_dir() / "fragments" / "1GID_A195-202.pdb";
const std::filesystem::path k1GID = rubric_test::tests_dir() / "data" / "1GID.cif.gz";

// C3'-endo and C2'-endo delta bands, and the P-perpendicular cutoff (Richardson / MolProbity).
bool in_c3_band(double d) { return d >= 60.0 && d <= 105.0; }
bool in_c2_band(double d) { return d >= 125.0 && d <= 165.0; }
constexpr double kPperpCut = 2.9;

/// Copy of a PDB text without the atoms of one residue whose name is in `atoms`
/// (or with them renamed to `rename_to`, if given).
std::string edit_atom(const std::string& pdb, int res_num, const std::string& atom,
                      const std::string& rename_to = "") {
    std::istringstream in(pdb);
    std::ostringstream out;
    std::string line;
    while (std::getline(in, line)) {
        const bool is_atom = line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0;
        if (is_atom && line.size() > 26 && std::stoi(line.substr(22, 4)) == res_num) {
            std::string name = line.substr(12, 4);
            name.erase(0, name.find_first_not_of(' '));
            name.erase(name.find_last_not_of(' ') + 1);
            if (name == atom) {
                if (rename_to.empty()) continue;
                std::string padded = " " + rename_to;
                padded.resize(4, ' ');
                line = line.substr(0, 12) + padded + line.substr(16);
                if (line.size() >= 78) line = line.substr(0, 76) + " F" + line.substr(78);
            }
        }
        out << line << "\n";
    }
    return out.str();
}

std::map<std::string, pairfinder::scoring::BackboneTorsion> suites_by_residue(
        const pairfinder::scoring::StructureScore& s) {
    std::map<std::string, pairfinder::scoring::BackboneTorsion> m;
    for (const auto& t : s.backbone_torsions)
        if (t.scored) m[t.res_id] = t;
    return m;
}

}  // namespace

/** The fragment's first nucleotide has no suite; the other seven do. */
TEST(fragment_suites) {
    const auto r = run_engine(kFragment);
    const auto suites = suites_by_residue(r.score);
    CHECK_EQ(suites.size(), std::size_t{7});
    CHECK(suites.count("A-G-195") == 0);
    CHECK_EQ(suites.at("A-U-202").conformer, std::string("1a"));
    CHECK(suites.at("A-U-202").suiteness > 0.8);
    CHECK(suites.at("A-G-200").is_outlier);
    CHECK_NEAR(suites.at("A-G-200").suiteness, 0.0, 1e-12);
    CHECK(r.score.residues_score > 0.0 && r.score.residues_score < 100.0);
}

/** Recommendation labels follow their definitions, and targets move toward the conformer centre. */
TEST(recommendation_rules) {
    pairfinder::scoring::RichardsonClassifier suites(pairfinder::resources::file("richardson_suites.json"));
    int fixable = 0, review = 0, flag_only = 0;
    for (const auto& path : {kFragment, k1GID}) {
        const auto r = run_engine(path);
        const auto by_res = suites_by_residue(r.score);
        for (const auto& rs : r.score.residue_scores) {
            const auto it = by_res.find(rs.res_id);
            if (it == by_res.end()) continue;
            const auto& t = it->second;
            if (!rs.recommendation) {
                // Only assigned suites with suiteness of at least 0.5 go without a recommendation.
                CHECK_MSG(!t.is_outlier && t.suiteness >= 0.5, rs.res_id);
                continue;
            }
            const auto& rec = *rs.recommendation;
            if (rec.tier == "fixable") {
                ++fixable;
                CHECK_MSG(!t.is_outlier && t.suiteness < 0.5 && !rec.deviations.empty(), rs.res_id);
                CHECK_EQ(rec.target_conformer, t.conformer);
            } else if (rec.tier == "review") {
                ++review;
                CHECK_MSG(!t.is_outlier && t.suiteness < 0.5 && rec.deviations.empty(), rs.res_id);
            } else {
                CHECK_EQ(rec.tier, std::string("flag_only"));
                ++flag_only;
                CHECK_MSG(t.is_outlier && t.suiteness == 0.0, rs.res_id);
                CHECK_MSG(suites.center_for(rec.target_conformer) != nullptr, rs.res_id);
            }
            for (const auto& d : rec.deviations) {
                CHECK(d.conformer_center.has_value());
                const double before = std::fabs(pairfinder::scoring::signed_angle_gap(d.value, *d.conformer_center));
                const double after = std::fabs(pairfinder::scoring::signed_angle_gap(d.target, *d.conformer_center));
                CHECK_MSG(after <= before + 1e-6, rs.res_id << " " << d.angle << ": target moves away");
                CHECK_NEAR(before, std::fabs(d.gap), 0.06);
            }
        }
    }
    CHECK_MSG(fixable > 0 && review > 0 && flag_only > 0,
              "fixable " << fixable << ", review " << review << ", flag_only " << flag_only);
}

/** U199 is flagged: C3'-endo by delta but C2'-endo by its P-perpendicular distance. */
TEST(pucker_outlier_u199) {
    const auto r = run_engine(kFragment);
    bool found = false;
    for (const auto& rs : r.score.residue_scores) {
        if (rs.res_id != "A-U-199") {
            CHECK_MSG(!rs.pucker_outlier, rs.res_id);
            continue;
        }
        found = true;
        CHECK(rs.pucker_outlier.has_value());
        CHECK_EQ(rs.pucker_outlier->pucker_by_delta, std::string("C3'-endo"));
        CHECK_EQ(rs.pucker_outlier->pucker_by_pperp, std::string("C2'-endo"));
        CHECK(rs.pucker_outlier->pperp < kPperpCut);
    }
    CHECK(found);
}

/** Every pucker flag in 1GID is consistent with the delta bands and the P-perp cutoff. */
TEST(pucker_rules) {
    const auto r = run_engine(k1GID);
    int outliers = 0, unclassified = 0;
    for (const auto& rs : r.score.residue_scores) {
        CHECK_MSG(!(rs.pucker_outlier && rs.pucker_unclassified), rs.res_id);
        if (rs.pucker_outlier) {
            ++outliers;
            const auto& p = *rs.pucker_outlier;
            CHECK(p.delta.has_value());
            const bool c3_short = in_c3_band(*p.delta) && p.pperp < kPperpCut;
            const bool c2_long = in_c2_band(*p.delta) && p.pperp >= kPperpCut;
            CHECK_MSG(c3_short || c2_long, rs.res_id << " delta " << *p.delta << " pperp " << p.pperp);
        }
        if (rs.pucker_unclassified) {
            ++unclassified;
            const double d = rs.pucker_unclassified->delta;
            CHECK_MSG(!in_c3_band(d) && !in_c2_band(d), rs.res_id << " delta " << d);
        }
    }
    CHECK_EQ(outliers, 3);  // A-U-199, B-G-129, B-A-123
    CHECK(unclassified > 0);
}

/** A nucleotide without its 2'-O gets no suite, and neither does the suite that follows it. */
TEST(no_suite_without_2prime_oxygen) {
    const std::string pdb = rubric_test::read_text(kFragment);
    for (const std::string& rename : {std::string(""), std::string("F2'")}) {
        const auto path = rubric_test::write_temp(
            rename.empty() ? "no_o2prime.pdb" : "f2prime.pdb", edit_atom(pdb, 199, "O2'", rename));
        const auto suites = suites_by_residue(run_engine(path).score);
        CHECK_MSG(suites.count("A-U-199") == 0, "suite of the modified nucleotide is still scored");
        CHECK_MSG(suites.count("A-G-200") == 0, "the following suite is still scored");
        CHECK_MSG(suites.count("A-A-198") == 1, "an unaffected suite was dropped");
    }
}

/** A missing backbone atom drops only the suites that need it, without failing. */
TEST(missing_phosphorus) {
    const auto path = rubric_test::write_temp(
        "no_p200.pdb", edit_atom(rubric_test::read_text(kFragment), 200, "P"));
    const auto suites = suites_by_residue(run_engine(path).score);
    CHECK(suites.count("A-G-200") == 0);
    CHECK(suites.count("A-U-202") == 1);
}

RUBRIC_TEST_MAIN()
