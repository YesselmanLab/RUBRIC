/**
 * Regression panel: ten real structures chosen for the cases that are hard to get
 * right (listed in tests/data/panel/README.md). For each, the base pairs (residues,
 * class, ambiguity, score, tier), structure scores, every suite (conformer,
 * suiteness), backbone recommendations and pucker flags are compared with the
 * stored results in tests/data/panel/expected/. Scores and suiteness are compared
 * within small tolerances so the test passes on every platform and compiler.
 *
 * After an intended change to RUBRIC's results, rewrite the stored results with
 *   RUBRIC_UPDATE_PANEL=1 ctest --test-dir build -R test_panel
 * and review the difference in git before committing it.
 */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

namespace fs = std::filesystem;
using nlohmann::json;

namespace {

const fs::path kPanel = rubric_test::tests_dir() / "data" / "panel";
constexpr double kScoreTol = 0.05;      // pair and structure scores (0-100)
constexpr double kSuitenessTol = 0.002; // per-suite and structure suiteness (0-1)

json round_to(double v, int digits) {
    const double f = std::pow(10.0, digits);
    return std::round(v * f) / f;
}

/// Everything the panel compares, in a stable order.
json summarize(const rubric_test::EngineResult& r) {
    json pairs = json::array();
    for (const auto& c : r.selected) {
        json score = nullptr, tier = nullptr;
        for (const auto& p : r.score.pair_scores) {
            if (p.res_id1 == c.res_id1 && p.res_id2 == c.res_id2) {
                score = round_to(p.score, 2);
                tier = p.tier();
            }
        }
        pairs.push_back({c.res_id1, c.res_id2, c.lw_class, c.is_ambiguous, score, tier});
    }
    std::sort(pairs.begin(), pairs.end());

    json suites = json::array();
    for (const auto& t : r.score.backbone_torsions)
        if (t.scored) suites.push_back({t.res_id, t.conformer, round_to(t.suiteness, 3)});
    std::sort(suites.begin(), suites.end());

    json recs = json::array(), puckers = json::array(), unclassified = json::array();
    for (const auto& rs : r.score.residue_scores) {
        if (rs.recommendation) recs.push_back({rs.res_id, rs.recommendation->tier, rs.recommendation->target_conformer});
        if (rs.pucker_outlier) puckers.push_back(rs.res_id);
        if (rs.pucker_unclassified) unclassified.push_back(rs.res_id);
    }
    std::sort(recs.begin(), recs.end());
    std::sort(puckers.begin(), puckers.end());
    std::sort(unclassified.begin(), unclassified.end());

    const auto t = r.score.tier_summary();
    return {
        {"n_pairs", r.selected.size()},
        {"skipped_pairs", r.score.skipped_pairs},
        {"pairs_score", r.score.pair_scores.empty() ? json(nullptr) : round_to(r.score.pairs_score, 2)},
        {"suiteness", r.score.residue_scores.empty() ? json(nullptr) : round_to(r.score.residues_score / 100.0, 3)},
        {"tiers", {t.preferred, t.acceptable, t.review}},
        {"pairs", pairs},
        {"suites", suites},
        {"recommendations", recs},
        {"pucker_outliers", puckers},
        {"pucker_unclassified", unclassified},
    };
}

bool close(const json& a, const json& b, double tol) {
    if (a.is_null() || b.is_null()) return a.is_null() && b.is_null();
    return std::fabs(a.get<double>() - b.get<double>()) <= tol;
}

/// Compare a got/expected list of rows whose numeric column `num` is compared within tol.
void compare_rows(const std::string& id, const std::string& what, const json& got, const json& exp,
                  int num, double tol) {
    CHECK_MSG(got.size() == exp.size(), id << " " << what << ": " << got.size() << " vs " << exp.size());
    for (std::size_t i = 0; i < got.size(); ++i) {
        for (std::size_t k = 0; k < got[i].size(); ++k) {
            if (static_cast<int>(k) == num) {
                CHECK_MSG(close(got[i][k], exp[i][k], tol), id << " " << what << " " << exp[i].dump() << " got " << got[i].dump());
            } else {
                CHECK_MSG(got[i][k] == exp[i][k], id << " " << what << " " << exp[i].dump() << " got " << got[i].dump());
            }
        }
    }
}

std::vector<fs::path> panel_files() {
    // 1GID lives in tests/data (shared with the package check); the others in tests/data/panel.
    std::vector<fs::path> files{rubric_test::tests_dir() / "data" / "1GID.cif.gz"};
    for (const auto& f : fs::directory_iterator(kPanel))
        if (f.path().string().size() > 7 && f.path().string().ends_with(".cif.gz")) files.push_back(f.path());
    std::sort(files.begin(), files.end());
    return files;
}

}  // namespace

/** The panel has its ten structures and a stored result for each. */
TEST(panel_complete) {
    const auto files = panel_files();
    CHECK_EQ(files.size(), std::size_t{10});
    if (std::getenv("RUBRIC_UPDATE_PANEL")) return;
    for (const auto& f : files) {
        const std::string id = f.filename().string().substr(0, 4);
        CHECK_MSG(fs::exists(kPanel / "expected" / (id + ".json")), id);
    }
}

/** Every panel structure gives the stored results. */
TEST(panel_matches_expected) {
    const bool update = std::getenv("RUBRIC_UPDATE_PANEL") != nullptr;
    for (const auto& f : panel_files()) {
        const std::string id = f.filename().string().substr(0, 4);
        const json got = summarize(rubric_test::run_engine(f));
        const auto expected_path = kPanel / "expected" / (id + ".json");
        if (update) {
            std::ofstream(expected_path) << got.dump(1) << "\n";
            std::cout << "  wrote " << expected_path.filename().string() << "\n";
            continue;
        }
        const json exp = json::parse(rubric_test::read_text(expected_path));
        for (const char* key : {"n_pairs", "skipped_pairs", "tiers", "recommendations", "pucker_outliers",
                                "pucker_unclassified"})
            CHECK_MSG(got[key] == exp[key], id << " " << key << ": " << got[key].dump() << " vs " << exp[key].dump());
        CHECK_MSG(close(got["pairs_score"], exp["pairs_score"], kScoreTol), id << " pairs_score");
        CHECK_MSG(close(got["suiteness"], exp["suiteness"], kSuitenessTol), id << " suiteness");
        compare_rows(id, "pair", got["pairs"], exp["pairs"], 4, kScoreTol);
        compare_rows(id, "suite", got["suites"], exp["suites"], 2, kSuitenessTol);
    }
}

RUBRIC_TEST_MAIN()
