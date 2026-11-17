/**
 * The `rubric` command line: every user option, the JSON fields it documents,
 * the PyMOL script, error exits, and the developer `dump-*` commands. Each test
 * runs the built program on tests/data/1GID.cif.gz, offline.
 */
#include <cstdio>
#include <filesystem>
#include <cmath>
#include <string>

#ifndef _WIN32
#include <sys/wait.h>
#endif

#include <nlohmann/json.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

namespace fs = std::filesystem;
using nlohmann::json;

namespace {

struct Run {
    int code = -1;
    std::string out;  ///< standard output and standard error
};

/// Run the rubric program with arguments; returns its exit code and output.
Run rubric(const std::string& args) {
    const std::string cmd = std::string("\"") + RUBRIC_EXE + "\" " + args + " 2>&1";
#ifdef _WIN32
    FILE* pipe = _popen(("\"" + cmd + "\"").c_str(), "r");
#else
    FILE* pipe = popen(cmd.c_str(), "r");
#endif
    Run r;
    if (!pipe) return r;
    char buf[4096];
    std::size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, pipe)) > 0) r.out.append(buf, n);
#ifdef _WIN32
    r.code = _pclose(pipe);
#else
    const int status = pclose(pipe);
    r.code = WIFEXITED(status) ? WEXITSTATUS(status) : -1;
#endif
    return r;
}

std::string q(const fs::path& p) { return "\"" + p.string() + "\""; }

const fs::path k1GID = rubric_test::tests_dir() / "data" / "1GID.cif.gz";

json run_json(const std::string& extra, const std::string& name) {
    const auto out = rubric_test::write_temp(name, "");
    const auto r = rubric(q(k1GID) + " --no-download " + extra + " --out " + q(out));
    CHECK_MSG(r.code == 0, r.out);
    return json::parse(rubric_test::read_text(out));
}

}  // namespace

/** --version prints the project version. */
TEST(version) {
    const auto r = rubric("--version");
    CHECK_EQ(r.code, 0);
    CHECK_MSG(r.out.find(RUBRIC_VERSION_EXPECTED) != std::string::npos, r.out);
}

/** --help describes the usage. */
TEST(help) {
    const auto r = rubric("--help");
    CHECK_EQ(r.code, 0);
    CHECK(r.out.find("--out") != std::string::npos);
    CHECK(r.out.find("--pml") != std::string::npos);
}

/** The report carries every documented top-level field, with consistent counts. */
TEST(report_fields) {
    const json d = run_json("", "report.json");
    for (const char* key : {"pdb_id", "n_pairs", "pairs", "pairs_score", "backbone_suiteness",
                            "tier_summary", "backbone_residues", "pucker_outliers", "clashes",
                            "n_clashes", "candidates_total", "candidates_valid"})
        CHECK_MSG(d.contains(key), key);
    CHECK_EQ(d["pdb_id"].get<std::string>(), std::string("1GID"));
    CHECK_EQ(d["n_pairs"].get<std::size_t>(), d["pairs"].size());
    CHECK_EQ(d["n_clashes"].get<std::size_t>(), d["clashes"].size());
    const auto& t = d["tier_summary"];
    for (const char* tier : {"Preferred", "Secondary review", "Priority review"}) CHECK_MSG(t.contains(tier), tier);
    for (const auto& p : d["pairs"]) {
        for (const char* key : {"res_id1", "res_id2", "bp_type", "classification", "score", "tier",
                                "issues", "issue_details", "num_hbonds"})
            CHECK_MSG(p.contains(key), key);
    }
    // Suiteness is reported with at most two decimals.
    const double s = d["backbone_suiteness"].get<double>();
    CHECK_NEAR(s * 100.0, std::round(s * 100.0), 1e-9);
}

/** --details adds the per-pair H-bonds and the full per-suite table without changing the results. */
TEST(details_adds_fields_only) {
    const json plain = run_json("", "plain.json");
    const json det = run_json("--details", "details.json");
    CHECK(det.contains("backbone_torsions"));
    CHECK(det.contains("score_details"));
    CHECK(det["pairs"].at(0).contains("hbonds"));
    CHECK_EQ(det["pairs_score"], plain["pairs_score"]);
    CHECK_EQ(det["backbone_suiteness"], plain["backbone_suiteness"]);
    CHECK_EQ(det["tier_summary"], plain["tier_summary"]);
    CHECK_EQ(det["pairs"].size(), plain["pairs"].size());
}

/** --no-score reports the annotation without scores. */
TEST(no_score) {
    const json d = run_json("--no-score", "noscore.json");
    CHECK(d.contains("pairs"));
    CHECK(!d.contains("pairs_score"));
    CHECK(!d["pairs"].at(0).contains("score"));
}

/** Without --out the report goes to standard output; --pml FILE writes a PyMOL script. */
TEST(stdout_and_pml) {
    const auto pml = rubric_test::write_temp("view.pml", "");
    const auto r = rubric(q(k1GID) + " --no-download --pml " + q(pml));
    CHECK_EQ(r.code, 0);
    CHECK(r.out.find("\"pairs_score\"") != std::string::npos);
    const std::string script = rubric_test::read_text(pml);
    CHECK(script.find("cmd.load(") != std::string::npos);
    CHECK(script.find("select ") != std::string::npos);
    CHECK(script.find("rubric_priority") != std::string::npos);
}

/** Missing inputs, unknown options and an unavailable PDB ID fail with exit code 1 and a message. */
TEST(error_exits) {
    const auto missing = rubric(q(fs::path("/nonexistent/missing.cif")) + " --no-download");
    CHECK_EQ(missing.code, 1);
    CHECK(missing.out.find("error") != std::string::npos);
    const auto no_download = rubric("9ZZZ --no-download --cache-dir " + q(rubric_test::write_temp("empty_cache.txt", "").parent_path() / "no_cache"));
    CHECK_EQ(no_download.code, 1);
    const auto bogus = rubric("--not-an-option");
    CHECK_EQ(bogus.code, 1);
}

/** Every developer dump-* command runs on a real structure. */
TEST(dump_commands_run) {
    for (const char* cmd : {"dump-parse", "dump-frames", "dump-hbonds", "dump-candidates", "dump-lwclass",
                            "dump-classified", "dump-scored", "dump-pairs", "dump-selection", "dump-scores",
                            "dump-structure", "dump-torsions", "dump-pucker"}) {
        const auto r = rubric(std::string(cmd) + " " + q(k1GID));
        CHECK_MSG(r.code == 0, cmd << ": " << r.out.substr(0, 200));
        CHECK_MSG(!r.out.empty(), cmd);
    }
}

/** Rebuilding the four precomputed resource caches from their sources reproduces the
 *  bundled caches byte for byte (so the shipped caches are current and the rebuild is
 *  deterministic). */
TEST(resource_caches_rebuild_identically) {
    const fs::path out = fs::temp_directory_path() / "rubric_tests" / "caches";
    fs::remove_all(out);
    for (const char* sub : {"basepair-idealized", "config", "reference"}) fs::create_directories(out / sub);
    const auto r = rubric("dump-caches " + q(out));
    CHECK_MSG(r.code == 0, r.out);
    const fs::path bundled = fs::path(RUBRIC_TESTS_DIR).parent_path() / "resources";
    for (const char* file : {"basepair-idealized/hbond_patterns_cache.json", "config/base_typing_cache.json",
                             "config/hbond_chemistry_cache.json", "reference/scoring_tables_cache.bin"}) {
        CHECK_MSG(fs::exists(out / file), file);
        CHECK_MSG(rubric_test::read_text(out / file) == rubric_test::read_text(bundled / file),
                  file << " differs from the bundled cache");
    }
}

RUBRIC_TEST_MAIN()
