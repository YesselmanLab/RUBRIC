/**
 * Base pairs: H-bond finding, Leontis-Westhof classification and pair scoring,
 * on the idealized pair templates bundled in resources/basepair-idealized and on
 * deliberately distorted copies of them.
 */
#include <filesystem>
#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/scoring/scorer.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

using rubric_test::run_engine;
namespace fs = std::filesystem;

namespace {

fs::path idealized(const std::string& lw, const std::string& seq) {
    return pairfinder::resources::basepair_idealized() / lw / (seq + ".pdb");
}

/// The two residues of a two-nucleotide file, in file order.
std::pair<const pairfinder::core::Residue*, const pairfinder::core::Residue*> two_residues(
        const pairfinder::PreparedStructure& prep) {
    std::vector<const pairfinder::core::Residue*> res;
    for (const auto& ch : prep.structure.chains())
        for (const auto& r : ch.residues()) res.push_back(&r);
    return {res.at(0), res.at(1)};
}

/// Atom-name pairs of the H-bonds found between the two residues of a template.
std::set<std::string> hbond_atoms(const fs::path& path, double max_distance = 4.0) {
    const auto prep = rubric_test::pipeline().prepare(path.string());
    auto [a, b] = two_residues(prep);
    pairfinder::algorithms::hbond::HBondFinder finder(rubric_test::pipeline().chem(), max_distance);
    std::set<std::string> out;
    for (const auto& h : finder.find_between(*a, *b)) {
        CHECK_MSG(h.distance > 2.4 && h.distance < 3.6, h.donor_atom << "-" << h.acceptor_atom << " " << h.distance);
        out.insert(h.donor_atom + "-" + h.acceptor_atom);
    }
    return out;
}

/// C1'-C1' unit vector of a template (from the first residue to the second).
pairfinder::geometry::Vector3d c1_axis(const fs::path& path) {
    const auto prep = rubric_test::pipeline().prepare(path.string());
    auto [a, b] = two_residues(prep);
    return (b->get_atom("C1'")->coords - a->get_atom("C1'")->coords).normalized();
}

}  // namespace

/** Watson-Crick G-C has its three H-bonds and A-U its two. */
TEST(watson_crick_hbonds) {
    const auto gc = hbond_atoms(idealized("cWW", "GC"));
    CHECK_EQ(gc.size(), std::size_t{3});
    CHECK(gc.count("N1-N3") == 1);  // G N1 donates to C N3
    CHECK(gc.count("N2-O2") == 1);
    CHECK(gc.count("N4-O6") == 1);
    const auto au = hbond_atoms(idealized("cWW", "AU"));
    CHECK_EQ(au.size(), std::size_t{2});
    CHECK(au.count("N6-O4") == 1);
    CHECK(au.count("N3-N1") == 1);
}

/** No H-bond is reported once the bases are pulled 3 A apart, or with a 2 A distance limit. */
TEST(hbond_distance_limits) {
    const auto path = idealized("cWW", "GC");
    CHECK(hbond_atoms(path, 2.0).empty());
    const auto axis = c1_axis(path);
    const auto far = rubric_test::write_temp(
        "gc_far.pdb", rubric_test::shift_chain(rubric_test::read_text(path), 'B',
                                               3.0 * axis.x, 3.0 * axis.y, 3.0 * axis.z));
    CHECK(hbond_atoms(far).empty());
}

/** Every idealized template that forms a pair is classified as its own LW class. */
TEST(templates_classified_as_their_class) {
    int paired = 0;
    std::set<std::string> classes_seen;
    for (const auto& dir : fs::directory_iterator(pairfinder::resources::basepair_idealized())) {
        if (!dir.is_directory()) continue;
        const std::string lw = dir.path().filename().string();
        for (const auto& f : fs::directory_iterator(dir.path())) {
            if (f.path().extension() != ".pdb") continue;
            const auto r = run_engine(f.path(), false);
            CHECK_MSG(r.selected.size() <= 1, f.path());
            if (r.selected.empty()) continue;  // base-only templates without the sugar H-bond atoms
            ++paired;
            classes_seen.insert(lw);
            CHECK_MSG(!r.selected[0].is_ambiguous, f.path());
            CHECK_MSG(r.selected[0].lw_class == lw, f.path() << " classified " << r.selected[0].lw_class);
        }
    }
    CHECK_MSG(paired >= 130, paired << " of 178 templates formed a pair");
    CHECK_EQ(classes_seen.size(), std::size_t{12});
}

/** Idealized pairs score at least 90, and every cWW template scores 100 with no issue. */
TEST(idealized_pairs_score_high) {
    for (const std::string seq : {"AA", "AC", "AG", "AU", "CA", "CC", "CG", "CU",
                                  "GA", "GC", "GG", "GU", "UA", "UC", "UG", "UU"}) {
        const auto r = run_engine(idealized("cWW", seq), false);
        CHECK_EQ(r.score.pair_scores.size(), std::size_t{1});
        const auto& p = r.score.pair_scores[0];
        CHECK_NEAR(p.score, 100.0, 1e-9);
        CHECK(p.issues.empty());
        CHECK_EQ(p.tier(), std::string("Preferred"));
    }
    for (const std::string& lw : {std::string("tWH"), std::string("cWH"), std::string("tHS")}) {
        for (const auto& f : fs::directory_iterator(pairfinder::resources::basepair_idealized() / lw)) {
            const auto r = run_engine(f.path(), false);
            for (const auto& p : r.score.pair_scores) CHECK_MSG(p.score >= 90.0, f.path() << " " << p.score);
        }
    }
}

/** Pulling the bases of a G-C pair apart: small stretches are tolerated, larger ones lower the
 *  score step by step and name the stretch and H-bond distance issues, and beyond about 0.8 A
 *  the bases no longer count as paired. */
TEST(stretching_lowers_score) {
    const auto path = idealized("cWW", "GC");
    const auto axis = c1_axis(path);
    const std::string pdb = rubric_test::read_text(path);
    auto stretched = [&](double d) {
        return run_engine(rubric_test::write_temp(
            "gc_stretch.pdb", rubric_test::shift_chain(pdb, 'B', d * axis.x, d * axis.y, d * axis.z)), false);
    };
    const auto small = stretched(0.2);
    CHECK_EQ(small.score.pair_scores.size(), std::size_t{1});
    CHECK_NEAR(small.score.pair_scores[0].score, 100.0, 1e-9);

    double previous = 100.0;
    for (double d : {0.45, 0.55, 0.65}) {
        const auto r = stretched(d);
        CHECK_EQ(r.score.pair_scores.size(), std::size_t{1});
        const auto& p = r.score.pair_scores[0];
        CHECK_MSG(p.score < previous, d << " A: " << p.score << " >= " << previous);
        previous = p.score;
        std::set<std::string> issues;
        for (const auto& i : p.issues) issues.insert(i.issue);
        CHECK_MSG(issues.count("stretch") == 1 && issues.count("distance") == 1, d << " A");
    }
    CHECK(previous < 95.0);
    CHECK(stretched(1.0).score.pair_scores.empty());
}

/** Pair-score bookkeeping on real data: bounds, penalties, issue order, tiers and the structure mean. */
TEST(score_bookkeeping) {
    const auto r = run_engine(rubric_test::tests_dir() / "data" / "1GID.cif.gz", false);
    const auto& s = r.score;
    CHECK(!s.pair_scores.empty());
    double sum = 0.0;
    for (const auto& p : s.pair_scores) {
        CHECK(p.score >= 0.0 && p.score <= 100.0);
        CHECK_NEAR(p.score, std::max(0.0, 100.0 - p.penalty), 0.01);
        double weights = 0.0;
        for (std::size_t i = 0; i < p.issues.size(); ++i) {
            const auto& is = p.issues[i];
            CHECK(is.severity >= 0.0 && is.severity <= 1.0);
            CHECK(is.weight >= 0.0);
            if (i > 0) CHECK(p.issues[i - 1].weight >= is.weight);
            weights += is.weight;
        }
        CHECK_NEAR(p.penalty, weights, 0.01);
        sum += p.score;
    }
    CHECK_NEAR(s.pairs_score, sum / s.pair_scores.size(), 0.01);
    const auto t = s.tier_summary();
    CHECK_EQ(static_cast<std::size_t>(t.preferred + t.acceptable + t.review), s.pair_scores.size());
}

/** Tier boundaries: 100 Preferred, 80 to below 100 Secondary review, below 80 Priority review. */
TEST(tier_boundaries) {
    using pairfinder::scoring::pair_tier;
    CHECK_EQ(pair_tier(100.0), std::string("Preferred"));
    CHECK_EQ(pair_tier(97.5), std::string("Secondary review"));
    CHECK_EQ(pair_tier(84.0), std::string("Secondary review"));
    CHECK_EQ(pair_tier(72.0), std::string("Priority review"));
    CHECK_EQ(pair_tier(0.0), std::string("Priority review"));
}

RUBRIC_TEST_MAIN()
