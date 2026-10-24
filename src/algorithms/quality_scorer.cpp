/**
 * @file quality_scorer.cpp
 * @brief Composite 0-1 quality score for selection ranking.
 */
#include <pairfinder/algorithms/quality_scorer.hpp>

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace pairfinder::algorithms::classification {

namespace {

constexpr double kIdealDistMin = 2.7;
constexpr double kIdealDistMax = 3.2;
constexpr double kCwwRmsdPerfect = 0.3, kCwwRmsdCeiling = 1.25;
constexpr double kNoncwwRmsdPerfect = 0.5, kNoncwwRmsdCeiling = 2.5;
constexpr double kCwwGeomBoost = 0.85, kNoncwwGeomBoost = 1.8;
constexpr double kCwwLenFull = 0.5, kCwwLenZero = 0.8;
constexpr double kNoncwwLenFull = 1.0, kNoncwwLenZero = 2.0;

constexpr int kNoncwwDefaultHbondCount = 1;

int expected_count_for(const std::string& lw, const std::string& seq) {
    static const std::unordered_map<std::string, int> by_class = {
        {"tWW", 2}, {"cWH", 2}, {"tWH", 2}};
    static const std::unordered_map<std::string, int> by_seq = {
        {"GC", 3}, {"CG", 3}, {"AU", 2}, {"UA", 2}, {"GU", 2}, {"UG", 2}};
    if (!lw.empty() && lw != "cWW") {
        auto it = by_class.find(lw);
        return it == by_class.end() ? kNoncwwDefaultHbondCount : it->second;
    }
    auto it = by_seq.find(seq);
    return it == by_seq.end() ? 2 : it->second;
}

double round3(double v) { return std::round(v * 1000.0) / 1000.0; }

double rmsd_score(double rmsd, bool non_cww) {
    const double perfect = non_cww ? kNoncwwRmsdPerfect : kCwwRmsdPerfect;
    const double ceiling = non_cww ? kNoncwwRmsdCeiling : kCwwRmsdCeiling;
    if (rmsd <= perfect) return 1.0;
    if (rmsd >= ceiling) return 0.0;
    return 1.0 - (rmsd - perfect) / (ceiling - perfect);
}

double geometry_implied_coverage(double rmsd, bool non_cww) {
    if (non_cww) {
        if (rmsd <= 0.8) return 0.9;
        if (rmsd <= 1.2) return 0.75;
        if (rmsd <= 1.8) return 0.5;
    } else {
        if (rmsd <= 0.3) return 0.9;
        if (rmsd <= 0.5) return 0.75;
        if (rmsd <= 0.85) return 0.5;
    }
    return 0.0;
}

double geometry_leniency(double rmsd, bool non_cww) {
    const double full = non_cww ? kNoncwwLenFull : kCwwLenFull;
    const double zero = non_cww ? kNoncwwLenZero : kCwwLenZero;
    if (rmsd <= full) return 1.0;
    if (rmsd >= zero) return 0.0;
    return 1.0 - (rmsd - full) / (zero - full);
}

double distance_score(double dist, double leniency) {
    if (kIdealDistMin <= dist && dist <= kIdealDistMax) return 1.0;
    if (dist < kIdealDistMin) return std::max(0.5, 1.0 - (kIdealDistMin - dist) / 0.5);
    const double lenient_max = kIdealDistMax + 1.0 * leniency;
    if (dist <= lenient_max) return 1.0;
    return std::max(0.0, 1.0 - (dist - lenient_max) / 0.5);
}

double alignment_score(double alignment) {
    if (alignment <= 1.0) return 1.0;
    if (alignment >= 2.0) return 0.0;
    return 1.0 - (alignment - 1.0);
}

double pattern_coverage(const std::vector<const ScoringHBond*>& found,
                        const std::vector<HBondPattern>& expected) {
    int matched = 0;
    for (const auto& pat : expected) {
        for (const ScoringHBond* hb : found) {
            if (is_hbond_match({hb->donor_atom, hb->acceptor_atom}, {pat})) {
                ++matched;
                break;
            }
        }
    }
    return static_cast<double>(matched) / expected.size();
}

double count_coverage(int found_count, int expected_count) {
    if (expected_count == 0) return 0.0;
    return std::min(static_cast<double>(found_count) / expected_count, 1.0);
}

double hbond_quality(const std::vector<const ScoringHBond*>& hbonds, double rmsd, bool non_cww) {
    if (hbonds.empty()) return 0.0;
    const double leniency = geometry_leniency(rmsd, non_cww);
    double sum = 0.0;
    for (const ScoringHBond* hb : hbonds) {
        const double ds = distance_score(hb->distance, leniency);
        const double as = alignment_score(hb->alignment);
        sum += 0.7 * ds + 0.3 * as;
    }
    return sum / hbonds.size();
}

}  // namespace

double QualityScorer::compute_bp_score(const std::string& sequence, double rmsd,
                                       const std::vector<const ScoringHBond*>& found,
                                       const std::string& lw_class,
                                       const HBondPatterns& patterns,
                                       const core::BaseTyping& typing) const {
    const int expected_count = expected_count_for(lw_class, sequence);
    const bool non_cww = !lw_class.empty() && lw_class != "cWW";
    double rscore = rmsd_score(rmsd, non_cww);
    double coverage;
    {
        std::vector<HBondPattern> expected =
            lw_class.empty() ? std::vector<HBondPattern>{}
                             : patterns.get_expected_hbonds(lw_class, sequence, typing);
        if (!expected.empty())
            coverage = pattern_coverage(found, expected);
        else
            coverage = count_coverage(static_cast<int>(found.size()), expected_count);
    }
    double qscore = hbond_quality(found, rmsd, non_cww);

    const double boost_cutoff = non_cww ? kNoncwwGeomBoost : kCwwGeomBoost;
    if (rmsd <= boost_cutoff && !found.empty()) {
        coverage = std::max(coverage, geometry_implied_coverage(rmsd, non_cww));
        qscore = std::max(qscore, rscore * 0.8);
    }
    const double total =
        rmsd_weight_ * rscore + coverage_weight_ * coverage + quality_weight_ * qscore;
    return round3(total);
}

double QualityScorer::compute_score(const ValidationResult& validation,
                                    const std::string& sequence,
                                    const std::vector<ScoringHBond>& hbonds,
                                    std::optional<double> rmsd, const std::string& lw_class,
                                    const HBondPatterns& patterns,
                                    const core::BaseTyping& typing) const {
    if (!validation.is_valid) return 0.0;
    const double r = rmsd ? *rmsd : validation.quality_score / 10.0;
    std::vector<const ScoringHBond*> base_base;
    for (const auto& hb : hbonds)
        if (hb.base_base) base_base.push_back(&hb);
    return compute_bp_score(sequence, r, base_base, lw_class, patterns, typing);
}

}  // namespace pairfinder::algorithms::classification
