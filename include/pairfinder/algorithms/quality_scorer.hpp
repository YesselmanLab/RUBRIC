/**
 * @file quality_scorer.hpp
 * @brief Composite 0-1 quality score for pair ranking.
 *
 * This is the SELECTION ranking metric (candidate.quality_score), distinct from
 * the user-facing empirical Scorer. Combines RMSD-to-template, H-bond
 * coverage (pattern-aware), and H-bond quality.
 */
#ifndef PAIRFINDER_ALGORITHMS_QUALITY_SCORER_HPP
#define PAIRFINDER_ALGORITHMS_QUALITY_SCORER_HPP

#include <optional>
#include <string>
#include <vector>

#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/core/base_typing.hpp>

namespace pairfinder::algorithms::classification {

/// One H-bond in the scorer's dict form (context = base_base vs other).
struct ScoringHBond {
    std::string donor_atom;
    std::string acceptor_atom;
    double distance = 3.0;
    double alignment = 1.0;
    bool base_base = false;
};

class QualityScorer {
public:
    QualityScorer(double rmsd_weight = 0.35, double coverage_weight = 0.30,
                  double quality_weight = 0.35)
        : rmsd_weight_(rmsd_weight),
          coverage_weight_(coverage_weight),
          quality_weight_(quality_weight) {}

    /// 0-1 quality score. rmsd nullopt -> validation.quality_score/10.
    double compute_score(const ValidationResult& validation, const std::string& sequence,
                         const std::vector<ScoringHBond>& hbonds,
                         std::optional<double> rmsd, const std::string& lw_class,
                         const HBondPatterns& patterns, const core::BaseTyping& typing) const;

private:
    double compute_bp_score(const std::string& sequence, double rmsd,
                            const std::vector<const ScoringHBond*>& found,
                            const std::string& lw_class, const HBondPatterns& patterns,
                            const core::BaseTyping& typing) const;

    double rmsd_weight_, coverage_weight_, quality_weight_;
};

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_QUALITY_SCORER_HPP
