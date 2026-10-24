/**
 * @file confidence.hpp
 * @brief Intrinsic LW-class confidence.
 *
 *   confidence = sqrt(template_fit_score * edge_pattern_match_score)
 */
#ifndef PAIRFINDER_ALGORITHMS_CONFIDENCE_HPP
#define PAIRFINDER_ALGORITHMS_CONFIDENCE_HPP

#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms::classification {

struct ConfidenceResult {
    std::string display_label;          ///< single class or "ambiguous (X|Y)"
    std::vector<double> confidences;    ///< one per named class (rounded to 3 dp)
    bool is_ambiguous = false;
};

/// Intrinsic confidence for the classified pair (res1/res2 in post-swap order,
/// primary_class = final lw_class, hbonds = strict selected H-bonds).
ConfidenceResult compute_classification_confidence(
    const core::Residue& res1, const core::Residue& res2, const std::string& sequence,
    const std::string& primary_class, const std::vector<hbond::HBond>& hbonds,
    const TemplateAligner& aligner, const HBondPatterns& patterns,
    const core::BaseTyping& typing);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_CONFIDENCE_HPP
