/**
 * @file classify_pair.hpp
 * @brief Full LW classification of a pair: geometric classifier + template RMSD
 *        + ambiguity resolution.
 */
#ifndef PAIRFINDER_ALGORITHMS_CLASSIFY_PAIR_HPP
#define PAIRFINDER_ALGORITHMS_CLASSIFY_PAIR_HPP

#include <optional>
#include <string>
#include <vector>

#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/quality_scorer.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms::classification {

/// Final classification of a pair.
struct ClassifiedPair {
    std::string lw_class;                 ///< resolved class (or "X|Y" if ambiguous)
    bool swapped = false;                 ///< from classify_lw_class
    std::optional<double> template_rmsd;  ///< primary class RMSD (nullopt if no template)
    bool ambiguous = false;               ///< true when pipe form kept
    std::vector<std::string> ambig_classes;        ///< two classes when ambiguous
    std::vector<std::optional<double>> ambig_rmsds; ///< their RMSDs when ambiguous
};

/// Classify a pair end-to-end. ``voters`` = face-voter H-bonds (for the geometric
/// classifier); ``strict`` = selected H-bonds (for ambiguity resolution).
ClassifiedPair classify_pair(const core::ReferenceFrame& frame1,
                             const core::ReferenceFrame& frame2, const core::Residue& res1,
                             const core::Residue& res2,
                             const std::vector<hbond::HBond>& voters,
                             const std::vector<hbond::HBond>& strict,
                             const TemplateAligner& aligner, const core::BaseTyping& typing);

/// Selection-ranking quality score (candidate.quality_score) for a classified
/// pair: QualityScorer with the ambiguous-max + sugar/HH geometry fallbacks.
double compute_selection_quality(const ValidationResult& validation,
                                 const std::string& sequence,
                                 const std::vector<hbond::HBond>& strict,
                                 const ClassifiedPair& classified, const QualityScorer& qscorer,
                                 const HBondPatterns& patterns, const core::BaseTyping& typing);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_CLASSIFY_PAIR_HPP
