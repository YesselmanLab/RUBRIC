/**
 * @file score_candidates.hpp
 * @brief Full per-candidate scoring pass.
 *
 * For every valid candidate: find H-bonds, assign pair_category, run the LW
 * pipeline for base-base pairs (classify + template RMSD + ambiguity), compute
 * the selection quality_score / num_hbonds / has_strong_base_hbond / confidence,
 * then the carbonyl-O2 correction as a 2-pass post-step. Produces the scored
 * candidate set that selection consumes.
 */
#ifndef PAIRFINDER_ALGORITHMS_SCORE_CANDIDATES_HPP
#define PAIRFINDER_ALGORITHMS_SCORE_CANDIDATES_HPP

#include <string>
#include <vector>

#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/quality_scorer.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/structure.hpp>

namespace pairfinder::algorithms::classification {

/// A scored candidate pair: the fields downstream selection and output use.
struct ScoredCandidate {
    std::string res_id1;
    std::string res_id2;
    std::string res_name1;  ///< pipeline single-letter (DC->C)
    std::string res_name2;
    const core::ReferenceFrame* frame1 = nullptr;
    const core::ReferenceFrame* frame2 = nullptr;
    ValidationResult validation;
    std::string pair_category = "unknown";
    std::string lw_class;                  ///< "" if not base-base; may be "X|Y"
    std::optional<double> template_rmsd;
    double quality_score = 0.0;
    int num_hbonds = 0;
    std::vector<std::string> hbond_categories;  ///< per strict hbond, pre-swap order
    bool has_strong_base_hbond = false;
    bool is_ambiguous = false;
    std::string lw_class_display;
    std::vector<double> lw_class_confidence;
    std::string precorr_lw;                ///< set when O2 correction fired
};

/// Run the full scoring pass over all candidates (frames keyed by res_id;
/// residues taken from ``structure`` with pipeline base types already applied).
std::vector<ScoredCandidate> score_candidates(
    const core::Structure& structure,
    const std::vector<std::pair<std::string, core::ReferenceFrame>>& frames,
    const core::BaseTyping& typing, hbond::HBondChemistry& chem,
    const TemplateAligner& aligner, const HBondPatterns& patterns,
    const QualityScorer& qscorer, double max_hbond_distance = 3.5);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_SCORE_CANDIDATES_HPP
