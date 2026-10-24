/**
 * @file selection.hpp
 * @brief Two-phase pair selection.
 *
 * Phase 1: cWW canonical/wobble pairs via helix-priority
 * (helix continuity > raw score; isolated pairs stricter). Phase 2:
 * slot-aware maximum-weighted independent set over the
 * remaining candidates. Face-aware conflicts (a residue may pair on W and H
 * faces simultaneously). Pre-filtered by the stacking-signature filter.
 */
#ifndef PAIRFINDER_ALGORITHMS_SELECTION_HPP
#define PAIRFINDER_ALGORITHMS_SELECTION_HPP

#include <vector>

#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>

namespace pairfinder::algorithms::selection {

/// Final selected pairs (subset of the scored candidates, with selection-time
/// ambiguity resolution applied to lw_class). min_score from FinderConfig.
/// Per-pair helix membership: 1 if the pair sits in a detected helix segment
/// (>=2 stacked consecutive pairs), else 0. Same index order as `pairs`. Used
/// for the sugar-pucker flag.
std::vector<char> helix_membership(
    const std::vector<classification::ScoredCandidate>& pairs, const RNAChains& chains);

std::vector<classification::ScoredCandidate> select_pairs(
    std::vector<classification::ScoredCandidate> candidates, const RNAChains& chains,
    double min_score = 0.0);

/// Per-candidate selection disposition (why a candidate is / isn't in the final
/// output). Runs the REAL select_pairs to determine the selected set, then tags
/// every input candidate with the actual production cause using the same
/// predicates select_pairs uses - so downstream tools read the true reason
/// instead of reimplementing selection.
///   selected - in the final selected output
///   invalid - validation.is_valid == false
///   stacking - removed by the stacking-signature pre-filter
///   low_score - quality_score < min_score
///   adjacent - same-chain adjacent residues (phase-2 skip)
///   ineligible - fails class_valid (no strong base H-bond / SS / HH gate)
///   competition - eligible, but displaced by a better pair (helix phase or MWIS)
struct SelectionDisposition {
    std::string res_id1;
    std::string res_id2;
    std::string lw_class;
    std::string disposition;
};

std::vector<SelectionDisposition> select_pairs_dispositions(
    std::vector<classification::ScoredCandidate> candidates, const RNAChains& chains,
    double min_score = 0.0);

}  // namespace pairfinder::algorithms::selection

#endif  // PAIRFINDER_ALGORITHMS_SELECTION_HPP
