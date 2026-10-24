/**
 * @file edge_classifier.hpp
 * @brief Leontis-Westhof edge classification (geometry + H-bond face voting).
 *
 * The classification path: cis/trans orientation, centroid face
 * detection, atom-confirmation, cross-check, H-bond face signal (+ boundary
 * ambiguity), close-contact refinement, and the pyrimidine W/S tiebreaker.
 * Template-RMSD disambiguation, Saenger fallback, and confidence are later steps.
 */
#ifndef PAIRFINDER_ALGORITHMS_EDGE_CLASSIFIER_HPP
#define PAIRFINDER_ALGORITHMS_EDGE_CLASSIFIER_HPP

#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms::classification {

/// Canonical LW class (e.g. "cWW", "tHS", or ambiguous "cWH|cWW") plus a swap
/// flag - when true the caller must swap res1/res2 (and parallel fields).
struct LwResult {
    std::string lw_class;
    bool swapped = false;
};

/// Classify the LW edge class for a residue pair from frames + face-voter H-bonds.
LwResult classify_lw_class(const core::ReferenceFrame& frame1,
                           const core::ReferenceFrame& frame2,
                           const core::Residue& res1, const core::Residue& res2,
                           const std::vector<hbond::HBond>& hbonds,
                           const core::BaseTyping& typing);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_EDGE_CLASSIFIER_HPP
