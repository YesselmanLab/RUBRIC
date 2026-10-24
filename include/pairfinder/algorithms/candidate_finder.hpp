/**
 * @file candidate_finder.hpp
 * @brief Pair-candidate finding + geometric validation.
 *
 * For every residue pair whose frame origins are within max_distance, computes
 * the frame-based geometry (origin distance, vertical separation along the mean
 * helix axis, base-plane angle, glycosidic-N distance, direction cosines) and
 * the validation flags against ValidationThresholds.
 */
#ifndef PAIRFINDER_ALGORITHMS_CANDIDATE_FINDER_HPP
#define PAIRFINDER_ALGORITHMS_CANDIDATE_FINDER_HPP

#include <string>
#include <utility>
#include <vector>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/structure.hpp>

namespace pairfinder::algorithms {

/// Default geometric validation thresholds.
struct ValidationThresholds {
    double max_dorg = 15.0;
    double max_d_v = 2.5;
    double max_plane_angle = 70.0;
    double min_dNN = 4.5;
    double d_v_weight = 1.5;
    double plane_angle_divisor = 180.0;
};

/// Result of validating one pair's geometry.
struct ValidationResult {
    double dorg = 0.0;
    double d_v = 0.0;
    double plane_angle = 0.0;
    double dNN = 0.0;
    double dir_x = 0.0;
    double dir_y = 0.0;
    double dir_z = 0.0;
    double quality_score = 0.0;
    bool distance_check = false;
    bool d_v_check = false;
    bool plane_angle_check = false;
    bool dNN_check = false;
    bool is_valid = false;
    /// Ring-centroid stacking signature (class-independent): one base sits over
    /// the other (centroid offset larger ALONG the mean base normal than in-plane).
    /// A more discriminating stacking detector than the dNN-based
    /// is_stacking_signature; catches stacks in classes that one cannot
    /// (cWS/cSS/tWS/tWH), where real S-edge pairs have low dNN by design.
    bool ring_stacked = false;
};

/// A validated pair candidate.
struct PairCandidate {
    std::string res_id1;
    std::string res_id2;
    std::string res_name1;
    std::string res_name2;
    ValidationResult validation;
};

/// Validate one pair's geometry from frames + glycosidic-N positions.
ValidationResult validate_pair(const core::ReferenceFrame& f1,
                               const core::ReferenceFrame& f2,
                               const geometry::Vector3d& n1n9_1,
                               const geometry::Vector3d& n1n9_2,
                               const ValidationThresholds& thr = {});

/// All candidates (pairs within max_distance with both glycosidic Ns present),
/// each carrying its validation result (valid and invalid alike).
std::vector<PairCandidate> find_candidates(
    const core::Structure& structure,
    const std::vector<std::pair<std::string, core::ReferenceFrame>>& frames,
    const core::BaseTyping& typing, double max_distance = 15.0,
    const ValidationThresholds& thr = {});

}  // namespace pairfinder::algorithms

#endif  // PAIRFINDER_ALGORITHMS_CANDIDATE_FINDER_HPP
