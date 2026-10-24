/**
 * @file carbonyl_margin.hpp
 * @brief Carbonyl-margin Watson/Sugar correction for pyrimidine O2.
 *
 * A pyrimidine's O2 is shared between its Watson and Sugar edges; which edge is
 * engaged is decided by the DIRECTION the partner's H-bond enters O2 (toward N3 =
 * Watson). ``o2_watson_angle`` returns the signed in-plane angle (+ toward N3);
 * ``correct_cis_o2`` flips a cis pyrimidine S->W when that angle is decisively
 * Watson, guarded by partner-edge and committed-W-face checks.
 */
#ifndef PAIRFINDER_ALGORITHMS_CARBONYL_MARGIN_HPP
#define PAIRFINDER_ALGORITHMS_CARBONYL_MARGIN_HPP

#include <optional>
#include <string>
#include <unordered_set>
#include <utility>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms::classification {

constexpr double kDecisiveWatsonDeg = 20.0;

/// Signed in-plane angle (deg, + toward N3/Watson) of the closest partner base
/// atom into res_b's O2, or nullopt if no contact / geometry unavailable.
std::optional<double> o2_watson_angle(const core::Residue& res_b,
                                      const core::Residue& partner_res,
                                      const core::BaseTyping& typing);

/// Result of the guarded cis O2 S->W correction.
struct O2Correction {
    bool fired = false;
    std::string new_class;
    bool swapped = false;
};

/// Apply the guarded cis pyrimidine-O2 S->W correction (returns fired=false if
/// no change). committed_w_faces = residue ids whose W edge is already used.
O2Correction correct_cis_o2(const std::string& lw_class, const std::string& res_id1,
                            const std::string& res_id2, const core::Residue& res1,
                            const core::Residue& res2,
                            const std::unordered_set<std::string>& committed_w_faces,
                            const core::BaseTyping& typing);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_CARBONYL_MARGIN_HPP
