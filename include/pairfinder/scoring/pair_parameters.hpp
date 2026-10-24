/**
 * @file pair_parameters.hpp
 * @brief Tsukuba/X3DNA six base-pair parameters from two frames.
 */
#ifndef PAIRFINDER_SCORING_PAIR_PARAMETERS_HPP
#define PAIRFINDER_SCORING_PAIR_PARAMETERS_HPP

#include <pairfinder/core/reference_frame.hpp>

namespace pairfinder::scoring {

struct PairParameters {
    double shear = 0.0, stretch = 0.0, stagger = 0.0;
    double buckle = 0.0, propeller = 0.0, opening = 0.0;
};

/// Six base-pair parameters (convention: frame2 first, frame1 second).
PairParameters compute_pair_parameters(const core::ReferenceFrame& frame1,
                                       const core::ReferenceFrame& frame2);

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_PAIR_PARAMETERS_HPP
