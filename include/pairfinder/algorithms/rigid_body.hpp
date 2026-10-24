/**
 * @file rigid_body.hpp
 * @brief Tsukuba / 3DNA rigid-body base-pair parameters (CLI output block).
 *
 * NOTE: this is a DIFFERENT decomposition from the scorer's
 * `scoring::compute_pair_parameters` - that one
 * uses the hinge/half-angle convention; this one flips frame2 by diag(1,-1,-1),
 * halves the relative rotation, and expresses everything in the mid-frame
 * (Lu & Olson 2003 / 3DNA CEHS). They are not interchangeable.
 */
#ifndef PAIRFINDER_ALGORITHMS_RIGID_BODY_HPP
#define PAIRFINDER_ALGORITHMS_RIGID_BODY_HPP

#include <pairfinder/core/reference_frame.hpp>

namespace pairfinder::algorithms {

/// The six Tsukuba/3DNA base-pair parameters (translations in Angstroms,
/// rotations in degrees), in CLI emission order.
struct RigidBodyParameters {
    double shear = 0.0;
    double stretch = 0.0;
    double stagger = 0.0;
    double buckle = 0.0;
    double propeller = 0.0;
    double opening = 0.0;
};

RigidBodyParameters compute_rigid_body_parameters(const core::ReferenceFrame& frame1,
                                                  const core::ReferenceFrame& frame2);

}  // namespace pairfinder::algorithms

#endif  // PAIRFINDER_ALGORITHMS_RIGID_BODY_HPP
