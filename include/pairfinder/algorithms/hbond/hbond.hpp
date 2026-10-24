/**
 * @file hbond.hpp
 * @brief Hydrogen-bond value type.
 */
#ifndef PAIRFINDER_ALGORITHMS_HBOND_HBOND_HPP
#define PAIRFINDER_ALGORITHMS_HBOND_HBOND_HPP

#include <string>

namespace pairfinder::algorithms::hbond {

/// A single hydrogen bond with explicit donor/acceptor roles.
struct HBond {
    std::string donor_res_id;
    std::string acceptor_res_id;
    std::string donor_atom;
    std::string acceptor_atom;
    double distance = 0.0;
    int h_slot_idx = -1;
    int lp_slot_idx = -1;
    double alignment_score = 0.0;
};

}  // namespace pairfinder::algorithms::hbond

#endif  // PAIRFINDER_ALGORITHMS_HBOND_HBOND_HPP
