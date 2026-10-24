/**
 * @file saenger.hpp
 * @brief Saenger 28-type H-bond pattern matching, used for ambiguity resolution.
 */
#ifndef PAIRFINDER_ALGORITHMS_SAENGER_HPP
#define PAIRFINDER_ALGORITHMS_SAENGER_HPP

#include <vector>

#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::algorithms::classification {

/// Result of a Saenger match: the per-residue face letters (W/H/S).
struct SaengerMatch {
    bool found = false;
    char face1 = 0;
    char face2 = 0;
};

/// Match the pair's H-bonds against the 28 Saenger types. Returns the implied per-residue faces, or found=false.
SaengerMatch match_saenger(const core::Residue& res1, const core::Residue& res2,
                           const std::vector<hbond::HBond>& hbonds,
                           const core::BaseTyping& typing);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_SAENGER_HPP
