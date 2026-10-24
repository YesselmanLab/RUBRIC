/**
 * @file rna_chains.hpp
 * @brief Spatial RNA chain connectivity.
 *
 * Builds chains by following O3'-P backbone distance (< 2.75 A) rather than
 * sequence numbers, then exposes 3'/5' neighbor and position lookups used by
 * helix detection. Chain membership and internal 5'->3' order are connectivity-
 * determined (deterministic); only the chain *discovery order* (chain index)
 * depends on iteration order - and selection is insensitive to it.
 */
#ifndef PAIRFINDER_ALGORITHMS_RNA_CHAINS_HPP
#define PAIRFINDER_ALGORITHMS_RNA_CHAINS_HPP

#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pairfinder/core/structure.hpp>

namespace pairfinder::algorithms {

class RNAChains {
public:
    /// Build from a structure's residues (RNA = has C1' or P), structure order.
    static RNAChains from_structure(const core::Structure& structure);

    /// 3' neighbor res_id, or "" if at chain end.
    const std::string& next_residue(const std::string& res_id) const;
    /// 5' neighbor res_id, or "" if at chain start.
    const std::string& prev_residue(const std::string& res_id) const;
    /// (chain_idx, position) or nullopt.
    std::optional<std::pair<int, int>> position(const std::string& res_id) const;

private:
    std::unordered_map<std::string, std::string> next_;
    std::unordered_map<std::string, std::string> prev_;
    std::unordered_map<std::string, std::pair<int, int>> pos_;
    std::string empty_;
};

}  // namespace pairfinder::algorithms

#endif  // PAIRFINDER_ALGORITHMS_RNA_CHAINS_HPP
