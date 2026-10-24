/**
 * @file torsions.hpp
 * @brief Backbone torsion angles.
 *
 * Computes the six suite backbone torsions (alpha, beta, gamma, delta, epsilon,
 * zeta) per residue, sequence-ordered (chain, seq) with adjacent prev/next.
 * chi / sugar-pucker are not needed by suiteness scoring and are omitted.
 */
#ifndef PAIRFINDER_SCORING_TORSIONS_HPP
#define PAIRFINDER_SCORING_TORSIONS_HPP

#include <optional>
#include <string>
#include <unordered_map>

#include <pairfinder/core/structure.hpp>

namespace pairfinder::scoring {

struct ResidueTorsions {
    std::optional<double> alpha, beta, gamma, delta, epsilon, zeta;
};

/// Parsed res_id sequence component.
struct ResSeq {
    std::string chain;
    long num = 0;
    std::string icode;
};

/// Parse "chain-base-seq[icode]" into (chain, seq_num, insertion_code).
/// Handles negative sequence numbers ("G-C--8" -> {"G", -8, ""}) and
/// insertion / alt-conf codes ("A-G-20A" -> {"A", 20, "A"}) by taking the
/// full sequence field (everything after the second '-') and splitting a
/// signed integer prefix from trailing letters. Returns nullopt if there is
/// no integer seq.
std::optional<ResSeq> parse_res_seq(const std::string& res_id);

/// res_id -> torsions, for every residue (sequence-ordered prev/next neighbors).
std::unordered_map<std::string, ResidueTorsions> compute_all_torsions(
    const core::Structure& structure);

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_TORSIONS_HPP
