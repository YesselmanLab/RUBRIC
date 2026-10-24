/**
 * @file atom_type.hpp
 * @brief Strongly-typed atom names for O(1) typed lookups.
 *
 * Mirrors the x3dna pattern (typed AtomType enum instead of string compares).
 * Covers the RNA/DNA base, sugar and phosphate atoms used by the pair finder.
 * Prime atoms (O2', C1', ...) map to *_PRIME spellings.
 */
#ifndef PAIRFINDER_CORE_ATOM_TYPE_HPP
#define PAIRFINDER_CORE_ATOM_TYPE_HPP

#include <string>
#include <string_view>

namespace pairfinder::core {

enum class AtomType {
    UNKNOWN = 0,
    // Phosphate
    P, OP1, OP2, OP3,
    // Sugar (ribose / deoxyribose)
    O5_PRIME, C5_PRIME, C4_PRIME, O4_PRIME, C3_PRIME,
    O3_PRIME, C2_PRIME, O2_PRIME, C1_PRIME,
    // Purine ring + exocyclic
    N9, C8, N7, C5, C6, N6, O6, N1, C2, N2, N3, C4,
    // Pyrimidine exocyclic (ring atoms N1,C2,N3,C4,C5,C6 shared with purine spellings)
    O2, N4, O4, C5M, C7,
};

/// Convert a PDB atom name (e.g. "O2'", "C1*", "N9") to an AtomType.
/// Returns AtomType::UNKNOWN for anything not recognised.
AtomType atom_type_from_name(std::string_view pdb_name);

/// Canonical PDB-style name for a type (prime form, e.g. "O2'").
std::string atom_type_to_name(AtomType type);

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_ATOM_TYPE_HPP
