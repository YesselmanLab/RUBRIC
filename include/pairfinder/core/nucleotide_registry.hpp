/**
 * @file nucleotide_registry.hpp
 * @brief Residue-name -> nucleotide typing, loaded from modified_nucleotides.json.
 *
 * Follows x3dna's TypeRegistry::load_nucleotides: one authoritative
 * table (standard + DNA + ~400 modified codes) mapping a residue name to its
 * one-letter code, parent base, purine flag, and modified flag. Drives standard
 * base-template selection for the reference frame (Atomic_<X>.pdb uppercase for
 * standard, Atomic.<x>.pdb lowercase for modified).
 */
#ifndef PAIRFINDER_CORE_NUCLEOTIDE_REGISTRY_HPP
#define PAIRFINDER_CORE_NUCLEOTIDE_REGISTRY_HPP

#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

namespace pairfinder::core {

/// Typing facts for one residue name.
struct NucleotideInfo {
    char one_letter = '?';   ///< code from the table (lowercase => modified)
    char canonical = '?';    ///< parent base letter, uppercase (A/G/C/U/T/I/P)
    bool is_purine = false;
    bool is_modified = false;
};

/// Loads modified_nucleotides.json once and answers per-residue typing queries.
class NucleotideRegistry {
public:
    /// @param config_file path to modified_nucleotides.json
    explicit NucleotideRegistry(const std::filesystem::path& config_file);

    /// Lookup by residue name (as parsed, upper/stripped). nullopt if unknown.
    std::optional<NucleotideInfo> lookup(std::string_view res_name) const;

    /// Standard template filename for this residue, or "" if unknown.
    /// Uppercase Atomic_<X>.pdb for standard, Atomic.<x>.pdb for modified.
    std::string template_filename(std::string_view res_name) const;

    std::size_t size() const { return table_.size(); }

private:
    std::unordered_map<std::string, NucleotideInfo> table_;
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_NUCLEOTIDE_REGISTRY_HPP
