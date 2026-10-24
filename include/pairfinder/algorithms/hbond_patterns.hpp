/**
 * @file hbond_patterns.hpp
 * @brief Expected LW H-bond patterns (template REMARK 7 + curated dict) and
 *        matching, for the intrinsic confidence H-bond term.
 *
 * Prefers patterns parsed
 * from the idealized templates' REMARK 7 records (geometric fallback otherwise),
 * falling back to the hand-curated pattern table.
 */
#ifndef PAIRFINDER_ALGORITHMS_HBOND_PATTERNS_HPP
#define PAIRFINDER_ALGORITHMS_HBOND_PATTERNS_HPP

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/core/base_typing.hpp>

namespace pairfinder::algorithms::classification {

using HBondPattern = std::pair<std::string, std::string>;  // (donor_atom, acceptor_atom)

/// Loads/caches the template-derived pattern table and exposes the expected
/// donor-acceptor atom pairs for a (class, sequence).
class HBondPatterns {
public:
    /// @param idealized_dir basepair-idealized root; @param chem for geometric fallback.
    /// @param use_cache when true, load the precomputed ``hbond_patterns_cache.json``
    ///   from @p idealized_dir if present (an exact serialization of the scanned
    ///   table), falling back to scanning the .pdb files when it is missing or
    ///   unreadable. Pass false to force a fresh scan (used to regenerate the cache).
    HBondPatterns(std::filesystem::path idealized_dir, const hbond::HBondChemistry& chem,
                  bool use_cache = true);

    /// Expected H-bonds for (lw_class, sequence): template patterns first, then
    /// the curated dict. Sequence is normalized to its 2-letter parent form.
    std::vector<HBondPattern> get_expected_hbonds(const std::string& lw_class,
                                                  const std::string& sequence,
                                                  const core::BaseTyping& typing) const;

    /// Serialize the loaded pattern table to JSON, so it can be reloaded without
    /// re-parsing every idealized .pdb. Written by the ``dump-hbond-patterns``
    /// tool; consumed by the ``use_cache`` constructor path.
    void save_cache(const std::filesystem::path& path) const;

private:
    /// Populate ``template_table_`` from a cache file. Returns false (leaving the
    /// table untouched) if the file is absent or malformed, so the caller scans.
    bool load_cache(const std::filesystem::path& path);

    std::unordered_map<std::string, std::unordered_map<std::string, std::vector<HBondPattern>>>
        template_table_;
};

/// True if ``found`` (donor,acceptor) matches any expected pattern (bidirectional).
bool is_hbond_match(const HBondPattern& found, const std::vector<HBondPattern>& expected);

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_HBOND_PATTERNS_HPP
