/**
 * @file base_typing.hpp
 * @brief Modified-base -> parent normalization and glycosidic-N selection.
 *
 * Provides ``normalize_base_type`` / ``is_purine`` / ``is_pyrimidine`` /
 * ``get_glycosidic_n``.
 * Lookup order: standard {A,G,C,U,T,I,P} -> DNA prefix (DA->A) -> curated
 * modified-base map -> extended map built from ligand_hbond_db.json (parent
 * field, with atom-signature inference for parent=null) -> unchanged.
 *
 * Distinct from core/nucleotide_registry (modified_nucleotides.json), which
 * drives standard frame-template selection; this one is the typing used by
 * candidate finding, H-bonds and classification.
 */
#ifndef PAIRFINDER_CORE_BASE_TYPING_HPP
#define PAIRFINDER_CORE_BASE_TYPING_HPP

#include <filesystem>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace pairfinder::core {

class BaseTyping {
public:
    /// Builds the curated maps and extends from the ligand DB at @p ligand_db_path
    /// (silently skipped if unreadable).
    /// @param use_cache when true, load the precomputed ``base_typing_cache.json``
    ///   from the ligand DB's directory instead of parsing the 1.5 MB DB (an exact
    ///   serialization of ``extended_``); falls back to parsing when it is missing.
    ///   Pass false to force a fresh parse (used to regenerate the cache).
    explicit BaseTyping(const std::filesystem::path& ligand_db_path, bool use_cache = true);

    /// Normalize a residue name to its parent base (single letter for known
    /// nucleotides; original name for unrecognized ligands).
    std::string normalize(std::string_view residue_name) const;

    bool is_purine(std::string_view base_type) const;
    bool is_pyrimidine(std::string_view base_type) const;

    /// True when the residue is DNA, including modified DNA whose name carries
    /// no ``D`` prefix (2'-fluoro AF2/GF2/CFZ/UFT, BRU, ...). Data-driven from
    /// the ligand DB ``type`` field, with the literal DA/DC/DG/DT/DU names as a
    /// floor so the answer stays sane if the DB could not be read.
    bool is_dna(std::string_view residue_name) const;

    /// True when the residue name is a nucleotide the ligand DB knows: standard
    /// bases, the DNA prefix names, the curated modified map, or anything typed
    /// RNA/DNA in the DB. Used by the parser to decide whether a HETATM record
    /// is a nucleotide worth keeping.
    bool is_known_nucleotide(std::string_view residue_name) const;

    /// True for a free nucleotide/cofactor the ligand DB types LIGAND rather
    /// than RNA/DNA (ATP, GTP, UTP, CTP, ...). Such a residue can pair, but it
    /// is not a polymer residue and carries no RNA backbone suite.
    bool is_ligand(std::string_view residue_name) const;

    /// Glycosidic-N atom name (N9 purine / N1 pyrimidine / C5 pseudouridine),
    /// or "" if unknown, with the normalize fallback.
    std::string glycosidic_n_name(std::string_view base_type) const;

    /// Serialize ``extended_`` (the DB-derived name->parent map) to JSON, so it
    /// can be reloaded without re-parsing the ligand DB.
    void save_cache(const std::filesystem::path& path) const;

private:
    /// Populate ``extended_`` from a cache file; false (table untouched) if the
    /// file is absent or malformed, so the caller parses the DB instead.
    bool load_cache(const std::filesystem::path& path);

    std::unordered_map<std::string, std::string> extended_;  ///< from ligand DB
    std::unordered_set<std::string> dna_codes_;              ///< ligand DB type == "DNA"
    std::unordered_set<std::string> ligand_codes_;           ///< ligand DB type == "LIGAND"
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_BASE_TYPING_HPP
