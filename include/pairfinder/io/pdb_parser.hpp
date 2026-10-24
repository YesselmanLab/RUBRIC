/**
 * @file pdb_parser.hpp
 * @brief Parse a structure file (PDB / mmCIF / .gz) into a Structure.
 *
 * GEMMI-backed loader (first model, alt-loc ' '/'A'/'1', modified-nucleotide
 * HETATM filter, atom-name normalization OP1->O1P etc.). Format auto-detected
 * by extension. res_id = "<chain>-<resname>-<seq+inscode>", base_type = upper.
 */
#ifndef PAIRFINDER_IO_PDB_PARSER_HPP
#define PAIRFINDER_IO_PDB_PARSER_HPP

#include <filesystem>

#include <pairfinder/core/structure.hpp>

namespace pairfinder::io {

/// Parse a PDB file into a Structure (residues grouped into chains in
/// first-seen order). Throws std::runtime_error if the file cannot be opened.
core::Structure parse_pdb(const std::filesystem::path& path);

}  // namespace pairfinder::io

#endif  // PAIRFINDER_IO_PDB_PARSER_HPP
