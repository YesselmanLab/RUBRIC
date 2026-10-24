/**
 * @file differential_commands.hpp
 * @brief Hidden `dump-*` subcommands used by the differential tests.
 *
 * Each emits a sorted TSV of one pipeline stage's output so it can be compared
 * byte-for-byte between builds. They are not
 * part of the user-facing CLI; they live apart from main.cpp to keep the entry
 * point focused on the production find/classify/score path.
 */
#ifndef PAIRFINDER_APP_DIFFERENTIAL_COMMANDS_HPP
#define PAIRFINDER_APP_DIFFERENTIAL_COMMANDS_HPP

#include <string>
#include <vector>

namespace pfcli {

int dump_parse(const std::string& path);       ///< parsed atoms
int dump_frames(const std::string& path);      ///< per-residue reference frames
int dump_hbonds(const std::string& path);      ///< H-bonds (find_between)
int dump_candidates(const std::string& path);  ///< candidate pairs + validation
int dump_lwclass(const std::string& path);     ///< LW edge classification
int dump_classified(const std::string& path);  ///< full classification + confidence
int dump_scored(const std::string& path);      ///< full per-candidate scoring pass
int dump_pairs(const std::string& path);       ///< selected pairs
int dump_selection(const std::string& path);   ///< per-candidate selection disposition
int dump_scores(const std::string& path);      ///< empirical per-pair scores
int dump_structure(const std::string& path);   ///< whole-structure score
int dump_torsions(const std::string& path);    ///< backbone suite torsions
int dump_pucker(const std::string& path);      ///< per-residue delta + Pperp
int classify_suite(const std::vector<std::string>& angles);  ///< suiteness of 7 angles

}  // namespace pfcli

#endif  // PAIRFINDER_APP_DIFFERENTIAL_COMMANDS_HPP
