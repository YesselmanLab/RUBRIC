/**
 * @file cli_support.hpp
 * @brief Resource-path resolution shared by the CLI and the debug commands.
 *
 * Each path defaults to the bundled resources/ tree (resolved exe-relative) and
 * is overridable by its environment variable for relocated/advanced setups.
 */
#ifndef PAIRFINDER_APP_CLI_SUPPORT_HPP
#define PAIRFINDER_APP_CLI_SUPPORT_HPP

#include <string>

#include <pairfinder/pipeline.hpp>

namespace pfcli {

std::string templates_dir();        ///< Atomic_<X>.pdb base templates
std::string nucleotides_config();   ///< modified_nucleotides.json
std::string ligand_db();            ///< ligand_hbond_db.json
std::string idealized_dir();        ///< basepair-idealized templates
std::string exemplar_dir();         ///< basepair-catalog-exemplars
std::string reference_dir();        ///< empirical-scorer reference data dir
std::string richardson_suites();    ///< richardson_suites.json

/// Pipeline resource paths with each env-var override applied.
pairfinder::PipelineConfig cli_pipeline_config();

}  // namespace pfcli

#endif  // PAIRFINDER_APP_CLI_SUPPORT_HPP
