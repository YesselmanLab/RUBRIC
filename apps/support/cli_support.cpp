/**
 * @file cli_support.cpp
 * @brief Resource-path resolution (see cli_support.hpp).
 */
#include "support/cli_support.hpp"

#include <cstdlib>
#include <filesystem>

#include <pairfinder/core/resource_locator.hpp>

namespace pfcli {
namespace {

namespace resources = pairfinder::resources;

// Env-var override, else the bundled-resources fallback.
std::string env_or(const char* var, const std::filesystem::path& fallback) {
    if (const char* env = std::getenv(var)) return env;
    return fallback.string();
}

}  // namespace

std::string templates_dir() {
    return env_or("PAIRFINDER_TEMPLATES_DIR", resources::templates_dir());
}
std::string nucleotides_config() {
    return env_or("PAIRFINDER_NT_CONFIG", resources::config("modified_nucleotides.json"));
}
std::string ligand_db() {
    return env_or("PAIRFINDER_LIGAND_DB", resources::config("ligand_hbond_db.json"));
}
std::string idealized_dir() {
    return env_or("PAIRFINDER_IDEALIZED_DIR", resources::basepair_idealized());
}
std::string exemplar_dir() {
    return env_or("PAIRFINDER_EXEMPLAR_DIR", resources::basepair_exemplars());
}
std::string reference_dir() {
    return env_or("PAIRFINDER_REFERENCE_DIR", resources::reference_dir());
}
std::string richardson_suites() {
    return env_or("PAIRFINDER_RICHARDSON_SUITES", resources::file("richardson_suites.json"));
}

pairfinder::PipelineConfig cli_pipeline_config() {
    pairfinder::PipelineConfig c;
    c.templates_dir = templates_dir();
    c.nucleotides_config = nucleotides_config();
    c.ligand_db = ligand_db();
    c.idealized_dir = idealized_dir();
    c.exemplar_dir = exemplar_dir();
    c.reference_dir = reference_dir();
    c.richardson_suites = richardson_suites();
    return c;
}

}  // namespace pfcli
