/**
 * @file resource_locator.hpp
 * @brief Locate the bundled `resources/` directory (templates, configs, scorer
 *        distributions, idealized templates, richardson suites).
 *
 * Resolution order (computed once):
 *   1. $RUBRIC_RESOURCES (explicit override)
 *   2. executable-relative: <exe>/resources, <exe>/../resources,
 *      <exe>/../share/rubric/resources   (deployed/installed layout)
 *   3. cwd-relative: resources, ../resources, ../../resources
 *   4. compile-time source path PAIRFINDER_SOURCE_RESOURCES (dev fallback)
 *
 * This makes the tool self-contained: it carries its own data. Follows x3dna's
 * config::ResourceLocator.
 */
#ifndef PAIRFINDER_CORE_RESOURCE_LOCATOR_HPP
#define PAIRFINDER_CORE_RESOURCE_LOCATOR_HPP

#include <filesystem>
#include <string>

namespace pairfinder::resources {

/// Resolved resources root (empty if not found). Computed once, then cached.
const std::filesystem::path& root();

/// resources/<relative-path>.
std::filesystem::path file(const std::string& relative);

/// resources/templates
std::filesystem::path templates_dir();
/// resources/config/<name>
std::filesystem::path config(const std::string& name);
/// resources/reference
std::filesystem::path reference_dir();
/// resources/basepair-idealized
std::filesystem::path basepair_idealized();
/// resources/basepair-catalog-exemplars
std::filesystem::path basepair_exemplars();

}  // namespace pairfinder::resources

#endif  // PAIRFINDER_CORE_RESOURCE_LOCATOR_HPP
