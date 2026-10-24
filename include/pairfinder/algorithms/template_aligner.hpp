/**
 * @file template_aligner.hpp
 * @brief Template-RMSD lookup + Kabsch alignment.
 */
#ifndef PAIRFINDER_ALGORITHMS_TEMPLATE_ALIGNER_HPP
#define PAIRFINDER_ALGORITHMS_TEMPLATE_ALIGNER_HPP

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms::classification {

/// Loads/caches base-pair template PDBs and aligns residue pairs to them.
class TemplateAligner {
public:
    TemplateAligner(std::filesystem::path idealized_dir, std::filesystem::path exemplar_dir)
        : idealized_dir_(std::move(idealized_dir)), exemplar_dir_(std::move(exemplar_dir)) {}

    /// Template RMSD for a pair against ``lw_class`` (forward + symmetric-reverse
    /// lookup, canonicalizing non-canonical input). nullopt if no template found.
    std::optional<double> compute_template_rmsd(const core::Residue& res1,
                                                const core::Residue& res2,
                                                const std::string& lw_class,
                                                const core::BaseTyping& typing) const;

private:
    using AtomMap = std::unordered_map<std::string, geometry::Vector3d>;
    using Template = std::pair<AtomMap, AtomMap>;  // (res1 atoms, res2 atoms)

    std::optional<std::filesystem::path> find_template(const std::string& lw_class,
                                                       const std::string& sequence) const;
    std::optional<std::filesystem::path> locate_template(const std::string& lw_class,
                                                         const std::string& sequence) const;
    const Template* load_template(const std::filesystem::path& path) const;
    double align_to_template(const core::Residue& res1, const core::Residue& res2,
                             const std::filesystem::path& path) const;

    std::filesystem::path idealized_dir_;
    std::filesystem::path exemplar_dir_;
    mutable std::unordered_map<std::string, Template> cache_;
    // find_template answers, keyed lw_class + '\x01' + sequence: the templates on
    // disk do not change during a run, so each (class, sequence) is looked up once.
    mutable std::unordered_map<std::string, std::optional<std::filesystem::path>> found_;
};

}  // namespace pairfinder::algorithms::classification

#endif  // PAIRFINDER_ALGORITHMS_TEMPLATE_ALIGNER_HPP
