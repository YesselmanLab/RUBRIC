/**
 * @file base_frame_calculator.hpp
 * @brief Compute each residue's base reference frame.
 *
 * As in 3DNA: least-squares fit of the
 * per-base standard template ring atoms (Atomic_<X>.pdb) onto the residue's
 * observed ring atoms. Each base (A,G,C,U,T,I) has its own standard ring
 * geometry, so templates are loaded per base.
 */
#ifndef PAIRFINDER_ALGORITHMS_BASE_FRAME_CALCULATOR_HPP
#define PAIRFINDER_ALGORITHMS_BASE_FRAME_CALCULATOR_HPP

#include <filesystem>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pairfinder/core/atom_type.hpp>
#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/core/structure.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms {

/// Loads and caches the standard base templates (Atomic_<X>.pdb / Atomic.<x>.pdb)
/// and exposes their ring-atom coordinates (in the standard reference frame).
/// The per-residue template is selected by the NucleotideRegistry, as in 3DNA:
/// uppercase template for standard bases, lowercase for modified.
class BaseTemplateLibrary {
public:
    BaseTemplateLibrary(std::filesystem::path templates_dir,
                        const core::NucleotideRegistry& registry)
        : templates_dir_(std::move(templates_dir)), registry_(&registry) {}

    /// Ring-atom coordinates of the standard template for ``res_name``,
    /// keyed by AtomType. nullptr if the residue/template is unavailable.
    const std::unordered_map<core::AtomType, geometry::Vector3d>*
    ring_atoms(const std::string& res_name) const;

private:
    std::filesystem::path templates_dir_;
    const core::NucleotideRegistry* registry_;
    // Cached ring geometry keyed by template filename (shared across residue
    // names that resolve to the same template).
    mutable std::unordered_map<std::string,
                               std::unordered_map<core::AtomType, geometry::Vector3d>>
        cache_;
    mutable std::unordered_map<std::string, bool> missing_;
};

/// Base frame for one residue using ``lib``. ``valid`` false if no template or
/// fewer than 3 shared ring atoms.
core::ReferenceFrame calculate_base_frame(const core::Residue& residue,
                                          const BaseTemplateLibrary& lib);

/// Frames for every residue with a valid frame, keyed by res_id (in order).
std::vector<std::pair<std::string, core::ReferenceFrame>>
calculate_all_frames(const core::Structure& structure,
                     const BaseTemplateLibrary& lib);

}  // namespace pairfinder::algorithms

#endif  // PAIRFINDER_ALGORITHMS_BASE_FRAME_CALCULATOR_HPP
