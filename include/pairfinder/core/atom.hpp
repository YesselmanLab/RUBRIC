/**
 * @file atom.hpp
 * @brief A single atom: name, type, element, coordinates.
 *
 * Value type.
 */
#ifndef PAIRFINDER_CORE_ATOM_HPP
#define PAIRFINDER_CORE_ATOM_HPP

#include <string>

#include <pairfinder/core/atom_type.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::core {

struct Atom {
    std::string name;                         ///< Raw PDB atom name (e.g. "O2'").
    AtomType type = AtomType::UNKNOWN;         ///< Typed name for O(1) lookups.
    std::string element;                       ///< Element symbol (may be empty).
    geometry::Vector3d coords;                 ///< Cartesian coordinates (Angstrom).
    char alt_loc = ' ';                        ///< Alternate-location indicator.

    Atom() = default;
    Atom(std::string name_, geometry::Vector3d coords_, std::string element_ = "",
         char alt_loc_ = ' ')
        : name(name_),
          type(atom_type_from_name(name_)),
          element(std::move(element_)),
          coords(coords_),
          alt_loc(alt_loc_ == '\0' ? ' ' : alt_loc_) {}

    double distance_to(const Atom& other) const {
        return coords.distance_to(other.coords);
    }
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_ATOM_HPP
