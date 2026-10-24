/**
 * @file atom_type.cpp
 * @brief AtomType <-> PDB name mapping.
 */
#include <pairfinder/core/atom_type.hpp>

#include <array>
#include <unordered_map>

namespace pairfinder::core {

namespace {

// Normalize a PDB atom name: trim spaces and map the deoxy/RNA prime markers
// ("'" and "*") to a single canonical "'" form so "O2'" and "O2*" both match.
std::string canonical(std::string_view pdb_name) {
    std::string s;
    s.reserve(pdb_name.size());
    for (char c : pdb_name) {
        if (c == ' ') continue;
        s.push_back(c == '*' ? '\'' : c);
    }
    return s;
}

const std::unordered_map<std::string, AtomType>& name_table() {
    static const std::unordered_map<std::string, AtomType> table = {
        {"P", AtomType::P}, {"OP1", AtomType::OP1}, {"OP2", AtomType::OP2},
        {"OP3", AtomType::OP3}, {"O1P", AtomType::OP1}, {"O2P", AtomType::OP2},
        {"O5'", AtomType::O5_PRIME}, {"C5'", AtomType::C5_PRIME},
        {"C4'", AtomType::C4_PRIME}, {"O4'", AtomType::O4_PRIME},
        {"C3'", AtomType::C3_PRIME}, {"O3'", AtomType::O3_PRIME},
        {"C2'", AtomType::C2_PRIME}, {"O2'", AtomType::O2_PRIME},
        {"C1'", AtomType::C1_PRIME},
        {"N9", AtomType::N9}, {"C8", AtomType::C8}, {"N7", AtomType::N7},
        {"C5", AtomType::C5}, {"C6", AtomType::C6}, {"N6", AtomType::N6},
        {"O6", AtomType::O6}, {"N1", AtomType::N1}, {"C2", AtomType::C2},
        {"N2", AtomType::N2}, {"N3", AtomType::N3}, {"C4", AtomType::C4},
        {"O2", AtomType::O2}, {"N4", AtomType::N4}, {"O4", AtomType::O4},
        {"C5M", AtomType::C5M}, {"C7", AtomType::C7},
    };
    return table;
}

}  // namespace

AtomType atom_type_from_name(std::string_view pdb_name) {
    const auto& t = name_table();
    const auto it = t.find(canonical(pdb_name));
    return it == t.end() ? AtomType::UNKNOWN : it->second;
}

std::string atom_type_to_name(AtomType type) {
    switch (type) {
        case AtomType::P: return "P";
        case AtomType::OP1: return "OP1";
        case AtomType::OP2: return "OP2";
        case AtomType::OP3: return "OP3";
        case AtomType::O5_PRIME: return "O5'";
        case AtomType::C5_PRIME: return "C5'";
        case AtomType::C4_PRIME: return "C4'";
        case AtomType::O4_PRIME: return "O4'";
        case AtomType::C3_PRIME: return "C3'";
        case AtomType::O3_PRIME: return "O3'";
        case AtomType::C2_PRIME: return "C2'";
        case AtomType::O2_PRIME: return "O2'";
        case AtomType::C1_PRIME: return "C1'";
        case AtomType::N9: return "N9";
        case AtomType::C8: return "C8";
        case AtomType::N7: return "N7";
        case AtomType::C5: return "C5";
        case AtomType::C6: return "C6";
        case AtomType::N6: return "N6";
        case AtomType::O6: return "O6";
        case AtomType::N1: return "N1";
        case AtomType::C2: return "C2";
        case AtomType::N2: return "N2";
        case AtomType::N3: return "N3";
        case AtomType::C4: return "C4";
        case AtomType::O2: return "O2";
        case AtomType::N4: return "N4";
        case AtomType::O4: return "O4";
        case AtomType::C5M: return "C5M";
        case AtomType::C7: return "C7";
        case AtomType::UNKNOWN: return "UNKNOWN";
    }
    return "UNKNOWN";
}

}  // namespace pairfinder::core
