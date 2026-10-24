/**
 * @file nucleotide_registry.cpp
 * @brief Load modified_nucleotides.json into the residue typing table.
 */
#include <pairfinder/core/nucleotide_registry.hpp>

#include <cctype>
#include <fstream>
#include <stdexcept>

#include <nlohmann/json.hpp>

namespace pairfinder::core {

namespace {

// Parent base "type" string -> canonical uppercase one-letter.
char canonical_from_type(const std::string& type) {
    if (type == "ADENINE") return 'A';
    if (type == "GUANINE") return 'G';
    if (type == "CYTOSINE") return 'C';
    if (type == "URACIL") return 'U';
    if (type == "THYMINE") return 'T';
    if (type == "INOSINE") return 'I';
    if (type == "PSEUDOURIDINE") return 'P';
    return '?';
}

}  // namespace

NucleotideRegistry::NucleotideRegistry(const std::filesystem::path& config_file) {
    std::ifstream in(config_file);
    if (!in.is_open()) {
        throw std::runtime_error("NucleotideRegistry: cannot open " +
                                 config_file.string());
    }
    const nlohmann::json j = nlohmann::json::parse(in);
    const auto& cats = j.at("modified_nucleotides");
    for (const auto& [category, entries] : cats.items()) {
        for (const auto& [name, info] : entries.items()) {
            NucleotideInfo ni;
            const std::string code = info.value("code", std::string{});
            ni.one_letter = code.empty() ? '?' : code[0];
            ni.canonical = canonical_from_type(info.value("type", std::string{}));
            ni.is_purine = info.value("is_purine", false);
            ni.is_modified =
                std::islower(static_cast<unsigned char>(ni.one_letter)) != 0;
            table_[name] = ni;  // later categories override (matches load order)
        }
    }
}

std::optional<NucleotideInfo>
NucleotideRegistry::lookup(std::string_view res_name) const {
    const auto it = table_.find(std::string(res_name));
    if (it == table_.end()) return std::nullopt;
    return it->second;
}

std::string NucleotideRegistry::template_filename(std::string_view res_name) const {
    const auto info = lookup(res_name);
    if (!info || info->canonical == '?') return "";
    const char upper = info->canonical;
    if (info->is_modified) {
        std::string f = "Atomic.";
        f += static_cast<char>(std::tolower(static_cast<unsigned char>(upper)));
        f += ".pdb";
        return f;
    }
    return std::string("Atomic_") + upper + ".pdb";
}

}  // namespace pairfinder::core
