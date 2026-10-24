/**
 * @file base_typing.cpp
 * @brief Base-type normalization + glycosidic-N selection.
 */
#include <pairfinder/core/base_typing.hpp>

#include <array>
#include <cctype>
#include <fstream>
#include <set>

#include <nlohmann/json.hpp>

namespace pairfinder::core {

namespace {

std::string upper_strip(std::string_view s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == ' ' || c == '\t') continue;
        out.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(c))));
    }
    return out;
}

bool in(const std::string& s, std::initializer_list<const char*> set) {
    for (const char* x : set)
        if (s == x) return true;
    return false;
}

// Curated modified-base map - takes precedence over the DB map.
const std::unordered_map<std::string, std::string>& curated_map() {
    static const std::unordered_map<std::string, std::string> m = {
        {"2MG", "G"}, {"7MG", "G"}, {"M2G", "G"}, {"OMG", "G"}, {"YG", "G"},
        {"PSU", "U"}, {"H2U", "U"}, {"5MU", "U"}, {"4SU", "U"},
        {"5MC", "C"}, {"OMC", "C"}, {"1MA", "A"}, {"MIA", "A"}, {"6MA", "A"},
        {"INO", "I"}};
    return m;
}

// DNA prefix map: DT -> T here (NOT U).
const std::unordered_map<std::string, std::string>& dna_prefix_map() {
    static const std::unordered_map<std::string, std::string> m = {
        {"DA", "A"}, {"DG", "G"}, {"DC", "C"}, {"DT", "T"}};
    return m;
}

// Glycosidic-N atom per base.
const std::unordered_map<std::string, std::string>& glycosidic_map() {
    static const std::unordered_map<std::string, std::string> m = {
        {"A", "N9"}, {"G", "N9"}, {"I", "N9"}, {"C", "N1"}, {"U", "N1"},
        {"T", "N1"}, {"P", "C5"}, {"DA", "N9"}, {"DG", "N9"}, {"DC", "N1"},
        {"DT", "N1"}};
    return m;
}

// Parent base inferred from a ring-discriminating donor/acceptor atom signature.
std::string infer_parent(const nlohmann::json& entry) {
    std::set<std::string> atoms;
    for (const char* role : {"donors", "acceptors"}) {
        if (!entry.contains(role)) continue;
        for (const auto& a : entry[role]) atoms.insert(a.value("atom", std::string{}));
    }
    if (atoms.count("N7") || atoms.count("N9"))
        return (atoms.count("O6") || atoms.count("N2")) ? "G" : "A";
    if (atoms.count("N4")) return "C";
    if (atoms.count("O4")) return "U";
    return "";
}

}  // namespace

BaseTyping::BaseTyping(const std::filesystem::path& ligand_db_path, bool use_cache) {
    // Fast path: load the precomputed name->parent map (an exact serialization of
    // the table built below) instead of parsing the 1.5 MB ligand DB every startup.
    if (use_cache && load_cache(ligand_db_path.parent_path() / "base_typing_cache.json"))
        return;

    std::ifstream in_file(ligand_db_path);
    if (!in_file.is_open()) return;
    nlohmann::json db;
    try {
        db = nlohmann::json::parse(in_file);
    } catch (const std::exception&) {
        return;
    }
    static const std::unordered_map<std::string, std::string> dna_to_rna = {
        {"DA", "A"}, {"DG", "G"}, {"DC", "C"}, {"DT", "U"}, {"DU", "U"}};
    for (const auto& [code, entry] : db.items()) {
        if (code == "_metadata" || !entry.is_object()) continue;
        const std::string type = entry.value("type", std::string{});
        // Free nucleotides and cofactors (ATP, GTP, ...) are typed LIGAND. They
        // are recorded so the backbone scorer can refuse to build a suite on
        // them; they are deliberately NOT added to the parent map below.
        if (type == "LIGAND") {
            ligand_codes_.insert(code);
            continue;
        }
        if (type != "RNA" && type != "DNA") continue;
        // Record DNA-ness before the parent is collapsed to its RNA equivalent
        // below: normalize() wants DT -> U for base identity, but callers that
        // ask "is this DNA" must still be able to tell them apart.
        if (type == "DNA") dna_codes_.insert(code);

        if (entry.contains("parent") && entry["parent"].is_string()) {
            const std::string parent = entry["parent"].get<std::string>();
            if (in(parent, {"A", "G", "C", "U"})) {
                extended_[code] = parent;
            } else if (auto it = dna_to_rna.find(parent); it != dna_to_rna.end()) {
                extended_[code] = it->second;
            }
        } else {  // parent null/absent -> infer from atoms
            const std::string inferred = infer_parent(entry);
            if (!inferred.empty()) extended_[code] = inferred;
        }
    }
}

void BaseTyping::save_cache(const std::filesystem::path& path) const {
    nlohmann::json parents = nlohmann::json::object();
    for (const auto& [code, parent] : extended_) parents[code] = parent;
    nlohmann::json dna = nlohmann::json::array();
    for (const auto& code : std::set<std::string>(dna_codes_.begin(), dna_codes_.end()))
        dna.push_back(code);
    nlohmann::json lig = nlohmann::json::array();
    for (const auto& code : std::set<std::string>(ligand_codes_.begin(), ligand_codes_.end()))
        lig.push_back(code);
    nlohmann::json root = {{"parents", std::move(parents)},
                           {"dna", std::move(dna)},
                           {"ligand", std::move(lig)}};
    std::ofstream out(path);
    out << root.dump(2) << "\n";
}

bool BaseTyping::load_cache(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in.is_open()) return false;
    nlohmann::json root;
    try {
        in >> root;
    } catch (const nlohmann::json::exception&) {
        return false;
    }
    if (!root.is_object()) return false;
    // Cache layout is {"parents": {code: parent}, "dna": [code, ...]}. A cache
    // without the DNA set is rejected so the caller re-parses the ligand DB
    // rather than silently losing modified-DNA typing.
    if (!root.contains("parents") || !root.contains("dna") || !root.contains("ligand"))
        return false;
    if (!root["parents"].is_object() || !root["dna"].is_array() ||
        !root["ligand"].is_array())
        return false;

    std::unordered_map<std::string, std::string> table;
    for (auto it = root["parents"].begin(); it != root["parents"].end(); ++it) {
        if (!it->is_string()) return false;
        table[it.key()] = it->get<std::string>();
    }
    std::unordered_set<std::string> dna;
    for (const auto& v : root["dna"]) {
        if (!v.is_string()) return false;
        dna.insert(v.get<std::string>());
    }
    std::unordered_set<std::string> lig;
    for (const auto& v : root["ligand"]) {
        if (!v.is_string()) return false;
        lig.insert(v.get<std::string>());
    }
    extended_ = std::move(table);
    dna_codes_ = std::move(dna);
    ligand_codes_ = std::move(lig);
    return true;
}

bool BaseTyping::is_known_nucleotide(std::string_view residue_name) const {
    const std::string name = upper_strip(residue_name);
    if (in(name, {"A", "G", "C", "U", "T", "I", "P"})) return true;
    if (dna_prefix_map().count(name) || curated_map().count(name)) return true;
    return extended_.count(name) > 0;
}

bool BaseTyping::is_ligand(std::string_view residue_name) const {
    return ligand_codes_.count(upper_strip(residue_name)) > 0;
}

bool BaseTyping::is_dna(std::string_view residue_name) const {
    const std::string name = upper_strip(residue_name);
    // Floor: the literal DNA names, so DNA is still recognized even when the
    // ligand DB is unavailable.
    if (name.size() == 2 && name[0] == 'D' &&
        in(std::string(1, name[1]), {"A", "C", "G", "T", "U"}))
        return true;
    return dna_codes_.count(name) > 0;
}

std::string BaseTyping::normalize(std::string_view residue_name) const {
    const std::string name = upper_strip(residue_name);
    if (in(name, {"A", "G", "C", "U", "T", "I", "P"})) return name;
    if (auto it = dna_prefix_map().find(name); it != dna_prefix_map().end())
        return it->second;
    if (auto it = curated_map().find(name); it != curated_map().end())
        return it->second;
    if (auto it = extended_.find(name); it != extended_.end()) return it->second;
    return name;
}

bool BaseTyping::is_purine(std::string_view base_type) const {
    return in(normalize(base_type), {"A", "G", "I"});
}

bool BaseTyping::is_pyrimidine(std::string_view base_type) const {
    return in(normalize(base_type), {"C", "U", "T"});
}

std::string BaseTyping::glycosidic_n_name(std::string_view base_type) const {
    const std::string key = upper_strip(base_type);
    if (auto it = glycosidic_map().find(key); it != glycosidic_map().end())
        return it->second;
    const std::string parent = normalize(key);
    if (auto it = glycosidic_map().find(parent); it != glycosidic_map().end())
        return it->second;
    return "";
}

}  // namespace pairfinder::core
