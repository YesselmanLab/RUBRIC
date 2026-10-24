/**
 * @file hbond_patterns.cpp
 * @brief Expected LW H-bond patterns + matching.
 */
#include <pairfinder/algorithms/hbond_patterns.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <regex>
#include <set>
#include <sstream>

#include <nlohmann/json.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms::classification {

namespace {

using geometry::Vector3d;

constexpr double kHbondMaxDistance = 3.5;

bool is_standard_base(const std::string& b) {
    return b == "A" || b == "G" || b == "C" || b == "U" || b == "T" || b == "I" || b == "P";
}

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// Normalize a pair sequence to the key used by the curated pattern table.
std::string normalize_sequence(const std::string& sequence, const core::BaseTyping& typing) {
    const std::string seq = upper(sequence);
    if (seq.size() == 2 && is_standard_base(std::string(1, seq[0])) &&
        is_standard_base(std::string(1, seq[1])))
        return seq;
    for (std::size_t split = 1; split < seq.size(); ++split) {
        const std::string left = typing.normalize(seq.substr(0, split));
        const std::string right = typing.normalize(seq.substr(split));
        if (is_standard_base(left) && is_standard_base(right)) return left + right;
    }
    return seq;
}

// Curated H-bond patterns. lw -> seq -> list of (donor, acceptor).
const std::unordered_map<std::string, std::unordered_map<std::string, std::vector<HBondPattern>>>&
curated_patterns() {
    static const std::unordered_map<std::string,
                                    std::unordered_map<std::string, std::vector<HBondPattern>>>
        m = {
            {"cWW",
             {{"GC", {{"N1", "N3"}, {"N2", "O2"}, {"O6", "N4"}}},
              {"CG", {{"N4", "O6"}, {"N3", "N1"}, {"O2", "N2"}}},
              {"AU", {{"N6", "O4"}, {"N1", "N3"}}},
              {"UA", {{"O4", "N6"}, {"N3", "N1"}}},
              {"AT", {{"N6", "O4"}, {"N1", "N3"}}},
              {"TA", {{"O4", "N6"}, {"N3", "N1"}}},
              {"GU", {{"N1", "O2"}}},
              {"UG", {{"O2", "N1"}}},
              {"GT", {{"N1", "O2"}}},
              {"TG", {{"O2", "N1"}}},
              {"GA", {{"O6", "N6"}, {"N1", "N1"}}},
              {"AG", {{"N6", "O6"}, {"N1", "N1"}}},
              {"AA", {{"N6", "N1"}, {"N1", "N6"}}},
              {"GG", {{"O6", "N1"}, {"N1", "O6"}}},
              {"UC", {{"N3", "N3"}}},
              {"CU", {{"N3", "N3"}}},
              {"UU", {{"N3", "O4"}, {"O4", "N3"}}},
              {"CC", {{"N4", "N3"}, {"N3", "N4"}}}}},
            {"tWW",
             {{"GC", {{"N2", "O2"}, {"N1", "N3"}}},
              {"CG", {{"O2", "N2"}, {"N3", "N1"}}},
              {"AU", {{"N6", "O4"}}},
              {"UA", {{"O4", "N6"}}}}},
            {"cWH", {{"UA", {{"N3", "N7"}, {"N6", "O4"}}}}},
            {"tWH",
             {{"UA", {{"N3", "N7"}, {"N6", "O2"}}},
              {"GU", {{"N1", "O6"}, {"O6", "N7"}}},
              {"UG", {{"N3", "N7"}}},
              {"GG", {{"N1", "N7"}, {"O6", "N2"}}}}},
            {"cWS",
             {{"GC", {{"N2", "O2"}}},
              {"CG", {{"O2", "N2"}}},
              {"AU", {{"N6", "O2"}}},
              {"UA", {{"O2", "N6"}}}}},
            {"tWS",
             {{"GU", {{"N2", "O2"}}},
              {"UG", {{"O2", "N2"}}},
              {"GA", {{"N2", "N3"}}},
              {"AG", {{"N3", "N2"}}}}},
            {"cHH", {{"AA", {{"N6", "N7"}}}, {"GG", {{"N2", "N7"}, {"O6", "N7"}}}}},
            {"tHH", {{"AA", {{"N6", "N7"}}}, {"GG", {{"N2", "N7"}}}}},
            {"cHS", {{"GA", {{"N2", "N3"}}}, {"AG", {{"N7", "N6"}}}, {"GG", {{"N2", "N3"}}}}},
            {"tHS", {{"AA", {{"N6", "N1"}}}, {"GG", {{"N2", "N7"}}}}},
            {"cSS", {{"GC", {{"N2", "O2"}}}, {"AU", {{"O2", "O2"}}}}},
            {"tSS", {{"GA", {{"N3", "N3"}}}, {"AG", {{"N3", "N3"}}}}}};
    return m;
}

std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(" \t");
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(" \t");
    return s.substr(a, b - a + 1);
}

// Patterns from a template's REMARK 7 records.
std::vector<HBondPattern> parse_remark7(const std::string& text) {
    static const std::regex re(R"(^REMARK\s+7\s+([A-Z0-9']+)-([A-Z0-9']+)\s*:\s*(?:d=)?(\d+\.\d+))");
    std::set<HBondPattern> seen;
    std::vector<HBondPattern> out;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        std::smatch m;
        if (!std::regex_search(line, m, re)) continue;
        if (std::stod(m[3].str()) > kHbondMaxDistance) continue;
        HBondPattern key{m[1].str(), m[2].str()};
        if (seen.insert(key).second) out.push_back(key);
    }
    return out;
}

// Patterns derived from template geometry (fallback when no REMARK 7).
std::vector<HBondPattern> derive_from_geometry(const std::string& text,
                                               const hbond::HBondChemistry& chem) {
    std::unordered_map<std::string, Vector3d> aa, ab;
    std::string base_a, base_b;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("ATOM", 0) != 0 || line.size() < 54) continue;
        const std::string name = trim(line.substr(12, 4));
        const std::string res = trim(line.substr(17, 3));
        const std::string chain = trim(line.substr(21, 1));
        Vector3d xyz{std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                     std::stod(line.substr(46, 8))};
        if (chain == "A") { aa[name] = xyz; base_a = res; }
        else if (chain == "B") { ab[name] = xyz; base_b = res; }
    }
    if (aa.empty() || ab.empty()) return {};
    std::set<HBondPattern> seen;
    std::vector<HBondPattern> out;
    auto scan = [&](const std::unordered_map<std::string, Vector3d>& don, const std::string& db,
                    const std::unordered_map<std::string, Vector3d>& acc, const std::string& ab2) {
        for (const auto& [dn, dx] : don) {
            if (!chem.donor_capacity(db, dn)) continue;
            for (const auto& [an, ax] : acc) {
                if (!chem.acceptor_capacity(ab2, an)) continue;
                if ((dx - ax).norm() <= kHbondMaxDistance) {
                    HBondPattern key{dn, an};
                    if (seen.insert(key).second) out.push_back(key);
                }
            }
        }
    };
    scan(aa, base_a, ab, base_b);
    scan(ab, base_b, aa, base_a);
    return out;
}

}  // namespace

HBondPatterns::HBondPatterns(std::filesystem::path idealized_dir,
                             const hbond::HBondChemistry& chem, bool use_cache) {
    // Fast path: load the precomputed table (an exact serialization of the scan
    // below) instead of parsing ~170 .pdb files with std::regex every startup.
    if (use_cache && load_cache(idealized_dir / "hbond_patterns_cache.json")) return;

    if (!std::filesystem::exists(idealized_dir)) return;
    for (const auto& class_entry : std::filesystem::directory_iterator(idealized_dir)) {
        if (!class_entry.is_directory()) continue;
        const std::string lw = class_entry.path().filename().string();
        std::unordered_map<std::string, std::vector<HBondPattern>> seq_table;
        for (const auto& f : std::filesystem::directory_iterator(class_entry.path())) {
            if (f.path().extension() != ".pdb") continue;
            std::string seq = f.path().stem().string();
            seq.erase(std::remove(seq.begin(), seq.end(), '_'), seq.end());
            seq = upper(seq);
            std::ifstream in(f.path());
            if (!in.is_open()) continue;
            std::stringstream ss;
            ss << in.rdbuf();
            const std::string text = ss.str();
            std::vector<HBondPattern> pats = parse_remark7(text);
            if (pats.empty()) pats = derive_from_geometry(text, chem);
            if (!pats.empty()) seq_table[seq] = std::move(pats);
        }
        if (!seq_table.empty()) template_table_[lw] = std::move(seq_table);
    }
}

void HBondPatterns::save_cache(const std::filesystem::path& path) const {
    nlohmann::json root = nlohmann::json::object();
    for (const auto& [lw, seq_table] : template_table_) {
        nlohmann::json classes = nlohmann::json::object();
        for (const auto& [seq, pats] : seq_table) {
            nlohmann::json arr = nlohmann::json::array();
            for (const auto& [donor, acceptor] : pats)
                arr.push_back(nlohmann::json::array({donor, acceptor}));
            classes[seq] = std::move(arr);
        }
        root[lw] = std::move(classes);
    }
    std::ofstream out(path);
    out << root.dump(2) << "\n";
}

bool HBondPatterns::load_cache(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in.is_open()) return false;
    nlohmann::json root;
    try {
        in >> root;
    } catch (const nlohmann::json::exception&) {
        return false;
    }
    if (!root.is_object()) return false;

    decltype(template_table_) table;
    for (auto lw_it = root.begin(); lw_it != root.end(); ++lw_it) {
        if (!lw_it->is_object()) return false;
        for (auto seq_it = lw_it->begin(); seq_it != lw_it->end(); ++seq_it) {
            if (!seq_it->is_array()) return false;
            std::vector<HBondPattern> pats;
            for (const auto& pair : *seq_it) {
                if (!pair.is_array() || pair.size() != 2) return false;
                pats.emplace_back(pair[0].get<std::string>(), pair[1].get<std::string>());
            }
            table[lw_it.key()][seq_it.key()] = std::move(pats);
        }
    }
    template_table_ = std::move(table);
    return true;
}

std::vector<HBondPattern> HBondPatterns::get_expected_hbonds(const std::string& lw_class,
                                                            const std::string& sequence,
                                                            const core::BaseTyping& typing) const {
    const std::string norm = normalize_sequence(sequence, typing);
    if (auto cit = template_table_.find(lw_class); cit != template_table_.end()) {
        if (auto sit = cit->second.find(norm); sit != cit->second.end()) return sit->second;
    }
    const auto& cur = curated_patterns();
    if (auto cit = cur.find(lw_class); cit != cur.end()) {
        if (auto sit = cit->second.find(norm); sit != cit->second.end()) return sit->second;
    }
    return {};
}

bool is_hbond_match(const HBondPattern& found, const std::vector<HBondPattern>& expected) {
    for (const auto& [ed, ea] : expected) {
        if (found.first == ed && found.second == ea) return true;
        if (found.first == ea && found.second == ed) return true;
    }
    return false;
}

}  // namespace pairfinder::algorithms::classification
