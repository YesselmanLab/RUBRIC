/**
 * @file template_aligner.cpp
 * @brief Template-RMSD lookup + alignment.
 */
#include <pairfinder/algorithms/template_aligner.hpp>

#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <vector>

#include <pairfinder/geometry/superposition.hpp>

namespace pairfinder::algorithms::classification {

namespace {

using geometry::Vector3d;

const std::array<const char*, 9> kRingAtoms = {"C2", "C4", "C5", "C6", "N1",
                                               "N3", "N7", "C8", "N9"};

int face_priority(char f) { return f == 'W' ? 0 : (f == 'H' ? 1 : 2); }

bool is_symmetric(const std::string& lw) {
    return lw == "cWW" || lw == "tWW" || lw == "cHH" || lw == "tHH" || lw == "cSS" ||
           lw == "tSS";
}

std::string trim(const std::string& s) {
    std::size_t a = s.find_first_not_of(' ');
    if (a == std::string::npos) return "";
    std::size_t b = s.find_last_not_of(' ');
    return s.substr(a, b - a + 1);
}

}  // namespace

std::optional<std::filesystem::path>
TemplateAligner::find_template(const std::string& lw_class, const std::string& seq) const {
    if (seq.size() < 2) return std::nullopt;
    const std::string memo_key = lw_class + '\x01' + seq;
    if (const auto it = found_.find(memo_key); it != found_.end()) return it->second;
    return found_[memo_key] = locate_template(lw_class, seq);
}

std::optional<std::filesystem::path>
TemplateAligner::locate_template(const std::string& lw_class, const std::string& seq) const {
    // Idealized: <dir>/<lw>/<seq>.pdb or <dir>/<lw>/<s0>_<lower s1>.pdb
    const auto base = idealized_dir_ / lw_class;
    const std::string p1 = seq + ".pdb";
    std::string p2;
    p2 += seq[0];
    p2 += '_';
    p2 += static_cast<char>(std::tolower(static_cast<unsigned char>(seq[1])));
    p2 += ".pdb";
    for (const auto& p : {p1, p2}) {
        const auto path = base / p;
        if (std::filesystem::exists(path)) return path;
    }
    // Exemplar (flat) fallbacks.
    const char s0 = seq[0], s1 = seq[1];
    const char s0l = static_cast<char>(std::tolower(static_cast<unsigned char>(s0)));
    std::vector<std::string> pats = {
        std::string() + s0 + '-' + s1 + '-' + lw_class + ".pdb",
        std::string() + s0 + "plus" + s1 + '-' + lw_class + ".pdb",
        std::string() + s0l + '-' + s1 + '-' + lw_class + ".pdb",
        seq + "-" + lw_class + ".pdb"};
    for (const auto& p : pats) {
        const auto path = exemplar_dir_ / p;
        if (std::filesystem::exists(path)) return path;
    }
    return std::nullopt;
}

const TemplateAligner::Template*
TemplateAligner::load_template(const std::filesystem::path& path) const {
    const std::string key = path.string();
    if (auto it = cache_.find(key); it != cache_.end()) return &it->second;
    Template tmpl;
    std::ifstream in(path);
    if (!in.is_open()) return nullptr;
    std::string line;
    while (std::getline(in, line)) {
        if (line.rfind("ATOM", 0) != 0) continue;
        if (line.size() < 54) continue;
        const std::string atom_name = trim(line.substr(12, 4));
        int res_seq = 0;
        try {
            res_seq = std::stoi(trim(line.substr(22, 4)));
        } catch (const std::exception&) {
            continue;
        }
        Vector3d c{std::stod(line.substr(30, 8)), std::stod(line.substr(38, 8)),
                   std::stod(line.substr(46, 8))};
        (res_seq == 1 ? tmpl.first : tmpl.second)[atom_name] = c;
    }
    auto [it, _] = cache_.emplace(key, std::move(tmpl));
    return &it->second;
}

double TemplateAligner::align_to_template(const core::Residue& res1, const core::Residue& res2,
                                          const std::filesystem::path& path) const {
    const Template* tmpl = load_template(path);
    if (tmpl == nullptr) return std::numeric_limits<double>::infinity();
    std::vector<Vector3d> tpts, qpts;
    auto collect = [&](const AtomMap& tatoms, const core::Residue& target) {
        for (const char* a : kRingAtoms) {
            auto it = tatoms.find(a);
            const core::Atom* qa = target.get_atom(a);
            if (it != tatoms.end() && qa != nullptr) {
                tpts.push_back(it->second);
                qpts.push_back(qa->coords);
            }
        }
    };
    collect(tmpl->first, res1);
    collect(tmpl->second, res2);
    if (tpts.size() < 4) return std::numeric_limits<double>::infinity();
    const auto fit = geometry::superpose(tpts, qpts);  // optimal RMSD == Kabsch
    return fit.valid ? fit.rms : std::numeric_limits<double>::infinity();
}

std::optional<double> TemplateAligner::compute_template_rmsd(const core::Residue& res1_in,
                                                            const core::Residue& res2_in,
                                                            const std::string& lw_in,
                                                            const core::BaseTyping& typing) const {
    const core::Residue* res1 = &res1_in;
    const core::Residue* res2 = &res2_in;
    std::string lw = lw_in;
    if (lw.size() == 3 && face_priority(lw[1]) > face_priority(lw[2])) {
        lw = std::string() + lw[0] + lw[2] + lw[1];
        std::swap(res1, res2);
    }
    const std::string seq = typing.normalize(res1->base_type()) + typing.normalize(res2->base_type());

    if (auto path = find_template(lw, seq)) return align_to_template(*res1, *res2, *path);
    if (is_symmetric(lw)) {
        const std::string rev = std::string() + seq[1] + seq[0];
        if (auto path = find_template(lw, rev)) return align_to_template(*res2, *res1, *path);
    }
    return std::nullopt;
}

}  // namespace pairfinder::algorithms::classification
