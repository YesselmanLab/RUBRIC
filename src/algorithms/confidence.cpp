/**
 * @file confidence.cpp
 * @brief Intrinsic LW confidence.
 */
#include <pairfinder/algorithms/confidence.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <unordered_map>
#include <unordered_set>

namespace pairfinder::algorithms::classification {

namespace {

using hbond::HBond;

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};

constexpr double kNoTemplateScore = 0.05;
constexpr double kPatternMatchFloor = 0.05;
constexpr double kTemplateRhoDefault = 2.5;

double template_rho(const std::string& lw) {
    static const std::unordered_map<std::string, double> rho = {
        {"cWW", 1.0},  {"tWW", 2.67}, {"cWH", 2.65}, {"tWH", 3.0},  {"cHW", 2.65},
        {"tHW", 3.0},  {"cWS", 2.93}, {"tWS", 2.38}, {"cSW", 2.93}, {"tSW", 2.38},
        {"cHH", 3.59}, {"tHH", 1.1},  {"cHS", 2.9},  {"tHS", 1.65}, {"cSH", 2.9},
        {"tSH", 1.65}, {"cSS", 2.81}, {"tSS", 2.71}};
    const auto it = rho.find(lw);
    return it == rho.end() ? kTemplateRhoDefault : it->second;
}

double template_fit_score(const std::optional<double>& template_rmsd, const std::string& lw) {
    if (!template_rmsd) return kNoTemplateScore;
    const double z = *template_rmsd / template_rho(lw);
    return std::exp(-z * z);
}

double edge_pattern_match_score(const std::vector<HBond>& hbonds, const std::string& lw,
                                const std::string& sequence, const HBondPatterns& patterns,
                                const core::BaseTyping& typing) {
    const auto expected = patterns.get_expected_hbonds(lw, sequence, typing);
    int n_detected = 0;
    std::vector<const HBond*> base_bonds;
    for (const auto& hb : hbonds) {
        if (kBaseAtoms.count(hb.donor_atom) && kBaseAtoms.count(hb.acceptor_atom)) {
            base_bonds.push_back(&hb);
            ++n_detected;
        }
    }
    if (expected.empty()) {
        if (n_detected == 0) return 0.0;
        return std::min(1.0, n_detected / 2.0);
    }
    if (n_detected == 0) return 0.0;
    int matched = 0;
    for (const HBond* hb : base_bonds)
        if (is_hbond_match({hb->donor_atom, hb->acceptor_atom}, expected)) ++matched;
    if (matched == 0) return kPatternMatchFloor;
    const double precision = std::min(1.0, static_cast<double>(matched) / n_detected);
    const double recall = std::min(1.0, static_cast<double>(matched) / expected.size());
    return 2.0 * precision * recall / (precision + recall);
}

double intrinsic_confidence(const std::vector<HBond>& hbonds,
                            const std::optional<double>& template_rmsd, const std::string& lw,
                            const std::string& sequence, const HBondPatterns& patterns,
                            const core::BaseTyping& typing) {
    const double s_geom = template_fit_score(template_rmsd, lw);
    const double s_hbond = edge_pattern_match_score(hbonds, lw, sequence, patterns, typing);
    const double product = s_geom * s_hbond;
    return product <= 0.0 ? 0.0 : std::sqrt(product);
}

double round3(double v) { return std::round(v * 1000.0) / 1000.0; }

const std::array<const char*, 12> kAllLwClasses = {"cWW", "cWH", "cWS", "cHH", "cHS", "cSS",
                                                   "tWW", "tWH", "tWS", "tHH", "tHS", "tSS"};

}  // namespace

ConfidenceResult compute_classification_confidence(
    const core::Residue& res1, const core::Residue& res2, const std::string& sequence,
    const std::string& primary_class, const std::vector<HBond>& hbonds,
    const TemplateAligner& aligner, const HBondPatterns& patterns,
    const core::BaseTyping& typing) {
    ConfidenceResult out;
    const char orientation = primary_class.empty() ? 'c' : primary_class[0];

    std::unordered_map<std::string, double> rmsd_by_class;
    int class_rmsd_count = 0;
    for (const char* cls : kAllLwClasses) {
        if (cls[0] != orientation) continue;
        const auto r = aligner.compute_template_rmsd(res1, res2, cls, typing);
        if (r) {
            rmsd_by_class[cls] = *r;
            ++class_rmsd_count;
        }
    }

    const bool geo_uncertain = primary_class.find('|') != std::string::npos;
    out.is_ambiguous = geo_uncertain && class_rmsd_count >= 2;

    std::vector<std::string> named;
    if (out.is_ambiguous) {
        const auto pipe = primary_class.find('|');
        std::vector<std::string> parts = {primary_class.substr(0, pipe),
                                          primary_class.substr(pipe + 1)};
        std::stable_sort(parts.begin(), parts.end(), [&](const std::string& a, const std::string& b) {
            const double ra = rmsd_by_class.count(a) ? rmsd_by_class.at(a)
                                                     : std::numeric_limits<double>::infinity();
            const double rb = rmsd_by_class.count(b) ? rmsd_by_class.at(b)
                                                     : std::numeric_limits<double>::infinity();
            return ra < rb;
        });
        named = {parts[0], parts[1]};
        out.display_label = "ambiguous (" + parts[0] + "|" + parts[1] + ")";
    } else {
        out.display_label = primary_class;
        named = {primary_class};
    }

    for (const auto& cls : named) {
        std::optional<double> tmpl;
        if (auto it = rmsd_by_class.find(cls); it != rmsd_by_class.end()) tmpl = it->second;
        out.confidences.push_back(
            round3(intrinsic_confidence(hbonds, tmpl, cls, sequence, patterns, typing)));
    }
    return out;
}

}  // namespace pairfinder::algorithms::classification
