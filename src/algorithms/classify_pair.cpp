/**
 * @file classify_pair.cpp
 * @brief Full LW classification orchestration.
 */
#include <pairfinder/algorithms/classify_pair.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_set>

#include <pairfinder/algorithms/edge_classifier.hpp>
#include <pairfinder/algorithms/saenger.hpp>

namespace pairfinder::algorithms::classification {

namespace {

using core::BaseTyping;
using core::Residue;
using hbond::HBond;

constexpr double kRmsdResolutionThreshold = 0.15;

int face_priority(char f) { return f == 'W' ? 0 : (f == 'H' ? 1 : 2); }

std::string canonicalize_lw(char orientation, char f1, char f2) {
    std::string o(1, orientation);
    if (face_priority(f1) <= face_priority(f2)) return o + f1 + f2;
    return o + f2 + f1;
}

std::optional<std::string> resolve_by_o2prime(const std::array<std::string, 2>& classes,
                                              const Residue& res1, const Residue& res2,
                                              const std::vector<HBond>& hbonds) {
    if (hbonds.empty()) return std::nullopt;
    bool r1 = false, r2 = false;
    for (const auto& hb : hbonds) {
        if (hb.donor_atom == "O2'" && hb.donor_res_id == res1.res_id()) r1 = true;
        else if (hb.acceptor_atom == "O2'" && hb.acceptor_res_id == res1.res_id()) r1 = true;
        if (hb.donor_atom == "O2'" && hb.donor_res_id == res2.res_id()) r2 = true;
        else if (hb.acceptor_atom == "O2'" && hb.acceptor_res_id == res2.res_id()) r2 = true;
    }
    if (!r1 && !r2) return std::nullopt;
    const char e00 = classes[0][1], e01 = classes[0][2];
    const char e10 = classes[1][1], e11 = classes[1][2];
    const std::array<std::pair<char, char>, 2> edges = {std::pair{e00, e10}, std::pair{e01, e11}};
    for (int pos = 0; pos < 2; ++pos) {
        if (edges[pos].first == edges[pos].second) continue;  // not a differing position
        const bool s0 = edges[pos].first == 'S';
        const bool s1 = edges[pos].second == 'S';
        if (s0 == s1) continue;  // not exactly one S
        if (r1 || r2) return s0 ? classes[0] : classes[1];
    }
    return std::nullopt;
}

std::optional<std::string> resolve_ambiguity(const std::array<std::string, 2>& classes,
                                             const std::array<std::optional<double>, 2>& rmsds,
                                             const Residue& res1, const Residue& res2,
                                             const std::vector<HBond>& hbonds,
                                             const BaseTyping& typing) {
    if (!hbonds.empty()) {
        const SaengerMatch sm = match_saenger(res1, res2, hbonds, typing);
        if (sm.found) {
            const char orientation = classes[0][0];
            const std::string sc = canonicalize_lw(orientation, sm.face1, sm.face2);
            if (sc == classes[0] || sc == classes[1]) return sc;
        }
    }
    if (auto o = resolve_by_o2prime(classes, res1, res2, hbonds)) return o;
    if (rmsds[0] && rmsds[1]) {
        const double diff = std::abs(*rmsds[0] - *rmsds[1]);
        if (diff > kRmsdResolutionThreshold) return *rmsds[0] < *rmsds[1] ? classes[0] : classes[1];
    }
    return std::nullopt;
}

// Sugar-edge / HH geometry fallback scores.
double round3(double v) { return std::round(v * 1000.0) / 1000.0; }
double sugar_fallback(double rmsd) {
    return rmsd >= 3.0 ? 0.0 : round3(0.30 * std::max(0.0, 1.0 - rmsd / 3.0));
}
double hh_fallback(double rmsd) {
    return rmsd >= 4.0 ? 0.0 : round3(0.35 * std::max(0.0, 1.0 - rmsd / 4.0));
}

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};

std::vector<ScoringHBond> to_scoring_hbonds(const std::vector<hbond::HBond>& strict) {
    std::vector<ScoringHBond> out;
    out.reserve(strict.size());
    for (const auto& hb : strict) {
        const bool bb = kBaseAtoms.count(hb.donor_atom) && kBaseAtoms.count(hb.acceptor_atom);
        out.push_back({hb.donor_atom, hb.acceptor_atom, hb.distance, hb.alignment_score, bb});
    }
    return out;
}

// score==0 geometry fallback for a class (sugar / HH).
double apply_fallback(double score, const std::string& cls, std::optional<double> rmsd) {
    if (score != 0.0 || !rmsd) return score;
    const bool has_sugar = cls.size() == 3 && cls.find('S') != std::string::npos;
    const bool has_hh = cls.size() == 3 && cls[1] == 'H' && cls[2] == 'H';
    if (has_sugar) return sugar_fallback(*rmsd);
    if (has_hh) return hh_fallback(*rmsd);
    return score;
}

}  // namespace

double compute_selection_quality(const ValidationResult& validation, const std::string& sequence,
                                 const std::vector<hbond::HBond>& strict,
                                 const ClassifiedPair& cp, const QualityScorer& qscorer,
                                 const HBondPatterns& patterns, const core::BaseTyping& typing) {
    const std::vector<ScoringHBond> hbonds = to_scoring_hbonds(strict);
    if (cp.ambiguous && cp.ambig_classes.size() == 2) {
        double best = 0.0;
        for (std::size_t i = 0; i < 2; ++i) {
            double s = qscorer.compute_score(validation, sequence, hbonds, cp.ambig_rmsds[i],
                                             cp.ambig_classes[i], patterns, typing);
            s = apply_fallback(s, cp.ambig_classes[i], cp.ambig_rmsds[i]);
            best = std::max(best, s);
        }
        return best;
    }
    double s = qscorer.compute_score(validation, sequence, hbonds, cp.template_rmsd, cp.lw_class,
                                     patterns, typing);
    if (!cp.lw_class.empty()) s = apply_fallback(s, cp.lw_class, cp.template_rmsd);
    return s;
}

ClassifiedPair classify_pair(const core::ReferenceFrame& frame1,
                             const core::ReferenceFrame& frame2, const core::Residue& res1_in,
                             const core::Residue& res2_in, const std::vector<HBond>& voters,
                             const std::vector<HBond>& strict, const TemplateAligner& aligner,
                             const BaseTyping& typing) {
    LwResult geo = classify_lw_class(frame1, frame2, res1_in, res2_in, voters, typing);
    ClassifiedPair out;
    out.swapped = geo.swapped;

    // Apply the swap to residue order for template/ambiguity (frames unused there).
    const Residue* r1 = &res1_in;
    const Residue* r2 = &res2_in;
    if (geo.swapped) std::swap(r1, r2);

    const std::string& lw = geo.lw_class;
    const auto pipe = lw.find('|');
    if (pipe != std::string::npos) {
        std::array<std::string, 2> classes = {lw.substr(0, pipe), lw.substr(pipe + 1)};
        std::array<std::optional<double>, 2> rmsds = {
            aligner.compute_template_rmsd(*r1, *r2, classes[0], typing),
            aligner.compute_template_rmsd(*r1, *r2, classes[1], typing)};
        auto resolved = resolve_ambiguity(classes, rmsds, *r1, *r2, strict, typing);
        // A pure W/H boundary ambiguity (no S edge) is genuinely undetermined,
        // so keep it rather than committing via template/Saenger. The W/S path
        // is resolved.
        if (classes[0].find('S') == std::string::npos &&
            classes[1].find('S') == std::string::npos)
            resolved = std::nullopt;
        if (resolved) {
            out.lw_class = *resolved;
            out.template_rmsd = rmsds[*resolved == classes[0] ? 0 : 1];
        } else {
            const bool cww0 = classes[0] == "cWW", cww1 = classes[1] == "cWW";
            const bool any_none = !rmsds[0] || !rmsds[1];
            if ((cww0 || cww1) && any_none) {
                const int cww_idx = cww0 ? 0 : 1;
                const int alt_idx = 1 - cww_idx;
                // If cWW's OWN template is missing but the alternative has a
                // good fit, commit the verifiable class instead of a cWW we
                // can't support.
                constexpr double kAltCommitRmsd = 1.5;
                if (!rmsds[cww_idx] && rmsds[alt_idx] &&
                    *rmsds[alt_idx] < kAltCommitRmsd) {
                    out.lw_class = classes[alt_idx];
                    out.template_rmsd = rmsds[alt_idx];
                } else {
                    out.lw_class = "cWW";
                    out.template_rmsd = rmsds[cww_idx];
                }
            } else {
                out.lw_class = lw;  // keep pipe form
                out.template_rmsd = rmsds[0];
                out.ambiguous = true;
                out.ambig_classes = {classes[0], classes[1]};
                out.ambig_rmsds = {rmsds[0], rmsds[1]};
            }
        }
    } else {
        out.lw_class = lw;
        out.template_rmsd = aligner.compute_template_rmsd(*r1, *r2, lw, typing);
    }
    return out;
}

}  // namespace pairfinder::algorithms::classification
