/**
 * @file score_candidates.cpp
 * @brief Full per-candidate scoring pass + carbonyl-O2 correction.
 */
#include <pairfinder/algorithms/score_candidates.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <pairfinder/algorithms/carbonyl_margin.hpp>
#include <pairfinder/algorithms/classify_pair.hpp>
#include <pairfinder/algorithms/confidence.hpp>

namespace pairfinder::algorithms::classification {

namespace {

using core::Residue;
using hbond::HBond;

const std::unordered_set<std::string> kBaseHbondAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::unordered_set<std::string> kSugarAtoms = {
    "O5'", "C3'", "O3'", "C5'", "C1'", "O2'", "C2'", "C4'", "O4'"};
const std::unordered_set<std::string> kPhosphateAtoms = {
    "O5'", "OP1", "OP2", "O2P", "P", "O1P"};
const std::array<const char*, 6> kPairPriority = {
    "base-base", "base-sugar", "base-phosphate",
    "sugar-sugar", "sugar-phosphate", "phosphate-phosphate"};

std::string atom_group(const std::string& a) {
    if (kBaseHbondAtoms.count(a)) return "base";
    if (kSugarAtoms.count(a)) return "sugar";
    if (kPhosphateAtoms.count(a)) return "phosphate";
    return "unknown";
}

std::string hbond_category(const std::string& donor, const std::string& acceptor) {
    std::string g1 = atom_group(donor), g2 = atom_group(acceptor);
    if (g1 > g2) std::swap(g1, g2);
    return g1 + "-" + g2;
}

std::string pair_category(const std::vector<HBond>& hbonds) {
    if (hbonds.empty()) return "unknown";
    std::unordered_set<std::string> cats;
    for (const auto& hb : hbonds) cats.insert(hbond_category(hb.donor_atom, hb.acceptor_atom));
    for (const char* c : kPairPriority)
        if (cats.count(c)) return c;
    return hbond_category(hbonds[0].donor_atom, hbonds[0].acceptor_atom);
}

// Adjacent residues: same chain, |seq diff| <= 1.
bool are_adjacent(const std::string& a, const std::string& b) {
    auto parse = [](const std::string& r, std::string& chain, long& num) -> bool {
        const auto p1 = r.find('-');
        const auto p2 = r.rfind('-');
        if (p1 == std::string::npos || p2 == p1) return false;
        chain = r.substr(0, p1);
        std::string seq = r.substr(p2 + 1);
        std::string digits;
        for (char c : seq)
            if (std::isdigit(static_cast<unsigned char>(c)) || c == '-') digits += c;
        if (digits.empty()) return false;
        try { num = std::stol(digits); } catch (const std::exception&) { return false; }
        return true;
    };
    std::string c1, c2;
    long n1 = 0, n2 = 0;
    if (!parse(a, c1, n1) || !parse(b, c2, n2)) return false;
    return c1 == c2 && std::abs(n1 - n2) <= 1;
}

bool strong_base_hbond(const std::vector<HBond>& strict, double max_dist) {
    for (const auto& hb : strict)
        if (kBaseHbondAtoms.count(hb.donor_atom) && kBaseHbondAtoms.count(hb.acceptor_atom) &&
            hb.distance <= max_dist)
            return true;
    return false;
}

}  // namespace

std::vector<ScoredCandidate> score_candidates(
    const core::Structure& structure,
    const std::vector<std::pair<std::string, core::ReferenceFrame>>& frames,
    const core::BaseTyping& typing, hbond::HBondChemistry& chem,
    const TemplateAligner& aligner, const HBondPatterns& patterns,
    const QualityScorer& qscorer, double max_hbond_distance) {
    std::unordered_map<std::string, const Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    std::unordered_map<std::string, const core::ReferenceFrame*> frame_by_id;
    for (const auto& [rid, f] : frames) frame_by_id[rid] = &f;

    const auto cands = find_candidates(structure, frames, typing);
    hbond::HBondFinder finder(chem);

    std::vector<ScoredCandidate> out;
    out.reserve(cands.size());
    for (const auto& c : cands) {
        ScoredCandidate sc;
        sc.res_id1 = c.res_id1;
        sc.res_id2 = c.res_id2;
        sc.res_name1 = c.res_name1;
        sc.res_name2 = c.res_name2;
        sc.frame1 = frame_by_id.count(c.res_id1) ? frame_by_id[c.res_id1] : nullptr;
        sc.frame2 = frame_by_id.count(c.res_id2) ? frame_by_id[c.res_id2] : nullptr;
        sc.validation = c.validation;
        if (!c.validation.is_valid) { out.push_back(std::move(sc)); continue; }

        const Residue* res1 = res_by_id.count(c.res_id1) ? res_by_id[c.res_id1] : nullptr;
        const Residue* res2 = res_by_id.count(c.res_id2) ? res_by_id[c.res_id2] : nullptr;
        if (!res1 || !res2) { out.push_back(std::move(sc)); continue; }

        std::vector<HBond> strict = finder.find_between(*res1, *res2);
        sc.hbond_categories.reserve(strict.size());
        for (const auto& hb : strict)
            sc.hbond_categories.push_back(hbond_category(hb.donor_atom, hb.acceptor_atom));
        sc.pair_category = pair_category(strict);
        if (sc.pair_category == "base-base" && are_adjacent(c.res_id1, c.res_id2))
            sc.pair_category = "adjacent-base-base";

        ClassifiedPair classified;
        if (sc.pair_category == "base-base") {
            const auto voters = finder.find_face_voters(*res1, *res2);
            classified = classify_pair(*sc.frame1, *sc.frame2, *res1, *res2, voters, strict,
                                       aligner, typing);
            sc.lw_class = classified.lw_class;
            sc.template_rmsd = classified.template_rmsd;
            if (classified.swapped) {
                std::swap(sc.res_id1, sc.res_id2);
                std::swap(sc.res_name1, sc.res_name2);
                std::swap(sc.frame1, sc.frame2);
                std::swap(res1, res2);
            }
        }
        const std::string sequence = sc.res_name1 + sc.res_name2;
        sc.quality_score = compute_selection_quality(c.validation, sequence, strict, classified,
                                                     qscorer, patterns, typing);
        sc.num_hbonds = static_cast<int>(strict.size());
        sc.has_strong_base_hbond = strong_base_hbond(strict, max_hbond_distance);
        if (sc.pair_category == "base-base" && !sc.lw_class.empty()) {
            const auto conf = compute_classification_confidence(*res1, *res2, sequence,
                                                                sc.lw_class, strict, aligner,
                                                                patterns, typing);
            sc.lw_class_display = conf.display_label;
            sc.lw_class_confidence = conf.confidences;
            sc.is_ambiguous = conf.is_ambiguous;
        }
        out.push_back(std::move(sc));
    }

    // --- Carbonyl-O2 correction (2-pass) ---
    std::unordered_set<std::string> committed;
    for (const auto& sc : out) {
        if (!sc.validation.is_valid || !sc.has_strong_base_hbond) continue;
        const std::string& lw = sc.lw_class;
        if (lw.size() != 3 || lw.find('|') != std::string::npos) continue;
        if (lw[1] == 'W') committed.insert(sc.res_id1);
        if (lw[2] == 'W') committed.insert(sc.res_id2);
    }
    for (auto& sc : out) {
        if (!sc.validation.is_valid || sc.pair_category != "base-base") continue;
        const Residue* res1 = res_by_id.count(sc.res_id1) ? res_by_id[sc.res_id1] : nullptr;
        const Residue* res2 = res_by_id.count(sc.res_id2) ? res_by_id[sc.res_id2] : nullptr;
        if (!res1 || !res2) continue;
        const O2Correction corr =
            correct_cis_o2(sc.lw_class, sc.res_id1, sc.res_id2, *res1, *res2, committed, typing);
        if (!corr.fired) continue;
        sc.precorr_lw = sc.lw_class;
        sc.lw_class = corr.new_class;
        if (corr.swapped) {
            std::swap(sc.res_id1, sc.res_id2);
            std::swap(sc.res_name1, sc.res_name2);
            std::swap(sc.frame1, sc.frame2);
            std::swap(res1, res2);
        }
        const std::string sequence = sc.res_name1 + sc.res_name2;
        sc.template_rmsd = aligner.compute_template_rmsd(*res1, *res2, corr.new_class, typing);
        const std::vector<HBond> strict = finder.find_between(*res1, *res2);
        ClassifiedPair cp;  // corrected classes are single (no ambiguity)
        cp.lw_class = corr.new_class;
        cp.template_rmsd = sc.template_rmsd;
        sc.quality_score =
            compute_selection_quality(sc.validation, sequence, strict, cp, qscorer, patterns, typing);
        const auto conf = compute_classification_confidence(*res1, *res2, sequence, corr.new_class,
                                                            strict, aligner, patterns, typing);
        sc.lw_class_display = conf.display_label;
        sc.lw_class_confidence = conf.confidences;
        sc.is_ambiguous = conf.is_ambiguous;
    }
    return out;
}

}  // namespace pairfinder::algorithms::classification
