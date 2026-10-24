/**
 * @file differential_commands.cpp
 * @brief Implementations of the hidden `dump-*` subcommands (see differential_commands.hpp).
 */
#include "differential/differential_commands.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <pairfinder/algorithms/base_frame_calculator.hpp>
#include <pairfinder/algorithms/candidate_finder.hpp>
#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/algorithms/classify_pair.hpp>
#include <pairfinder/algorithms/confidence.hpp>
#include <pairfinder/algorithms/edge_classifier.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/io/pdb_parser.hpp>
#include <pairfinder/pipeline.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <pairfinder/scoring/scorer.hpp>
#include <pairfinder/scoring/torsions.hpp>

#include "support/cli_support.hpp"

namespace pfcli {
namespace {

// Format the triggered issues of a PairScore as "issue:weight,..." for human-readable output.
std::string format_issues(const pairfinder::scoring::PairScore& ps) {
    if (ps.issues.empty()) return "";
    std::string out;
    for (const auto& ip : ps.issues) {
        if (!out.empty()) out += ',';
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%s:%.4f", ip.issue.c_str(), ip.weight);
        out += buf;
    }
    return out;
}

}  // namespace

// Emit "<res_id>\t<rms>\tR00..R22\ttx\tty\ttz" for diffing frames.
int dump_frames(const std::string& path) {
    const auto structure = pairfinder::io::parse_pdb(path);
    const pairfinder::core::NucleotideRegistry registry(nucleotides_config());
    pairfinder::algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = pairfinder::algorithms::calculate_all_frames(structure, lib);
    for (const auto& [res_id, f] : frames) {
        std::printf("%s\t%.6f", res_id.c_str(), f.rms_fit);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) std::printf("\t%.6f", f.rotation[i][j]);
        std::printf("\t%.6f\t%.6f\t%.6f\n", f.origin.x, f.origin.y, f.origin.z);
    }
    return 0;
}

// Full pipeline + empirical pair scoring (Cerny method). Emits per scorable selected pair:
// "<res_a>\t<res_b>\t<lw_class>\t<score>\t<penalty>\t<issues>" sorted.
// issues = comma-separated "param:weighted_penalty" in descending penalty order.
int dump_scores(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    auto& structure = prep.structure;
    auto& chem = pipeline.chem();
    const auto& typing = pipeline.typing();
    auto scored = pipeline.classify(prep);
    const auto selected = pipeline.select(std::move(scored), structure);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    auto scorer = pipeline.make_scorer();
    algorithms::hbond::HBondFinder finder(chem);

    std::vector<std::string> rows;
    for (const auto& sc : selected) {
        if (!scoring::Scorer::is_scorable(sc, typing)) continue;
        const core::Residue* r1 = res_by_id.count(sc.res_id1) ? res_by_id[sc.res_id1] : nullptr;
        const core::Residue* r2 = res_by_id.count(sc.res_id2) ? res_by_id[sc.res_id2] : nullptr;
        if (!r1 || !r2) continue;
        const auto ps = scorer.score_pair(sc, *r1, *r2, finder, chem, typing);
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%.2f\t%.4f",
                      ki.c_str(), kj.c_str(),
                      ps.lw_class.empty() ? "None" : ps.lw_class.c_str(),
                      ps.score, ps.penalty);
        std::string row = buf;
        const std::string iss = format_issues(ps);
        if (!iss.empty()) { row += '\t'; row += iss; }
        rows.emplace_back(std::move(row));
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Debug: classify a 7-angle suite and print its suiteness.
int classify_suite(const std::vector<std::string>& angles) {
    std::array<double, 7> s{};
    for (int i = 0; i < 7; ++i) s[i] = std::stod(angles[i]);
    pairfinder::scoring::RichardsonClassifier cl(richardson_suites());
    std::printf("%.6f\n", cl.suiteness(s));
    return 0;
}

// Debug: emit "<res_id>\t<alpha>.<zeta>" (suite torsions) for diffing.
int dump_torsions(const std::string& path) {
    using namespace pairfinder;
    auto structure = io::parse_pdb(path);
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);
    const auto tors = scoring::compute_all_torsions(structure);
    auto fmt = [](const std::optional<double>& v) {
        char b[24];
        if (v) std::snprintf(b, sizeof(b), "%.3f", *v);
        else std::snprintf(b, sizeof(b), "None");
        return std::string(b);
    };
    std::vector<std::string> rows;
    for (const auto& [rid, t] : tors)
        rows.push_back(rid + "\t" + fmt(t.alpha) + "\t" + fmt(t.beta) + "\t" + fmt(t.gamma) +
                       "\t" + fmt(t.delta) + "\t" + fmt(t.epsilon) + "\t" + fmt(t.zeta));
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Whole-structure score (Cerny method).
// Emits:
//   "SCORES\t<pairs_score>\t<residues_score>"
//   "PAIR\t<res_a>\t<res_b>\t<lw_class>\t<score>\t<penalty>\t<issues>"  (sorted)
//   "RES\t<res_id>\t<score>"                                             (sorted)
int dump_structure(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    auto& structure = prep.structure;
    auto& chem = pipeline.chem();
    const auto& typing = pipeline.typing();
    auto scored = pipeline.classify(prep);
    const auto selected = pipeline.select(std::move(scored), structure);

    auto scorer = pipeline.make_scorer();
    algorithms::hbond::HBondFinder finder(chem);
    const auto ss = scorer.score_structure(selected, structure, finder, chem, typing);

    std::printf("SCORES\t%.2f\t%.2f\n", ss.pairs_score, ss.residues_score);

    std::vector<std::string> rows;
    for (const auto& ps : ss.pair_scores) {
        const std::string ki = ps.res_id1 < ps.res_id2 ? ps.res_id1 : ps.res_id2;
        const std::string kj = ps.res_id1 < ps.res_id2 ? ps.res_id2 : ps.res_id1;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "PAIR\t%s\t%s\t%s\t%.2f\t%.4f",
                      ki.c_str(), kj.c_str(),
                      ps.lw_class.empty() ? "None" : ps.lw_class.c_str(),
                      ps.score, ps.penalty);
        std::string row = buf;
        const std::string iss = format_issues(ps);
        if (!iss.empty()) { row += '\t'; row += iss; }
        rows.emplace_back(std::move(row));
    }
    for (const auto& rs : ss.residue_scores) {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "RES\t%s\t%.2f", rs.res_id.c_str(), rs.score);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Per-residue sugar pucker: emit "<res_id>\t<delta>\t<pperp>\t<helical>" for every
// scored residue that has both - the full-set delta-vs-Pperp distribution, at C++ speed.
// `helical` = 1 if the residue sits in a HELICAL cWW G-C/A-U pair (the sugar-pucker
// check's context; A-form -> C3'-endo), else 0. (residue_scores are already RNA-only.)
int dump_pucker(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    auto& structure = prep.structure;
    auto& chem = pipeline.chem();
    const auto& typing = pipeline.typing();
    auto scored = pipeline.classify(prep);
    const auto selected = pipeline.select(std::move(scored), structure);

    auto scorer = pipeline.make_scorer();
    algorithms::hbond::HBondFinder finder(chem);
    const auto ss = scorer.score_structure(selected, structure, finder, chem, typing);

    // Residues in a helical cWW canonical (G-C/A-U) pair - the pucker check's context.
    static const std::unordered_set<std::string> kWc = {"G-C", "C-G", "A-U", "U-A"};
    const auto chains = algorithms::RNAChains::from_structure(structure);
    const auto in_helix = algorithms::selection::helix_membership(selected, chains);
    std::unordered_set<std::string> helical_wc_res;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto& sc = selected[i];
        if (!in_helix[i] || sc.lw_class != "cWW") continue;
        const std::string bp = typing.normalize(sc.res_name1) + "-" + typing.normalize(sc.res_name2);
        if (!kWc.count(bp)) continue;
        helical_wc_res.insert(sc.res_id1);
        helical_wc_res.insert(sc.res_id2);
    }

    for (const auto& rs : ss.residue_scores) {
        if (!rs.delta || !rs.pperp) continue;
        const int helical = helical_wc_res.count(rs.res_id) ? 1 : 0;
        std::printf("%s\t%.2f\t%.3f\t%d\n", rs.res_id.c_str(), *rs.delta, *rs.pperp, helical);
    }
    return 0;
}

// Full pipeline through selection: emit the SELECTED pairs as
// "<res_a>\t<res_b>\t<lw_class>" sorted, for diffing.
int dump_pairs(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    auto scored = pipeline.classify(prep);
    const auto selected = pipeline.select(std::move(scored), prep.structure);

    std::vector<std::string> rows;
    for (const auto& sc : selected) {
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s", ki.c_str(), kj.c_str(),
                      sc.lw_class.empty() ? "None" : sc.lw_class.c_str());
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Per-candidate selection disposition: WHY each candidate is / isn't in the
// final selected output, decided by the REAL production selection (no proxy).
// Emits "<res_a>\t<res_b>\t<lw>\t<disposition>" sorted. Dispositions:
//   selected | invalid | stacking | low_score | adjacent | ineligible | competition
int dump_selection(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    auto scored = pipeline.classify(prep);
    const auto chains = algorithms::RNAChains::from_structure(prep.structure);
    const auto disp = algorithms::selection::select_pairs_dispositions(std::move(scored), chains);

    std::vector<std::string> rows;
    for (const auto& d : disp) {
        // LW edges refer to res_id1 and res_id2 in this order. Sorting the
        // residue IDs without swapping the edges reverses asymmetric classes.
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%s", d.res_id1.c_str(), d.res_id2.c_str(),
                      d.lw_class.empty() ? "None" : d.lw_class.c_str(), d.disposition.c_str());
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Full scoring pass (classify + score + O2 correction) over all valid candidates.
// Emits "<res_a>\t<res_b>\t<pair_cat>\t<lw>\t<rmsd>\t<quality>\t<num_hb>\t<strong>\t
//        <display>\t<conf_csv>\t<is_amb>\t<precorr>" sorted by valid base-base candidates.
int dump_scored(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto prep = pipeline.prepare(path);
    const auto scored = pipeline.classify(prep);

    std::vector<std::string> rows;
    for (const auto& sc : scored) {
        if (!sc.validation.is_valid) continue;
        // Canonical key (sorted res_id) for order-independent comparison.
        const std::string ki = sc.res_id1 < sc.res_id2 ? sc.res_id1 : sc.res_id2;
        const std::string kj = sc.res_id1 < sc.res_id2 ? sc.res_id2 : sc.res_id1;
        std::string rmsd_s = "None";
        if (sc.template_rmsd) {
            char b[32];
            if (std::isinf(*sc.template_rmsd)) std::snprintf(b, sizeof(b), "inf");
            else std::snprintf(b, sizeof(b), "%.4f", *sc.template_rmsd);
            rmsd_s = b;
        }
        std::string conf_csv;
        for (std::size_t k = 0; k < sc.lw_class_confidence.size(); ++k) {
            char cb[16];
            std::snprintf(cb, sizeof(cb), "%.3f", sc.lw_class_confidence[k]);
            if (k) conf_csv += ",";
            conf_csv += cb;
        }
        char buf[640];
        std::snprintf(buf, sizeof(buf),
                      "%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%s\t%s\t%d\t%s", ki.c_str(),
                      kj.c_str(), sc.pair_category.c_str(),
                      sc.lw_class.empty() ? "None" : sc.lw_class.c_str(), rmsd_s.c_str(),
                      sc.quality_score, sc.num_hbonds, sc.has_strong_base_hbond ? 1 : 0,
                      sc.lw_class_display.empty() ? "None" : sc.lw_class_display.c_str(),
                      conf_csv.empty() ? "None" : conf_csv.c_str(), sc.is_ambiguous ? 1 : 0,
                      sc.precorr_lw.empty() ? "None" : sc.precorr_lw.c_str());
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Full classification (geometry + template RMSD + ambiguity) per valid candidate.
// Emits "<res_a>\t<res_b>\t<lw_class>\t<rmsd|None|inf>\t<swapped>" sorted.
int dump_classified(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto& chem = pipeline.chem();
    const auto& typing = pipeline.typing();
    const auto& aligner = pipeline.aligner();
    const auto& patterns = pipeline.patterns();
    const auto& qscorer = pipeline.quality_scorer();

    auto structure = io::parse_pdb(path);
    const auto frames = algorithms::calculate_all_frames(structure, pipeline.templates());
    const auto cands = algorithms::find_candidates(structure, frames, typing);
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    std::unordered_map<std::string, const core::ReferenceFrame*> frame_by_id;
    for (const auto& [rid, f] : frames) frame_by_id[rid] = &f;

    algorithms::hbond::HBondFinder finder(chem);

    std::vector<std::string> rows;
    for (const auto& c : cands) {
        if (!c.validation.is_valid) continue;
        std::string ra = c.res_id1, rb = c.res_id2;
        if (rb < ra) std::swap(ra, rb);
        const core::Residue* r1 = res_by_id[ra];
        const core::Residue* r2 = res_by_id[rb];
        const core::ReferenceFrame* f1 = frame_by_id[ra];
        const core::ReferenceFrame* f2 = frame_by_id[rb];
        if (!r1 || !r2 || !f1 || !f2) continue;
        const auto voters = finder.find_face_voters(*r1, *r2);
        const auto strict = finder.find_between(*r1, *r2);
        const auto res = algorithms::classification::classify_pair(*f1, *f2, *r1, *r2, voters,
                                                                   strict, aligner, typing);
        std::string rmsd_s = "None";
        if (res.template_rmsd) {
            const double v = *res.template_rmsd;
            char rb2[32];
            if (std::isinf(v)) std::snprintf(rb2, sizeof(rb2), "inf");
            else std::snprintf(rb2, sizeof(rb2), "%.4f", v);
            rmsd_s = rb2;
        }
        // Confidence: res order is post-swap (lw_class[1] describes cr1).
        const core::Residue* cr1 = res.swapped ? r2 : r1;
        const core::Residue* cr2 = res.swapped ? r1 : r2;
        const std::string sequence = cr1->base_type() + cr2->base_type();
        const auto conf = algorithms::classification::compute_classification_confidence(
            *cr1, *cr2, sequence, res.lw_class, strict, aligner, patterns, typing);
        std::string conf_csv;
        for (std::size_t k = 0; k < conf.confidences.size(); ++k) {
            char cb[16];
            std::snprintf(cb, sizeof(cb), "%.3f", conf.confidences[k]);
            if (k) conf_csv += ",";
            conf_csv += cb;
        }
        const double quality = algorithms::classification::compute_selection_quality(
            c.validation, sequence, strict, res, qscorer, patterns, typing);
        char buf[512];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%.3f", ra.c_str(),
                      rb.c_str(), res.lw_class.c_str(), rmsd_s.c_str(), res.swapped ? 1 : 0,
                      conf.display_label.c_str(), conf_csv.c_str(), conf.is_ambiguous ? 1 : 0,
                      quality);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Classify each valid candidate (canonical res_id order) and emit
// "<res_a>\t<res_b>\t<lw_class>\t<swapped>" sorted, for diffing the LW classifier.
int dump_lwclass(const std::string& path) {
    using namespace pairfinder;
    Pipeline pipeline(cli_pipeline_config());
    auto& chem = pipeline.chem();
    const auto& typing = pipeline.typing();

    auto structure = io::parse_pdb(path);
    const auto frames = algorithms::calculate_all_frames(structure, pipeline.templates());
    const auto cands = algorithms::find_candidates(structure, frames, typing);

    // Switch residues to pipeline base types (DC->C, DT->T) - frames
    // were already computed from the original names.
    for (auto& chain : structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);

    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    std::unordered_map<std::string, const core::ReferenceFrame*> frame_by_id;
    for (const auto& [rid, f] : frames) frame_by_id[rid] = &f;

    algorithms::hbond::HBondFinder finder(chem);

    std::vector<std::string> rows;
    for (const auto& c : cands) {
        if (!c.validation.is_valid) continue;
        // Canonical (res_a < res_b) so the classifier's input order is deterministic.
        std::string ra = c.res_id1, rb = c.res_id2;
        if (rb < ra) std::swap(ra, rb);
        const core::Residue* r1 = res_by_id[ra];
        const core::Residue* r2 = res_by_id[rb];
        const core::ReferenceFrame* f1 = frame_by_id[ra];
        const core::ReferenceFrame* f2 = frame_by_id[rb];
        if (!r1 || !r2 || !f1 || !f2) continue;
        const auto voters = finder.find_face_voters(*r1, *r2);
        const auto res = algorithms::classification::classify_lw_class(*f1, *f2, *r1, *r2,
                                                                       voters, typing);
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%s\t%s\t%s\t%d", ra.c_str(), rb.c_str(),
                      res.lw_class.c_str(), res.swapped ? 1 : 0);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Emit each pair candidate (within 15A) as a sorted TSV row for diffing:
// "<res_i>\t<res_j>\t dorg d_v plane dNN dx dy dz qscore valid".
int dump_candidates(const std::string& path) {
    using namespace pairfinder;
    const auto structure = io::parse_pdb(path);
    const core::NucleotideRegistry registry(nucleotides_config());
    algorithms::BaseTemplateLibrary lib(templates_dir(), registry);
    const auto frames = algorithms::calculate_all_frames(structure, lib);
    const core::BaseTyping typing(ligand_db());
    const auto cands = algorithms::find_candidates(structure, frames, typing);

    std::vector<std::string> rows;
    rows.reserve(cands.size());
    for (const auto& c : cands) {
        const std::string& a = c.res_id1;
        const std::string& b = c.res_id2;
        const std::string& ki = a < b ? a : b;
        const std::string& kj = a < b ? b : a;
        const auto& v = c.validation;
        char buf[640];
        std::snprintf(buf, sizeof(buf),
                      "%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.5f\t%.5f\t%.5f\t%d",
                      ki.c_str(), kj.c_str(), v.dorg, v.d_v, v.plane_angle, v.dNN,
                      v.dir_x, v.dir_y, v.dir_z, v.quality_score,
                      v.is_valid ? 1 : 0);
        rows.emplace_back(buf);
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Run find_between on every residue pair and emit each H-bond as
// "<res_i>\t<res_j>\t<donor_res>\t<donor_atom>\t<acc_res>\t<acc_atom>\t<dist>\t
//  <h_idx>\t<lp_idx>\t<align>" - sorted - for diffing.
int dump_hbonds(const std::string& path) {
    using namespace pairfinder;
    const auto structure = io::parse_pdb(path);
    algorithms::hbond::HBondChemistry chem(ligand_db());
    algorithms::hbond::HBondFinder finder(chem);

    std::vector<const core::Residue*> residues;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) residues.push_back(&res);
    // Canonical order (sorted res_id) so find_between's argument order is
    // deterministic (selection can depend on it).
    std::sort(residues.begin(), residues.end(),
              [](const core::Residue* a, const core::Residue* b) {
                  return a->res_id() < b->res_id();
              });

    std::vector<std::string> rows;
    for (std::size_t i = 0; i < residues.size(); ++i) {
        for (std::size_t j = i + 1; j < residues.size(); ++j) {
            const auto hbonds = finder.find_between(*residues[i], *residues[j]);
            for (const auto& hb : hbonds) {
                const std::string& a = residues[i]->res_id();
                const std::string& b = residues[j]->res_id();
                const std::string& key_i = a < b ? a : b;
                const std::string& key_j = a < b ? b : a;
                char buf[512];
                std::snprintf(buf, sizeof(buf),
                              "%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%d\t%d\t%.4f",
                              key_i.c_str(), key_j.c_str(), hb.donor_res_id.c_str(),
                              hb.donor_atom.c_str(), hb.acceptor_res_id.c_str(),
                              hb.acceptor_atom.c_str(), hb.distance, hb.h_slot_idx,
                              hb.lp_slot_idx, hb.alignment_score);
                rows.emplace_back(buf);
            }
        }
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

// Emit "<res_id>\t<atom>\t%.3f\t%.3f\t%.3f", sorted, for diffing.
int dump_parse(const std::string& path) {
    const auto structure = pairfinder::io::parse_pdb(path);
    std::vector<std::string> rows;
    for (const auto& chain : structure.chains()) {
        for (const auto& res : chain.residues()) {
            for (const auto& atom : res.atoms()) {
                char buf[256];
                std::snprintf(buf, sizeof(buf), "%s\t%s\t%.3f\t%.3f\t%.3f",
                              res.res_id().c_str(), atom.name.c_str(),
                              atom.coords.x, atom.coords.y, atom.coords.z);
                rows.emplace_back(buf);
            }
        }
    }
    std::sort(rows.begin(), rows.end());
    for (const auto& r : rows) std::cout << r << '\n';
    return 0;
}

}  // namespace pfcli
