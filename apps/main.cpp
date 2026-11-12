/**
 * @file main.cpp
 * @brief rubric CLI entry point.
 *
 * Production path: resolve the input (a 4-char PDB id fetched from RCSB, or a
 * local .cif/.pdb(.gz)), run the find + classify + score pipeline, and emit the
 * JSON report. The hidden `dump-*` debugging subcommands live in
 * differential_commands.cpp; resource-path resolution lives in cli_support.cpp.
 */
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/rigid_body.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/geometry/vector3d.hpp>
#include <pairfinder/pipeline.hpp>
#include <pairfinder/scoring/hbond_angles.hpp>
#include <pairfinder/scoring/pair_parameters.hpp>
#include <pairfinder/scoring/scorer.hpp>

#include "support/cli_support.hpp"
#include "differential/differential_commands.hpp"

namespace {
#ifndef RUBRIC_VERSION
#define RUBRIC_VERSION "unknown"
#endif
constexpr const char* kVersion = RUBRIC_VERSION;  // from project(VERSION) in CMakeLists.txt

// Maps raw parameter keys to human-readable report group names.
const std::unordered_map<std::string, std::string> kReportGroup = {
    {"shear",                 "Pair In-Plane Displacement"},
    {"stretch",               "Pair In-Plane Displacement"},
    {"stagger",               "Pair Non-Coplanarity"},
    {"buckle",                "Pair Non-Coplanarity"},
    {"propeller",             "Propeller Twist"},
    {"opening",               "Pair Opening"},
    {"distance",              "H-Bond Distance"},
    {"hbond_angles",          "H-Bond Angles"},
    {"incorrect_hbond_count", "Incorrect H-Bond Count"},
    {"clash",                 "Steric Clash"},
};

// ---- full `find` pipeline -> JSON ----
using json = nlohmann::json;

double round_to(double v, int n) {
    const double f = std::pow(10.0, n);
    return std::round(v * f) / f;
}

json vec3_json(const pairfinder::geometry::Vector3d& v) {
    return json::array({round_to(v.x, 4), round_to(v.y, 4), round_to(v.z, 4)});
}

// Readable pair label from two res_ids ("<chain>-<resname>-<seq+icode>"):
// same chain -> "B:U9-A33", cross-chain -> "A:U9-B:A33". Uses the parsed
// residue name, so modified bases stay as-is (e.g. "B:CBR12-G40").
std::string pair_label(const std::string& id1, const std::string& id2) {
    auto split = [](const std::string& id) {
        const auto a = id.find('-');
        const auto b = a == std::string::npos ? a : id.find('-', a + 1);  // seq may be negative
        if (b == std::string::npos) return std::pair<std::string, std::string>{"", id};
        return std::pair<std::string, std::string>{
            id.substr(0, a), id.substr(a + 1, b - a - 1) + id.substr(b + 1)};
    };
    const auto [c1, r1] = split(id1);
    const auto [c2, r2] = split(id2);
    if (c1 == c2) return (c1.empty() ? "" : c1 + ":") + r1 + "-" + r2;
    return c1 + ":" + r1 + "-" + c2 + ":" + r2;
}

// Filename stem with structure extensions stripped (.gz, then .cif/.pdb/.ent).
std::string clean_pdb_id(const std::string& path) {
    std::string s = std::filesystem::path(path).filename().string();
    auto ends_with = [&](const std::string& ext) {
        if (s.size() < ext.size()) return false;
        std::string tail = s.substr(s.size() - ext.size());
        for (char& c : tail) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return tail == ext;
    };
    if (ends_with(".gz")) s.erase(s.size() - 3);
    for (const char* ext : {".cif", ".pdb", ".ent"})
        if (ends_with(ext)) { s.erase(s.size() - 4); break; }
    return s;
}

json frame_json(const pairfinder::core::ReferenceFrame& f) {
    return json{{"origin", vec3_json(f.origin)},
                {"x_axis", vec3_json(f.x_axis())},
                {"y_axis", vec3_json(f.y_axis())},
                {"z_axis", vec3_json(f.z_axis())}};
}

// ---- --pml: a standalone PyMOL script reproducing the plug-in's display ------
// Mirrors integrations/pymol/rubric_pymol.py (score_to_rgb, _build_overlay,
// _build_pucker) so the script and the plug-in draw the same picture without
// needing the plug-in or the binary at view time.

constexpr double kReviewMax = 80.0;          // priority review: score < 80
constexpr double kSecondaryG0 = 0.867;       // secondary review green channel near 100
constexpr double kSecondaryG1 = 0.45;        //                          ... near 80
constexpr std::size_t kResiPerLine = 60;     // residues per selection line

// PyMOL object/selection name: letters, digits and '_' only.
std::string pymol_name(const std::string& s) {
    std::string n;
    for (char c : s) n += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
    return n.empty() ? "rubric_model" : n;
}

// res_id "<chain>-<resname>-<seq+icode>" -> (chain, resi); seq may be negative.
std::pair<std::string, std::string> chain_resi(const std::string& res_id) {
    const auto a = res_id.find('-');
    const auto b = a == std::string::npos ? a : res_id.find('-', a + 1);
    if (b == std::string::npos) return {"", res_id};
    std::string seq = res_id.substr(b + 1);
    if (!seq.empty() && seq[0] == '-') seq = "\\" + seq;  // PyMOL escapes a leading '-'
    return {res_id.substr(0, a), seq};
}

// "(chain A and resi 1+2+3)" clauses, at most kResiPerLine residues each, so a
// ribosome-sized selection is built over several short lines.
std::vector<std::string> resi_clauses(const std::vector<std::string>& res_ids) {
    std::vector<std::string> order;
    std::unordered_map<std::string, std::vector<std::string>> by_chain;
    for (const auto& id : res_ids) {
        auto [ch, resi] = chain_resi(id);
        if (!by_chain.count(ch)) order.push_back(ch);
        auto& v = by_chain[ch];
        if (std::find(v.begin(), v.end(), resi) == v.end()) v.push_back(resi);
    }
    std::vector<std::string> out;
    for (const auto& ch : order) {
        const auto& v = by_chain[ch];
        for (std::size_t i = 0; i < v.size(); i += kResiPerLine) {
            std::string c = "(" + (ch.empty() ? std::string() : "chain " + ch + " and ") + "resi ";
            for (std::size_t j = i; j < std::min(v.size(), i + kResiPerLine); ++j)
                c += (j > i ? "+" : "") + v[j];
            out.push_back(c + ")");
        }
    }
    return out;
}

// Write `select name, ...`: one line for a single clause, else extended over
// as many lines as needed (empty -> "none").
void write_selection(std::ostream& o, const std::string& name, const std::string& scope,
                     const std::vector<std::string>& res_ids) {
    const auto clauses = resi_clauses(res_ids);
    if (clauses.empty()) { o << "select " << name << ", none\n"; return; }
    o << "select " << name << ", " << scope << " and " << clauses[0] << "\n";
    for (std::size_t i = 1; i < clauses.size(); ++i)
        o << "select " << name << ", " << name << " or (" << scope << " and " << clauses[i] << ")\n";
}

// Python string literal for the embedded block (chain/resi are plain tokens,
// but escape defensively).
std::string py_str(const std::string& s) {
    std::string q = "\"";
    for (char c : s) q += (c == '"' || c == '\\') ? std::string("\\") + c : std::string(1, c);
    return q + "\"";
}

std::string fmt3(double v) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.3f", v);
    return buf;
}

// Emit the PyMOL script for a scored report. `structure_path` is the file the
// engine read; the script loads that same file, so selections match its chains.
void write_pml(const json& out, const std::string& structure_path, const std::string& pml_path) {
    const std::string obj = pymol_name(out.value("pdb_id", std::string()));
    const std::string abs_path = std::filesystem::absolute(structure_path).string();

    // Flagged pairs (score < 100), worst first; each residue takes the lowest
    // score of the flagged pairs it belongs to (as in the plug-in overlay).
    std::vector<const json*> flagged;
    for (const auto& p : out["pairs"])
        if (p.contains("score") && p["score"].is_number() && p["score"].get<double>() < 100.0)
            flagged.push_back(&p);
    std::stable_sort(flagged.begin(), flagged.end(), [](const json* a, const json* b) {
        return (*a)["score"].get<double>() < (*b)["score"].get<double>();
    });
    std::vector<std::string> res_order;
    std::unordered_map<std::string, double> lo_score;
    for (const json* p : flagged) {
        const double s = (*p)["score"].get<double>();
        for (const char* k : {"res_id1", "res_id2"}) {
            const std::string id = (*p)[k].get<std::string>();
            auto it = lo_score.find(id);
            if (it == lo_score.end()) { lo_score[id] = s; res_order.push_back(id); }
            else it->second = std::min(it->second, s);
        }
    }
    // Position on the plug-in's yellow(0) -> red(1) spectrum, whose green
    // channel is 1 - b: priority review -> 1 (red); secondary review -> the
    // yellow (near 100) to orange (near 80) stretch of the same ramp.
    auto ramp_b = [](double s) {
        if (s < kReviewMax) return 1.0;
        const double t = (100.0 - s) / (100.0 - kReviewMax);
        return 1.0 - (kSecondaryG0 + t * (kSecondaryG1 - kSecondaryG0));
    };

    std::vector<std::string> pucker_res;
    if (out.contains("pucker_outliers"))
        for (const auto& u : out["pucker_outliers"]) pucker_res.push_back(u["res_id"].get<std::string>());

    std::ostringstream o;
    auto num_or = [](const json& v, const char* fmt) -> std::string {
        if (!v.is_number()) return "n/a";
        char buf[32];
        std::snprintf(buf, sizeof buf, fmt, v.get<double>());
        return buf;
    };
    const json ts = out.value("tier_summary", json::object());
    o << "# RUBRIC " << kVersion << " - PyMOL view of " << out.value("pdb_id", std::string()) << "\n"
      << "# structure: " << abs_path << "\n"
      << "# base-pair score " << num_or(out.value("pairs_score", json()), "%.2f")
      << "  |  backbone suiteness " << num_or(out.value("backbone_suiteness", json()), "%.3f")
      << "  |  pairs: " << ts.value("Preferred", 0) << " preferred, "
      << ts.value("Secondary review", 0) << " secondary review, "
      << ts.value("Priority review", 0) << " priority review\n"
      << "# colors: gray = preferred | yellow -> orange = secondary review |"
         " red = priority review | cyan ribose = pucker outlier (delta vs Pperp)\n"
      << "# Open with: pymol " << std::filesystem::path(pml_path).filename().string()
      << "   (or @file.pml inside PyMOL). Coordinates are never modified;\n"
      << "# re-run RUBRIC with --pml after editing the model to refresh this view.\n\n"
      << "set ignore_case_chain, 0\n"
      << "/cmd.load(r\"\"\"" << abs_path << "\"\"\", \"" << obj << "\")\n"
      << "hide everything, " << obj << "\n"
      << "show cartoon, " << obj << "\n"
      << "color gray70, " << obj << "\n\n";

    o << "# ---- flagged base pairs (score < 100), worst first ----\n";
    if (flagged.empty()) o << "# none\n";
    for (std::size_t i = 0; i < flagged.size(); ++i) {
        const json& p = *flagged[i];
        std::string lw = "-", issues;
        if (p.contains("classification") && p["classification"].is_object())
            lw = p["classification"].value("lw_class", std::string("-"));
        for (const auto& is : p.value("issues", json::array()))
            issues += (issues.empty() ? "" : ", ") + is.get<std::string>();
        char buf[64];
        std::snprintf(buf, sizeof buf, "# %4zu  %6.2f  ", i + 1, p["score"].get<double>());
        o << buf << p.value("tier", std::string()) << "  " << p.value("pair", std::string())
          << "  " << lw << (issues.empty() ? "" : "  [" + issues + "]") << "\n";
    }
    if (!pucker_res.empty()) {
        o << "# ---- ribose pucker outliers (delta and Pperp disagree) ----\n";
        for (const auto& u : out["pucker_outliers"])
            o << "#   " << u.value("res_id", std::string()) << "  delta "
              << num_or(u.value("delta", json()), "%.1f") << " ("
              << u.value("pucker_by_delta", std::string("?")) << ") vs Pperp "
              << num_or(u.value("pperp", json()), "%.2f") << " A ("
              << u.value("pucker_by_pperp", std::string("?")) << ")\n";
    }
    o << "\n";

    const std::string nuc = "(" + obj + ") and (polymer.nucleic or byres ((" + obj + ") and name C1'))";
    if (!res_order.empty()) {
        // As in the plug-in's _build_overlay: stamp each flagged residue's ramp
        // position into the B-factor of a COPY, drop the rest, and color with one
        // spectrum call. One pass over the atoms, so a ribosome opens quickly; the
        // model itself is untouched.
        o << "# ---- flagged pairs as sticks on a copy, colored by pair score ----\n"
          << "python\n"
          << "_rubric_b = {\n";
        for (std::size_t i = 0; i < res_order.size(); ++i) {  // 8 per line: PyMOL echoes each line
            const std::string& id = res_order[i];
            const auto a = id.find('-');
            const auto b = a == std::string::npos ? a : id.find('-', a + 1);
            const std::string ch = b == std::string::npos ? "" : id.substr(0, a);
            const std::string seq = b == std::string::npos ? id : id.substr(b + 1);
            o << (i % 8 == 0 ? "    " : " ") << "(" << py_str(ch) << ", " << py_str(seq)
              << "): " << fmt3(ramp_b(lo_score[id])) << ","
              << (i % 8 == 7 || i + 1 == res_order.size() ? "\n" : "");
        }
        o << "}\n"
          << "cmd.create(\"rubric_pairs\", " << py_str(nuc) << ")\n"
          << "cmd.alter(\"rubric_pairs\", \"b = _rubric_b.get((chain, resi), -1.0)\",\n"
          << "          space={\"_rubric_b\": _rubric_b})\n"
          << "cmd.remove(\"rubric_pairs and b < -0.5\")\n"
          << "del _rubric_b\n"
          << "python end\n"
          << "hide everything, rubric_pairs\n"
          << "show sticks, rubric_pairs\n"
          << "set stick_radius, 0.18, rubric_pairs\n"
          << "spectrum b, yellow red, rubric_pairs, minimum=0, maximum=1\n"
          << "select rubric_priority, rubric_pairs and b > 0.999\n"
          << "select rubric_secondary, rubric_pairs and b < 0.999\n\n";

        // One named selection per priority-review pair (on the model), in worklist order.
        int k = 0;
        for (const json* p : flagged) {
            if ((*p)["score"].get<double>() >= kReviewMax) continue;
            char nm[32];
            std::snprintf(nm, sizeof nm, "rubric_pr_%03d", ++k);
            write_selection(o, nm, obj, {(*p)["res_id1"].get<std::string>(),
                                         (*p)["res_id2"].get<std::string>()});
        }
        if (k > 0) o << "group rubric_worklist, rubric_pr_*\n\n";
    }
    if (!pucker_res.empty()) {
        o << "# ---- pucker outliers: gray sticks with a cyan ribose ring ----\n";
        write_selection(o, "rubric_pucker_outliers", obj, pucker_res);
        o << "create rubric_pucker, (" << obj << ") and rubric_pucker_outliers\n"
          << "hide everything, rubric_pucker\n"
          << "show sticks, rubric_pucker\n"
          << "set stick_radius, 0.30, rubric_pucker\n"
          << "color gray60, rubric_pucker\n"
          << "set_color rubric_pucker_col, [0.000, 0.750, 0.850]\n"
          << "color rubric_pucker_col, rubric_pucker and name C1'+C2'+C3'+C4'+O4'\n\n";
    }
    o << "deselect\nzoom " << obj << "\n";

    std::filesystem::path pp(pml_path);
    std::error_code ec;
    if (pp.has_parent_path()) std::filesystem::create_directories(pp.parent_path(), ec);
    std::ofstream(pp) << o.str();
    std::fprintf(stderr, "Wrote %s\n", pp.string().c_str());
}

// Run the full pipeline and emit the report JSON (default + --details).
int find_json(const std::string& path, bool score, bool details, bool emit_classified,
              const std::string& out_path, const std::string& pml_path) {
    using namespace pairfinder;
    // Per-stage wall-clock when PAIRFINDER_PROFILE is set (-> stderr).
    const bool prof = std::getenv("RUBRIC_PROFILE") != nullptr;
    auto t0 = std::chrono::steady_clock::now();
    auto lap = [&](const char* name) {
        if (!prof) return;
        const auto now = std::chrono::steady_clock::now();
        std::fprintf(stderr, "[profile] %-22s %8.3f s\n", name,
                     std::chrono::duration<double>(now - t0).count());
        t0 = now;
    };

    Pipeline pipeline(pfcli::cli_pipeline_config());
    lap("load resources");
    auto prep = pipeline.prepare(path);
    lap("prepare (parse + frames + pipeline_view)");
    auto scored = pipeline.classify(prep);
    lap("score_candidates (find+classify)");

    int candidates_total = static_cast<int>(scored.size());
    int candidates_valid = 0;
    for (const auto& sc : scored)
        if (sc.validation.is_valid) ++candidates_valid;

    // Optionally capture every valid classified candidate (pre-selection): whether
    // a pair was found + classified at all, independent of helix-priority selection.
    json classified = json::array();
    if (emit_classified) {
        for (const auto& sc : scored) {
            if (!sc.validation.is_valid) continue;
            if (sc.lw_class.empty() && sc.lw_class_display.empty()) continue;
            classified.push_back(json{
                {"res_id1", sc.res_id1},
                {"res_id2", sc.res_id2},
                {"lw_class", sc.lw_class_display.empty() ? sc.lw_class : sc.lw_class_display},
                {"pair_category", sc.pair_category},
                {"is_ambiguous", sc.is_ambiguous}});
        }
    }

    const auto selected = pipeline.select(std::move(scored), prep.structure);
    lap("rna_chains + selection");

    // Empirical scoring (per-pair score/issues + structure aggregates).
    scoring::StructureScore ss;
    std::unordered_map<std::string, const scoring::PairScore*> score_by_key;
    if (score) {
        auto scorer = pipeline.make_scorer();
        algorithms::hbond::HBondFinder sfinder(pipeline.chem());
        ss = scorer.score_structure(selected, prep.structure, sfinder, pipeline.chem(),
                                    pipeline.typing(), /*collect_backbone=*/details);
        for (const auto& ps : ss.pair_scores)
            score_by_key[ps.res_id1 + '\x01' + ps.res_id2] = &ps;
        lap("score_structure");
    }

    json out;
    out["pdb_id"] = clean_pdb_id(path);
    out["candidates_total"] = candidates_total;
    out["candidates_valid"] = candidates_valid;
    out["n_pairs"] = static_cast<int>(selected.size());

    // For --details per-pair H-bond geometry: residue lookup + a finder to
    // re-run between paired residues (the candidate keeps only counts).
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    if (details)
        for (const auto& chain : prep.structure.chains())
            for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;
    algorithms::hbond::HBondFinder hb_finder(pipeline.chem());

    json pairs = json::array();
    for (const auto& sc : selected) {
        json p;
        p["res_id1"] = sc.res_id1;
        p["res_id2"] = sc.res_id2;
        p["pair"] = pair_label(sc.res_id1, sc.res_id2);
        p["sequence"] = sc.res_name1 + sc.res_name2;
        // Normalized parent base-pair type ("C-G"), matching the scoring-table
        // key and pairs.csv bp_type (modified bases -> parent, e.g. CBR -> C).
        p["bp_type"] = pipeline.typing().normalize(sc.res_name1) + "-" +
                       pipeline.typing().normalize(sc.res_name2);
        p["pair_category"] = sc.pair_category;
        p["num_hbonds"] = sc.num_hbonds;
        p["hbond_categories"] = sc.hbond_categories;

        if (sc.lw_class_display.empty() && sc.lw_class.empty()) {
            p["classification"] = nullptr;
        } else {
            p["classification"] = json{
                {"lw_class", sc.lw_class_display.empty() ? sc.lw_class : sc.lw_class_display},
                {"is_ambiguous", sc.is_ambiguous}};
        }

        const auto& v = sc.validation;
        p["geometry"] = json{{"dorg", round_to(v.dorg, 4)},
                             {"d_v", round_to(v.d_v, 4)},
                             {"plane_angle", round_to(v.plane_angle, 4)},
                             {"dNN", round_to(v.dNN, 4)}};

        if (sc.frame1 && sc.frame2) {
            const auto rb = scoring::compute_pair_parameters(*sc.frame1, *sc.frame2);
            p["rigid_body"] = json{{"shear", round_to(rb.shear, 4)},
                                   {"stretch", round_to(rb.stretch, 4)},
                                   {"stagger", round_to(rb.stagger, 4)},
                                   {"buckle", round_to(rb.buckle, 4)},
                                   {"propeller", round_to(rb.propeller, 4)},
                                   {"opening", round_to(rb.opening, 4)}};
            p["frames"] = json{{"res1", frame_json(*sc.frame1)},
                               {"res2", frame_json(*sc.frame2)}};
        }

        if (details) {
            if (!sc.template_rmsd) p["template_rmsd"] = nullptr;
            else if (std::isinf(*sc.template_rmsd)) p["template_rmsd"] = "Infinity";
            else p["template_rmsd"] = round_to(*sc.template_rmsd, 4);
            p["quality_score"] = round_to(sc.quality_score, 4);

            // Per-pair H-bond geometry: re-run the finder between the paired
            // residues and emit every bond (distance + donor/acceptor angle +
            // slot alignment), so the per-pair feature tables can be
            // regenerated from the engine alone.
            if (scoring::Scorer::is_scorable(sc, pipeline.typing())) {
                const auto r1 = res_by_id.find(sc.res_id1);
                const auto r2 = res_by_id.find(sc.res_id2);
                if (r1 != res_by_id.end() && r2 != res_by_id.end()) {
                    const core::Residue& res1 = *r1->second;
                    const core::Residue& res2 = *r2->second;
                    json hbonds = json::array();
                    for (const auto& hb : hb_finder.find_between(res1, res2)) {
                        const core::Residue& dres = (hb.donor_res_id == sc.res_id1) ? res1 : res2;
                        const core::Residue& ares = (hb.donor_res_id == sc.res_id1) ? res2 : res1;
                        const core::Atom* da = dres.get_atom(hb.donor_atom);
                        const core::Atom* aa = ares.get_atom(hb.acceptor_atom);
                        std::optional<double> d_ang, a_ang;
                        if (da && aa) {
                            d_ang = scoring::hbond_angle(dres, hb.donor_atom, da->coords,
                                                         aa->coords, pipeline.chem());
                            a_ang = scoring::hbond_angle(ares, hb.acceptor_atom, aa->coords,
                                                         da->coords, pipeline.chem());
                        }
                        hbonds.push_back(json{
                            {"donor_atom", hb.donor_atom},
                            {"acceptor_atom", hb.acceptor_atom},
                            {"distance", round_to(hb.distance, 4)},
                            {"donor_angle", d_ang ? json(round_to(*d_ang, 4)) : json(nullptr)},
                            {"acceptor_angle", a_ang ? json(round_to(*a_ang, 4)) : json(nullptr)},
                            {"slot_alignment", round_to(hb.alignment_score, 4)}});
                    }
                    p["hbonds"] = hbonds;
                }
            }
        }

        if (score) {
            const auto it = score_by_key.find(sc.res_id1 + '\x01' + sc.res_id2);
            if (it != score_by_key.end()) {
                p["score"] = round_to(it->second->score, 2);
                p["tier"] = it->second->tier();
                // Issues as a flat list of key strings.
                json iss = json::array();
                // Per-issue Cerny provenance (ProSco / Z' / tier) for downstream
                // highlighting; same order as "issues" (descending penalty).
                json issue_rows = json::array();
                for (const auto& ip : it->second->issues) {
                    iss.push_back(ip.issue);
                    issue_rows.push_back(json{
                        {"issue", ip.issue},
                        {"severity", round_to(ip.severity, 4)},
                        {"prosco", ip.prosco ? json(round_to(*ip.prosco, 4)) : json(nullptr)},
                        {"zprime", ip.zprime ? json(round_to(*ip.zprime, 4)) : json(nullptr)},
                        {"tier", ip.tier}});
                }
                p["issues"] = std::move(iss);
                p["issue_details"] = std::move(issue_rows);
                // Raw overlap behind the `clash` severity: ungated and
                // unsaturated, so weight fitting can recompute severity.
                if (details) {
                    const auto& ov = it->second->clash_overlap;
                    p["clash_overlap"] = ov ? json(round_to(*ov, 4)) : json(nullptr);
                }
                // Report-only "unusual sugar pucker" (helical cWW C2'-endo).
                if (!it->second->unusual_pucker.empty()) {
                    json ups = json::array();
                    for (const auto& u : it->second->unusual_pucker) {
                        ups.push_back(json{
                            {"res_id", u.res_id},
                            {"pucker", "C2'-endo"},
                            {"expected", "C3'-endo"},
                            {"delta", u.delta ? json(round_to(*u.delta, 1)) : json(nullptr)},
                            {"pperp", u.pperp ? json(round_to(*u.pperp, 2)) : json(nullptr)},
                            {"pperp_confirms", u.pperp_confirms}});
                    }
                    p["unusual_pucker"] = std::move(ups);
                }
            } else {
                p["score"] = nullptr;
                p["tier"] = nullptr;
                p["issues"] = json::array();
                p["issue_details"] = json::array();
            }
        }
        pairs.push_back(std::move(p));
    }
    out["pairs"] = std::move(pairs);
    if (emit_classified) out["classified"] = std::move(classified);

    if (score) {
        // There is deliberately no combined "overall_score": pair-geometry and
        // backbone scores live on different scales / measure different axes, and
        // an "overall" label overclaims whole-structure quality (clashscore,
        // bond geometry, etc. are not considered). Reported separately below.
        // null, not 0.0, when nothing was scored - unmeasured, not worst.
        // Unscorable pairs (DNA, ambiguous) leave n_pairs > 0 with no scores,
        // and an empty mean is 0.0, which would pin the structure to the
        // bottom of every correlation. Mirrors backbone_suiteness.
        if (ss.pair_scores.empty())
            out["pairs_score"] = nullptr;
        else
            out["pairs_score"] = round_to(ss.pairs_score, 2);
        // Mean Richardson suiteness as a 0-1 fraction (wwPDB RNA-suiteness
        // convention). Internally residues_score is 0-100, so scale down here.
        // null when no residue could be scored: that is "unmeasured", not
        // "worst possible backbone", and reporting 0.0 pins such structures to
        // the bottom of every correlation. Mirrors pairs_score for n_pairs == 0.
        if (ss.residue_scores.empty())
            out["backbone_suiteness"] = nullptr;
        else
            out["backbone_suiteness"] = round_to(ss.residues_score / 100.0, 3);
        // Pair triage tier counts (Preferred/Secondary review/Priority review).
        // Always all three keys, so consumers never have to test for absence.
        const auto ts = ss.tier_summary();
        const json tier_summary = json{
            {"Preferred", ts.preferred},
            {"Secondary review", ts.acceptable},
            {"Priority review", ts.review}};
        out["tier_summary"] = tier_summary;
        // Whole-structure RNA-RNA steric clashes (report-only; includes
        // residues that are not a base pair). Protein contacts are excluded by
        // design - see scoring/clashes.hpp.
        json clashes = json::array();
        for (const auto& c : ss.clashes) {
            clashes.push_back(json{{"res_id1", c.res_id1}, {"atom1", c.atom1},
                                   {"res_id2", c.res_id2}, {"atom2", c.atom2},
                                   {"distance", round_to(c.distance, 2)},
                                   {"overlap", round_to(c.overlap, 2)}});
        }
        out["n_clashes"] = static_cast<int>(ss.clashes.size());
        out["clashes"] = std::move(clashes);
        // Per-residue backbone recommendations (flagged residues only): which
        // torsions to correct toward which conformer. Suiteness stays the score.
        json backbone_residues = json::array();
        for (const auto& rs : ss.residue_scores) {
            if (!rs.recommendation) continue;
            const auto& rec = *rs.recommendation;
            json devs = json::array();
            for (const auto& d : rec.deviations) {
                json jd = {{"angle", d.angle}, {"value", round_to(d.value, 1)},
                           {"target", round_to(d.target, 1)}, {"gap", round_to(d.gap, 1)}};
                if (d.prosco) jd["prosco"] = round_to(*d.prosco, 2);
                if (d.conformer_center)
                    jd["conformer_center"] = round_to(*d.conformer_center, 1);
                devs.push_back(std::move(jd));
            }
            backbone_residues.push_back(json{
                {"res_id", rs.res_id},
                {"suiteness", round_to(rs.suiteness, 3)},
                {"tier", rec.tier},
                {"target_conformer", rec.target_conformer},
                {"deviations", std::move(devs)}});
        }
        out["backbone_residues"] = std::move(backbone_residues);

        // Full per-residue backbone torsions (--details only): one row per
        // classifiable residue, so the full torsion table can be regenerated
        // from the engine alone.
        if (details) {
            auto ang = [](const std::optional<double>& v) {
                return v ? json(round_to(*v, 4)) : json(nullptr);
            };
            json bt = json::array();
            for (const auto& r : ss.backbone_torsions) {
                bt.push_back(json{
                    {"res_id", r.res_id},
                    {"base_type", r.base_type},
                    {"alpha", ang(r.alpha)}, {"beta", ang(r.beta)}, {"gamma", ang(r.gamma)},
                    {"delta", ang(r.delta)}, {"epsilon", ang(r.epsilon)}, {"zeta", ang(r.zeta)},
                    {"conformer", r.conformer},
                    {"suiteness", round_to(r.suiteness, 3)},
                    {"suite_bin", r.suite_bin.empty() ? json(nullptr) : json(r.suite_bin)},
                    {"is_outlier", r.is_outlier},
                    {"scored", r.scored}});
            }
            out["backbone_torsions"] = std::move(bt);
        }
        // Report-only ribose pucker OUTLIERS (delta disagrees with Pperp), structure-level.
        json pucker_outliers = json::array();
        for (const auto& rs : ss.residue_scores) {
            if (!rs.pucker_outlier) continue;
            const auto& o = *rs.pucker_outlier;
            pucker_outliers.push_back(json{
                {"res_id", rs.res_id},
                {"pperp", round_to(o.pperp, 2)},
                {"delta", o.delta ? json(round_to(*o.delta, 1)) : json(nullptr)},
                {"pucker_by_pperp", o.pucker_by_pperp},
                {"pucker_by_delta", o.pucker_by_delta}});
        }
        if (!pucker_outliers.empty()) out["pucker_outliers"] = std::move(pucker_outliers);
        // Report-only: delta outside both sieve bands, so no pucker from delta.
        json pucker_unclassified = json::array();
        for (const auto& rs : ss.residue_scores) {
            if (!rs.pucker_unclassified) continue;
            const auto& u = *rs.pucker_unclassified;
            pucker_unclassified.push_back(json{
                {"res_id", rs.res_id},
                {"delta", round_to(u.delta, 1)},
                {"pperp", u.pperp ? json(round_to(*u.pperp, 2)) : json(nullptr)},
                {"pucker_by_pperp", u.pperp ? json(u.pucker_by_pperp) : json(nullptr)}});
        }
        if (!pucker_unclassified.empty())
            out["pucker_unclassified"] = std::move(pucker_unclassified);
        if (details) {
            // Group by human-readable category names.
            json summary = json::object();
            for (const auto& ps : ss.pair_scores) {
                for (const auto& ip : ps.issues) {
                    const auto git = kReportGroup.find(ip.issue);
                    const std::string& grp = (git != kReportGroup.end()) ? git->second : ip.issue;
                    summary[grp] = summary.value(grp, 0) + 1;
                }
            }
            out["score_details"] = json{
                {"pairs_score", ss.pair_scores.empty()
                                   ? json(nullptr)
                                   : json(round_to(ss.pairs_score, 2))},
                {"backbone_suiteness", ss.residue_scores.empty()
                                           ? json(nullptr)
                                           : json(round_to(ss.residues_score / 100.0, 3))},
                {"tier_summary", tier_summary},
                {"issue_summary", summary},
                {"n_scored_pairs", static_cast<int>(ss.pair_scores.size())},
                {"n_residues_scored", static_cast<int>(ss.residue_scores.size())},
                {"n_skipped_pairs", ss.skipped_pairs}};
        }
    }

    const std::string text = out.dump(2);
    if (out_path.empty()) {
        std::cout << text << '\n';
    } else {
        // The output path is used exactly as the user gives it: a bare name
        // lands in the current directory, a relative or absolute path is honored
        // as-is. Any parent directories in the path are created if missing.
        std::filesystem::path outp(out_path);
        std::error_code out_ec;
        if (outp.has_parent_path())
            std::filesystem::create_directories(outp.parent_path(), out_ec);
        std::ofstream(outp) << text << '\n';
        std::fprintf(stderr, "Wrote %s\n", outp.string().c_str());
    }
    if (!pml_path.empty()) write_pml(out, path, pml_path);
    lap("json build + output");
    // find_json is main's last act, so end the process here: the output is flushed,
    // and tearing down every residue, cache and table one by one only to have the
    // OS reclaim the memory anyway costs ~4% of a large structure's run.
    std::cout.flush();
    std::fflush(nullptr);
    std::_Exit(0);
}

// ---- PDB-ID auto-fetch (input acquisition; upstream of the pipeline) ---------

// Cache dir for downloaded structures: --cache-dir, else $RUBRIC_CACHE_DIR,
// else ~/.cache/rubric (%USERPROFILE%\.cache\rubric on Windows), else <tmp>/rubric.
// Persistent + per-user.
std::filesystem::path cache_dir(const std::string& override_dir) {
    if (!override_dir.empty()) return override_dir;
    if (const char* e = std::getenv("RUBRIC_CACHE_DIR")) return e;
    if (const char* h = std::getenv("HOME"))
        return std::filesystem::path(h) / ".cache" / "rubric";
    if (const char* u = std::getenv("USERPROFILE"))
        return std::filesystem::path(u) / ".cache" / "rubric";
    return std::filesystem::temp_directory_path() / "rubric";
}

// A bare PDB id: short, all-alphanumeric (e.g. 6V3A, 1GID). Distinguishes an id
// from a file path (which has an extension and/or a directory separator).
bool looks_like_pdb_id(const std::string& s) {
    if (s.empty() || s.size() > 8) return false;
    for (char c : s)
        if (!std::isalnum(static_cast<unsigned char>(c))) return false;
    return true;
}

// Ensure a cached structure file for <ID> exists, downloading from RCSB if
// needed. Tries formats in priority order (gzipped mmCIF, then plain mmCIF) so a
// transient/missing variant falls back to the next. Returns the path, or empty
// on failure / when downloads are disabled.
std::filesystem::path fetch_structure(const std::string& id_raw,
                                      const std::filesystem::path& cdir, bool allow_download) {
    std::string id = id_raw;
    for (char& c : id) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    std::error_code ec;
    std::filesystem::create_directories(cdir, ec);

    // mmCIF is the canonical archive format (exists for every released entry,
    // incl. large structures with no legacy .pdb); .cif.gz first, .cif fallback.
    const std::array<const char*, 2> suffixes = {".cif.gz", ".cif"};

    // Cache hit: return the first already-downloaded variant.
    for (const char* suf : suffixes) {
        const std::filesystem::path p = cdir / (id + suf);
        if (std::filesystem::exists(p)) return p;
    }
    if (!allow_download) return {};

    // Download: try each variant until one succeeds.
    for (const char* suf : suffixes) {
        const std::filesystem::path dest = cdir / (id + suf);
        const std::filesystem::path part = dest.string() + ".part";
        const std::string url = "https://files.rcsb.org/download/" + id + suf;
        std::fprintf(stderr, "Fetching %s -> %s\n", url.c_str(), dest.string().c_str());
#if defined(_WIN32)
        // cmd.exe does not treat single quotes as quoting; curl ships with Windows 10+.
        const std::string cmd = "curl -fsSL -o \"" + part.string() + "\" \"" + url + "\"";
#else
        const std::string cmd = "curl -fsSL -o '" + part.string() + "' '" + url + "'";
#endif
        if (std::system(cmd.c_str()) == 0) {
            std::filesystem::rename(part, dest, ec);
            if (!ec) return dest;
        }
        std::filesystem::remove(part, ec);
    }
    return {};
}

// Resolve the CLI argument to an actual structure file: existing path used as-is;
// otherwise treated as a PDB id and fetched into the cache.
std::filesystem::path resolve_input(const std::string& arg,
                                    const std::filesystem::path& cdir, bool allow_download) {
    if (std::filesystem::exists(arg)) return arg;
    if (looks_like_pdb_id(arg)) {
        auto p = fetch_structure(arg, cdir, allow_download);
        if (!p.empty()) return p;
        throw std::runtime_error("could not obtain PDB id '" + arg +
                                 "' (download failed or --no-download set)");
    }
    throw std::runtime_error("input not found: " + arg +
                             " (give a structure file path or a PDB id)");
}

// Hidden debugging subcommands keyed by name. Each takes a single
// <input.pdb> path; classify-suite (7 angles) is handled apart.
using DumpFn = int (*)(const std::string&);
const std::pair<const char*, DumpFn> kDumpCommands[] = {
    {"dump-parse", pfcli::dump_parse},
    {"dump-frames", pfcli::dump_frames},
    {"dump-hbonds", pfcli::dump_hbonds},
    {"dump-candidates", pfcli::dump_candidates},
    {"dump-lwclass", pfcli::dump_lwclass},
    {"dump-classified", pfcli::dump_classified},
    {"dump-scored", pfcli::dump_scored},
    {"dump-pairs", pfcli::dump_pairs},
    {"dump-selection", pfcli::dump_selection},
    {"dump-scores", pfcli::dump_scores},
    {"dump-structure", pfcli::dump_structure},
    {"dump-torsions", pfcli::dump_torsions},
    {"dump-pucker", pfcli::dump_pucker},
};

}  // namespace

int main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty() || args[0] == "-h" || args[0] == "--help") {
        std::cout << "rubric " << kVersion << "\n"
                  << "usage:\n"
                  << "  rubric <pdb-id | structure-file> [options]\n"
                  << "      find + classify + score -> JSON (stdout, or a file with --out)\n"
                  << "      a bare PDB id (e.g. 6V3A) is downloaded from RCSB and cached;\n"
                  << "      an existing .cif/.pdb(.gz) path is used directly.\n"
                  << "  --out [FILE]  write JSON to a file (default: stdout). A bare --out\n"
                  << "                writes <id>.json in the current directory; any name or\n"
                  << "                path you give is used exactly as-is\n"
                  << "  --pml [FILE]  also write a PyMOL script showing the flagged pairs and\n"
                  << "                pucker outliers (bare --pml -> <id>.pml in the current\n"
                  << "                directory); open it with `pymol FILE`\n"
                  << "  options: --details  --no-score  --no-download  --cache-dir DIR\n"
                  << "  cache:   $RUBRIC_CACHE_DIR or ~/.cache/rubric\n"
                  << "  rubric dump-parse <input.pdb>   # parsed atoms (TSV)\n";
        return args.empty() ? 1 : 0;
    }
    if (args[0] == "--version") {
        std::cout << kVersion << "\n";
        return 0;
    }

    // Hidden dump-* commands: <command> <input.pdb>.
    for (const auto& [name, fn] : kDumpCommands) {
        if (args[0] == name) {
            if (args.size() < 2) {
                std::cerr << args[0] << ": missing <input.pdb>\n";
                return 1;
            }
            return fn(args[1]);
        }
    }
    if (args[0] == "classify-suite") {
        if (args.size() < 8) { std::cerr << "classify-suite: need 7 angles\n"; return 1; }
        return pfcli::classify_suite(std::vector<std::string>(args.begin() + 1, args.begin() + 8));
    }
    // Regenerate the precomputed resource caches (H-bond patterns from the
    // idealized templates; base-typing + H-bond-chemistry tables from the ligand
    // DB) so the runtime use_cache paths avoid re-parsing them every startup.
    // `rubric dump-caches [resources_root]`: with a root, writes into that source
    // tree; otherwise into the resolved (staged) resource locations.
    if (args[0] == "dump-caches" || args[0] == "dump-hbond-patterns") {
        using namespace pairfinder;
        const auto cfg = pfcli::cli_pipeline_config();
        const std::filesystem::path cfg_dir = std::filesystem::path(cfg.ligand_db).parent_path();
        std::filesystem::path patterns_out = std::filesystem::path(cfg.idealized_dir) /
                                             "hbond_patterns_cache.json";
        std::filesystem::path typing_out = cfg_dir / "base_typing_cache.json";
        std::filesystem::path chem_out = cfg_dir / "hbond_chemistry_cache.json";
        std::filesystem::path tables_out =
            std::filesystem::path(cfg.reference_dir) / "scoring_tables_cache.bin";
        if (args.size() >= 2) {  // explicit resources root (source tree)
            const std::filesystem::path root = args[1];
            patterns_out = root / "basepair-idealized" / "hbond_patterns_cache.json";
            typing_out = root / "config" / "base_typing_cache.json";
            chem_out = root / "config" / "hbond_chemistry_cache.json";
            tables_out = root / "reference" / "scoring_tables_cache.bin";
        }
        algorithms::hbond::HBondChemistry chem(cfg.ligand_db, /*use_cache=*/false);
        algorithms::classification::HBondPatterns patterns(cfg.idealized_dir, chem,
                                                           /*use_cache=*/false);
        core::BaseTyping typing(cfg.ligand_db, /*use_cache=*/false);
        patterns.save_cache(patterns_out);
        typing.save_cache(typing_out);
        chem.save_cache(chem_out);
        // Scoring tables: hashed against the reference JSON, so they stay valid
        // only while those files are unchanged (see Scorer::save_tables_cache).
        scoring::Scorer(cfg.reference_dir + "/parameter_distributions.json",
                        cfg.reference_dir + "/prosco_distributions.json",
                        cfg.reference_dir + "/penalty_weights.json", cfg.richardson_suites)
            .save_tables_cache(tables_out);
        std::cerr << "wrote " << patterns_out << ", " << typing_out << ", " << chem_out << ", "
                  << tables_out << "\n";
        return 0;
    }

    // Default: `rubric <pdb-id | structure-file> [options]` runs the full
    // pipeline and emits the report JSON. A bare PDB id is fetched
    // from RCSB (cached) unless an existing file path is given.
    const std::string input = args[0];
    bool score = true, details = false, allow_download = true, emit_classified = false;
    bool out_requested = false, pml_requested = false;
    std::string out_path, pml_path, cache_override;
    for (std::size_t i = 1; i < args.size(); ++i) {
        if (args[i] == "--no-score") score = false;
        else if (args[i] == "--score") score = true;
        else if (args[i] == "--details") details = true;
        else if (args[i] == "--classified") emit_classified = true;
        else if (args[i] == "--no-download") allow_download = false;
        else if (args[i] == "--out") {
            // Optional value: `--out name.json` uses that name; a bare `--out`
            // (nothing, or followed by another flag) defaults to <pdb_id>.json.
            out_requested = true;
            if (i + 1 < args.size() && !args[i + 1].empty() && args[i + 1][0] != '-')
                out_path = args[++i];
        }
        else if (args[i] == "--pml") {
            // Same optional-value rule as --out; a bare --pml -> <pdb_id>.pml.
            pml_requested = true;
            if (i + 1 < args.size() && !args[i + 1].empty() && args[i + 1][0] != '-')
                pml_path = args[++i];
        }
        else if (args[i] == "--cache-dir" && i + 1 < args.size()) cache_override = args[++i];
        else { std::cerr << "unknown option: " << args[i] << "\n"; return 1; }
    }
    if (pml_requested && !score) {
        std::cerr << "--pml needs scoring (it colors pairs by score); drop --no-score\n";
        return 1;
    }
    try {
        const std::filesystem::path file =
            resolve_input(input, cache_dir(cache_override), allow_download);
        // `--out` with no name defaults to <pdb_id>.json in the current directory.
        if (out_requested && out_path.empty())
            out_path = clean_pdb_id(file.string()) + ".json";
        if (pml_requested && pml_path.empty())
            pml_path = clean_pdb_id(file.string()) + ".pml";
        return find_json(file.string(), score, details, emit_classified, out_path, pml_path);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
