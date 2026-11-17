/**
 * Nucleotide chemistry and the bundled resources: residue typing (parent base,
 * purine/pyrimidine, DNA, glycosidic N), the modified-nucleotide table, which
 * pairs are scored (DNA and ambiguous pairs are not), and resource loading.
 */
#include <filesystem>
#include <sstream>
#include <string>

#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/pipeline.hpp>

#include "engine_runner.hpp"
#include "test_helper.hpp"

using rubric_test::run_engine;
namespace fs = std::filesystem;

namespace {

/// Copy of a PDB text with every residue of `chain` renamed to `name` (columns 18-20).
std::string rename_chain(const std::string& pdb, char chain, const std::string& name) {
    std::istringstream in(pdb);
    std::ostringstream out;
    std::string line;
    char field[4];
    std::snprintf(field, sizeof field, "%3s", name.c_str());
    while (std::getline(in, line)) {
        const bool atom = line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0;
        if (atom && line.size() > 21 && line[21] == chain) {
            line = "HETATM" + line.substr(6);
            line = line.substr(0, 17) + field + line.substr(20);
        }
        out << line << "\n";
    }
    return out.str();
}

fs::path gc_template() { return pairfinder::resources::basepair_idealized() / "cWW" / "GC.pdb"; }

}  // namespace

/** Standard and modified residue names map to their parent base. */
TEST(parent_base_of_residue_names) {
    const auto& t = rubric_test::pipeline().typing();
    CHECK_EQ(t.normalize("G"), std::string("G"));
    CHECK_EQ(t.normalize("U"), std::string("U"));
    CHECK_EQ(t.normalize("PSU"), std::string("U"));  // pseudouridine
    CHECK_EQ(t.normalize("OMG"), std::string("G"));  // 2'-O-methylguanosine
    CHECK_EQ(t.normalize("5MC"), std::string("C"));  // 5-methylcytidine
    CHECK_EQ(t.normalize("1MA"), std::string("A"));  // 1-methyladenosine
}

/** Purines, pyrimidines, DNA, and the glycosidic nitrogen of each base. */
TEST(base_properties) {
    const auto& t = rubric_test::pipeline().typing();
    CHECK(t.is_purine("A") && t.is_purine("G"));
    CHECK(t.is_pyrimidine("C") && t.is_pyrimidine("U"));
    CHECK(!t.is_purine("U") && !t.is_pyrimidine("G"));
    CHECK(t.is_dna("DA") && t.is_dna("DG") && t.is_dna("DT"));
    CHECK(!t.is_dna("A") && !t.is_dna("G") && !t.is_dna("PSU"));
    CHECK(t.is_known_nucleotide("PSU"));
    CHECK(!t.is_known_nucleotide("HOH"));
    CHECK_EQ(t.glycosidic_n_name("G"), std::string("N9"));
    CHECK_EQ(t.glycosidic_n_name("A"), std::string("N9"));
    CHECK_EQ(t.glycosidic_n_name("C"), std::string("N1"));
    CHECK_EQ(t.glycosidic_n_name("U"), std::string("N1"));
}

/** The modified-nucleotide table knows the standard and common modified residues. */
TEST(nucleotide_registry) {
    pairfinder::core::NucleotideRegistry reg(pairfinder::PipelineConfig::defaults().nucleotides_config);
    CHECK(reg.size() > 100);
    const auto omg = reg.lookup("OMG");  // 2'-O-methylguanosine, listed as a modified G
    CHECK(omg.has_value());
    CHECK(omg->is_modified && omg->is_purine);
    CHECK_EQ(omg->canonical, 'G');
    CHECK(reg.lookup("PSU").has_value());  // pseudouridine has its own code (P) in the table
    const auto g = reg.lookup("G");
    CHECK(g.has_value() && g->is_purine && !g->is_modified);
    CHECK(!reg.lookup("ZZZ").has_value());
}

/** A pair of 2'-O-methyl nucleotides is found, typed by its parents and scored. */
TEST(modified_pair_is_scored) {
    std::string pdb = rubric_test::read_text(gc_template());
    pdb = rename_chain(rename_chain(pdb, 'A', "OMG"), 'B', "OMC");
    const auto r = run_engine(rubric_test::write_temp("omg_omc.pdb", pdb), false);
    CHECK_EQ(r.selected.size(), std::size_t{1});
    CHECK_EQ(r.selected[0].lw_class, std::string("cWW"));
    CHECK_EQ(r.score.pair_scores.size(), std::size_t{1});
    CHECK_EQ(r.score.pair_scores[0].bp_type, std::string("G-C"));
}

/** A pair with a DNA nucleotide on either side is found and classified but not scored. */
TEST(dna_pairs_are_not_scored) {
    const std::string pdb = rubric_test::read_text(gc_template());
    for (const auto& variant : {rename_chain(pdb, 'A', "DG"),
                                rename_chain(rename_chain(pdb, 'A', "DG"), 'B', "DC")}) {
        const auto r = run_engine(rubric_test::write_temp("dna_pair.pdb", variant), false);
        CHECK_EQ(r.selected.size(), std::size_t{1});
        CHECK(r.score.pair_scores.empty());
        CHECK_EQ(r.score.skipped_pairs, 1);
    }
}

/** Ambiguous pairs are reported but never scored; every other base-base pair is. */
TEST(ambiguous_pairs_are_not_scored) {
    const auto r = run_engine(rubric_test::tests_dir() / "data" / "1GID.cif.gz", false);
    int ambiguous = 0;
    for (const auto& c : r.selected) {
        if (!c.is_ambiguous) continue;
        ++ambiguous;
        for (const auto& p : r.score.pair_scores)
            CHECK_MSG(!(p.res_id1 == c.res_id1 && p.res_id2 == c.res_id2), c.res_id1 << " " << c.res_id2);
    }
    CHECK(ambiguous > 0);
    CHECK_EQ(r.score.pair_scores.size() + static_cast<std::size_t>(r.score.skipped_pairs), r.selected.size());
}

/** The bundled resources are found and every configured path exists. */
TEST(resources_are_found) {
    CHECK(!pairfinder::resources::root().empty());
    CHECK(fs::is_directory(pairfinder::resources::templates_dir()));
    CHECK(fs::is_directory(pairfinder::resources::basepair_idealized()));
    CHECK(fs::is_regular_file(pairfinder::resources::file("richardson_suites.json")));
    const auto cfg = pairfinder::PipelineConfig::defaults();
    for (const auto& p : {cfg.templates_dir, cfg.idealized_dir, cfg.exemplar_dir, cfg.reference_dir})
        CHECK_MSG(fs::is_directory(p), p);
    for (const auto& p : {cfg.nucleotides_config, cfg.ligand_db, cfg.richardson_suites})
        CHECK_MSG(fs::is_regular_file(p), p);
}

/** A missing structure file is reported as an error, not a crash or an empty result. */
TEST(missing_file_throws) {
    CHECK_THROWS(rubric_test::pipeline().prepare("/nonexistent/path/missing.cif"));
}

RUBRIC_TEST_MAIN()
