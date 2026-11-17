/**
 * Awkward but legal input files. Each variant is made at run time from the 1GID
 * fragment (chain A 195-202) and must give the same backbone and pucker results
 * as the original, or a clean empty result when there is no RNA.
 */
#include <algorithm>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

#include "engine_runner.hpp"
#include "test_helper.hpp"

using rubric_test::run_engine;

namespace {

const auto kFragment = rubric_test::tests_dir() / "fragments" / "1GID_A195-202.pdb";

bool is_atom_line(const std::string& l) { return l.rfind("ATOM", 0) == 0 || l.rfind("HETATM", 0) == 0; }

/// Apply `edit` to every atom line of a PDB text.
std::string map_atoms(const std::string& pdb, const std::function<std::string(const std::string&)>& edit) {
    std::istringstream in(pdb);
    std::ostringstream out;
    std::string line;
    while (std::getline(in, line)) out << (is_atom_line(line) ? edit(line) : line) << "\n";
    return out.str();
}

/// Set residue number (columns 23-26) and insertion code (column 27) of an atom line.
std::string set_number(const std::string& line, int num, char icode = ' ') {
    char field[6];
    std::snprintf(field, sizeof field, "%4d%c", num, icode);
    return line.substr(0, 22) + field + line.substr(27);
}

int number_of(const std::string& line) { return std::stoi(line.substr(22, 4)); }

/// Sorted suiteness values and pucker-outlier count: the backbone result of a structure,
/// independent of residue numbering.
struct BackboneSummary {
    std::vector<double> suiteness;
    int pucker_outliers = 0;
    bool operator==(const BackboneSummary& o) const {
        return suiteness == o.suiteness && pucker_outliers == o.pucker_outliers;
    }
};

std::ostream& operator<<(std::ostream& os, const BackboneSummary& s) {
    os << s.suiteness.size() << " suites, " << s.pucker_outliers << " pucker outliers";
    return os;
}

BackboneSummary summary(const std::string& name, const std::string& pdb) {
    const auto r = run_engine(rubric_test::write_temp(name, pdb));
    BackboneSummary s;
    for (const auto& t : r.score.backbone_torsions)
        if (t.scored) s.suiteness.push_back(t.suiteness);
    std::sort(s.suiteness.begin(), s.suiteness.end());
    for (const auto& rs : r.score.residue_scores) s.pucker_outliers += rs.pucker_outlier ? 1 : 0;
    return s;
}

const std::string& original() {
    static const std::string text = rubric_test::read_text(kFragment);
    return text;
}

const BackboneSummary& reference() {
    static const BackboneSummary s = summary("original.pdb", original());
    return s;
}

}  // namespace

/** The unmodified fragment: seven suites and one pucker outlier. */
TEST(reference_fragment) {
    CHECK_EQ(reference().suiteness.size(), std::size_t{7});
    CHECK_EQ(reference().pucker_outliers, 1);
}

/** Negative residue numbers (-3 to 4) give the same result. */
TEST(negative_numbers) {
    const auto pdb = map_atoms(original(), [](const std::string& l) { return set_number(l, number_of(l) - 198); });
    CHECK_EQ(summary("negative.pdb", pdb), reference());
}

/** Insertion codes (199, 199A, 199B, 199C for 199-202) give the same result. */
TEST(insertion_codes) {
    const auto pdb = map_atoms(original(), [](const std::string& l) {
        const int n = number_of(l);
        return n <= 199 ? l : set_number(l, 199, static_cast<char>('A' + (n - 200)));
    });
    CHECK_EQ(summary("icodes.pdb", pdb), reference());
}

/** In a file with two models, only the first is analysed. */
TEST(first_model_only) {
    std::ostringstream two;
    two << "MODEL        1\n" << original() << "ENDMDL\nMODEL        2\n"
        << rubric_test::shift_chain(original(), 'A', 40.0, 0.0, 0.0) << "ENDMDL\n";
    CHECK_EQ(summary("two_models.pdb", two.str()), reference());
}

/** An alternate location B of an atom is ignored in favour of location A. */
TEST(alternate_locations) {
    std::ostringstream out;
    std::istringstream in(original());
    std::string line;
    while (std::getline(in, line)) {
        if (is_atom_line(line) && line.substr(12, 4) == " C4'") {
            std::string a = line, b = rubric_test::shift_chain(line + "\n", 'A', 1.5, 0.0, 0.0);
            a[16] = 'A';
            b[16] = 'B';
            out << a << "\n" << b.substr(0, b.size() - 1) << "\n";
        } else {
            out << line << "\n";
        }
    }
    CHECK_EQ(summary("altloc.pdb", out.str()), reference());
}

/** Protein, water and ions next to the RNA do not change the RNA result. */
TEST(non_rna_residues_ignored) {
    const std::string extra =
        "ATOM   9001  N   ALA B   1      30.000  40.000  50.000  1.00 20.00           N\n"
        "ATOM   9002  CA  ALA B   1      31.400  40.000  50.000  1.00 20.00           C\n"
        "ATOM   9003  C   ALA B   1      32.000  41.400  50.000  1.00 20.00           C\n"
        "ATOM   9004  O   ALA B   1      31.300  42.400  50.000  1.00 20.00           O\n"
        "HETATM 9005  O   HOH W   1      25.000  45.000  55.000  1.00 30.00           O\n"
        "HETATM 9006 MG    MG M   1      26.000  46.000  56.000  1.00 30.00          MG\n";
    CHECK_EQ(summary("with_protein.pdb", original() + extra), reference());
}

/** A file with no RNA, or no atoms at all, gives an empty result. */
TEST(no_rna) {
    const std::string protein =
        "ATOM      1  N   ALA A   1      30.000  40.000  50.000  1.00 20.00           N\n"
        "ATOM      2  CA  ALA A   1      31.400  40.000  50.000  1.00 20.00           C\n"
        "END\n";
    for (const auto& [name, text] : {std::pair<std::string, std::string>{"protein.pdb", protein},
                                     std::pair<std::string, std::string>{"empty.pdb", ""}}) {
        const auto r = run_engine(rubric_test::write_temp(name, text));
        CHECK_MSG(r.selected.empty(), name);
        CHECK_MSG(r.score.pair_scores.empty(), name);
        CHECK_MSG(r.score.residue_scores.empty(), name);
    }
}

RUBRIC_TEST_MAIN()
