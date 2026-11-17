/**
 * @file engine_runner.hpp
 * @brief Run the RUBRIC engine on a structure file inside a test.
 *
 * The same stages as the `rubric` command (apps/main.cpp): prepare (parse +
 * frames), classify, select, then score. The Pipeline and Scorer load the
 * bundled resources once per test program and are reused.
 */
#ifndef RUBRIC_TESTS_ENGINE_RUNNER_HPP
#define RUBRIC_TESTS_ENGINE_RUNNER_HPP

#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/pipeline.hpp>
#include <pairfinder/scoring/scorer.hpp>

namespace rubric_test {

namespace classification = pairfinder::algorithms::classification;

/// Shared pipeline and scorer (loading the resources is the slow part).
inline pairfinder::Pipeline& pipeline() {
    static pairfinder::Pipeline instance;
    return instance;
}

inline pairfinder::scoring::Scorer& scorer() {
    static pairfinder::scoring::Scorer instance = pipeline().make_scorer();
    return instance;
}

/// Everything one engine run produces.
struct EngineResult {
    std::unique_ptr<pairfinder::PreparedStructure> prep;  ///< owns the frames the pairs point to
    std::vector<classification::ScoredCandidate> selected;
    pairfinder::scoring::StructureScore score;
};

/// Run the engine on a structure file. collect_backbone fills the full
/// per-suite torsion table (as `rubric --details` does).
inline EngineResult run_engine(const std::filesystem::path& path, bool collect_backbone = true) {
    EngineResult r;
    r.prep = std::make_unique<pairfinder::PreparedStructure>(pipeline().prepare(path.string()));
    auto scored = pipeline().classify(*r.prep);
    r.selected = pipeline().select(std::move(scored), r.prep->structure);
    pairfinder::algorithms::hbond::HBondFinder finder(pipeline().chem());
    r.score = scorer().score_structure(r.selected, r.prep->structure, finder, pipeline().chem(),
                                       pipeline().typing(), collect_backbone);
    return r;
}

/// Write text to a file in the system temporary directory and return its path.
inline std::filesystem::path write_temp(const std::string& name, const std::string& text) {
    const auto dir = std::filesystem::temp_directory_path() / "rubric_tests";
    std::filesystem::create_directories(dir);
    const auto path = dir / name;
    std::ofstream(path, std::ios::binary) << text;
    return path;
}

/// Read a whole text file.
inline std::string read_text(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream os;
    os << in.rdbuf();
    return os.str();
}

/// Copy a PDB-format text, moving every atom of one chain by (dx, dy, dz) Å.
/// Columns follow the PDB format (x 31-38, y 39-46, z 47-54, chain 22).
inline std::string shift_chain(const std::string& pdb, char chain, double dx, double dy, double dz) {
    std::istringstream in(pdb);
    std::ostringstream out;
    std::string line;
    while (std::getline(in, line)) {
        const bool atom = line.rfind("ATOM", 0) == 0 || line.rfind("HETATM", 0) == 0;
        if (atom && line.size() >= 54 && line[21] == chain) {
            const double x = std::stod(line.substr(30, 8)) + dx;
            const double y = std::stod(line.substr(38, 8)) + dy;
            const double z = std::stod(line.substr(46, 8)) + dz;
            char coords[32];
            std::snprintf(coords, sizeof coords, "%8.3f%8.3f%8.3f", x, y, z);
            line = line.substr(0, 30) + coords + line.substr(54);
        }
        out << line << "\n";
    }
    return out.str();
}

}  // namespace rubric_test

#endif  // RUBRIC_TESTS_ENGINE_RUNNER_HPP
