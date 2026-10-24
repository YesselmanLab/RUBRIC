/**
 * @file pipeline.cpp
 * @brief Pipeline orchestration (see pipeline.hpp).
 */
#include <pairfinder/pipeline.hpp>

#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/io/pdb_parser.hpp>

namespace pairfinder {

PipelineConfig PipelineConfig::defaults() {
    PipelineConfig c;
    c.templates_dir = resources::templates_dir().string();
    c.nucleotides_config = resources::config("modified_nucleotides.json").string();
    c.ligand_db = resources::config("ligand_hbond_db.json").string();
    c.idealized_dir = resources::basepair_idealized().string();
    c.exemplar_dir = resources::basepair_exemplars().string();
    c.reference_dir = resources::reference_dir().string();
    c.richardson_suites = resources::file("richardson_suites.json").string();
    return c;
}

Pipeline::Pipeline(PipelineConfig config)
    : config_(std::move(config)),
      registry_(config_.nucleotides_config),
      lib_(config_.templates_dir, registry_),
      typing_(config_.ligand_db),
      chem_(config_.ligand_db),
      aligner_(config_.idealized_dir, config_.exemplar_dir),
      patterns_(config_.idealized_dir, chem_),
      qscorer_() {}

PreparedStructure Pipeline::prepare(const std::string& path) const {
    PreparedStructure prep;
    prep.structure = io::parse_pdb(path);
    prep.frames = algorithms::calculate_all_frames(prep.structure, lib_);
    // Switch residues to pipeline base types (DC->C, DT->T) AFTER frames were
    // computed from the original names.
    for (auto& chain : prep.structure.chains())
        for (auto& res : chain.residues()) core::apply_pipeline_view(res);
    return prep;
}

std::vector<algorithms::classification::ScoredCandidate> Pipeline::classify(
    const PreparedStructure& prep) {
    return algorithms::classification::score_candidates(
        prep.structure, prep.frames, typing_, chem_, aligner_, patterns_, qscorer_);
}

std::vector<algorithms::classification::ScoredCandidate> Pipeline::select(
    std::vector<algorithms::classification::ScoredCandidate> scored,
    const core::Structure& structure) const {
    const auto chains = algorithms::RNAChains::from_structure(structure);
    return algorithms::selection::select_pairs(std::move(scored), chains);
}

scoring::Scorer Pipeline::make_scorer() const {
    return scoring::Scorer(config_.reference_dir + "/parameter_distributions.json",
                           config_.reference_dir + "/prosco_distributions.json",
                           config_.reference_dir + "/penalty_weights.json",
                           config_.richardson_suites);
}

}  // namespace pairfinder
