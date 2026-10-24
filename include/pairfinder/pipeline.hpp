/**
 * @file pipeline.hpp
 * @brief End-to-end find + classify + select orchestration as a reusable object.
 *
 * Bundles the resource objects (templates, typing, H-bond chemistry, template
 * aligner, H-bond patterns, quality scorer) that the pipeline needs - loaded
 * once at construction - and exposes the stages so callers don't have to re-wire
 * the stages by hand:
 *
 *   PreparedStructure prep = pipeline.prepare(path);   // parse + frames + view
 *   auto scored   = pipeline.classify(prep);           // find + classify
 *   auto selected = pipeline.select(scored, prep.structure);  // selection
 *   auto scorer   = pipeline.make_scorer();            // empirical scorer
 *
 * The individual resource accessors are exposed for debug paths that need a
 * different stage ordering (e.g. per-candidate classification dumps).
 */
#ifndef PAIRFINDER_PIPELINE_HPP
#define PAIRFINDER_PIPELINE_HPP

#include <string>
#include <utility>
#include <vector>

#include <pairfinder/algorithms/base_frame_calculator.hpp>
#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond_patterns.hpp>
#include <pairfinder/algorithms/quality_scorer.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>
#include <pairfinder/algorithms/template_aligner.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/reference_frame.hpp>
#include <pairfinder/core/structure.hpp>
#include <pairfinder/scoring/scorer.hpp>

namespace pairfinder {

/// All resource paths the pipeline reads. Each defaults to the bundled
/// ``resources/`` tree via PipelineConfig::defaults(); the CLI layers its
/// environment-variable overrides on top before constructing the Pipeline.
struct PipelineConfig {
    std::string templates_dir;        ///< Atomic_<X>.pdb base templates
    std::string nucleotides_config;   ///< modified_nucleotides.json
    std::string ligand_db;            ///< ligand_hbond_db.json
    std::string idealized_dir;        ///< basepair-idealized templates
    std::string exemplar_dir;         ///< basepair-catalog-exemplars
    std::string reference_dir;        ///< empirical-scorer reference data dir
    std::string richardson_suites;    ///< richardson_suites.json

    /// Paths resolved from the bundled resources/ tree (no env-var overrides).
    static PipelineConfig defaults();
};

/// A parsed structure with its per-residue reference frames. The structure has
/// the pipeline base-type view applied (DC->C, DT->T), exactly as the classifier
/// and downstream stages expect; the frames were computed from the original
/// (pre-view) residue names.
struct PreparedStructure {
    core::Structure structure;
    std::vector<std::pair<std::string, core::ReferenceFrame>> frames;
};

/// Owns the loaded pipeline resources and runs the find/classify/select stages.
/// Construct once per process (loading is the expensive part), then reuse.
class Pipeline {
public:
    explicit Pipeline(PipelineConfig config = PipelineConfig::defaults());

    /// parse_pdb + calculate_all_frames + apply_pipeline_view (per residue).
    PreparedStructure prepare(const std::string& path) const;

    /// Full per-candidate scoring pass (find H-bonds + classify + quality).
    std::vector<algorithms::classification::ScoredCandidate> classify(
        const PreparedStructure& prep);

    /// Build RNA chains and run two-phase selection over the scored candidates.
    std::vector<algorithms::classification::ScoredCandidate> select(
        std::vector<algorithms::classification::ScoredCandidate> scored,
        const core::Structure& structure) const;

    /// Construct the empirical scorer from the configured reference data.
    scoring::Scorer make_scorer() const;

    // --- resource accessors (for debug paths with a custom stage order) -------
    const PipelineConfig& config() const { return config_; }
    const algorithms::BaseTemplateLibrary& templates() const { return lib_; }
    const core::BaseTyping& typing() const { return typing_; }
    algorithms::hbond::HBondChemistry& chem() { return chem_; }
    const algorithms::classification::TemplateAligner& aligner() const { return aligner_; }
    const algorithms::classification::HBondPatterns& patterns() const { return patterns_; }
    const algorithms::classification::QualityScorer& quality_scorer() const { return qscorer_; }

private:
    PipelineConfig config_;
    // Declaration order is load order; lib_ stores a pointer into registry_ and
    // patterns_ is built from chem_, so registry_/chem_ must precede them.
    core::NucleotideRegistry registry_;
    algorithms::BaseTemplateLibrary lib_;
    core::BaseTyping typing_;
    algorithms::hbond::HBondChemistry chem_;
    algorithms::classification::TemplateAligner aligner_;
    algorithms::classification::HBondPatterns patterns_;
    algorithms::classification::QualityScorer qscorer_;
};

}  // namespace pairfinder

#endif  // PAIRFINDER_PIPELINE_HPP
