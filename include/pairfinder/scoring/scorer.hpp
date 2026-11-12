/**
 * @file scorer.hpp
 * @brief Empirical 0-100 base-pair quality scorer:
 *        score = 100 - sum(pair_weight[issue] * severity) over the scored
 *        parameters, severities driven by ProSco density tables, the graded
 *        H-bond count and the steric clash overlap.
 *
 * Loads reference data (parameter_distributions.json + prosco_distributions.json
 * + penalty_weights.json), which the engine consumes but does not generate.
 */
#ifndef PAIRFINDER_SCORING_SCORER_HPP
#define PAIRFINDER_SCORING_SCORER_HPP

#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond/finder.hpp>
#include <pairfinder/algorithms/score_candidates.hpp>
#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/core/structure.hpp>
#include <pairfinder/scoring/clashes.hpp>

namespace pairfinder::scoring {

namespace classification = pairfinder::algorithms::classification;
namespace hbond = pairfinder::algorithms::hbond;

/// One triggered issue and its weighted penalty contribution (weight x severity).
/// The Cerny three-tier provenance (ProSco / Z' / tier) is carried for reporting
/// and downstream highlighting (e.g. PyMOL); it does not affect the score.
struct IssuePenalty {
    std::string issue;
    double weight = 0.0;            ///< pair_weight[issue] x severity
    double severity = 0.0;          ///< Cerny severity in [0,1]
    std::optional<double> prosco;   ///< ProSco of the worst contributing value (none if N/A)
    std::optional<double> zprime;   ///< Z' of the worst contributing value (none if N/A)
    std::string tier;               ///< "Preferred"/"Allowed"/"Of Concern"; "" if not ProSco-scored
};

/// Pair-score triage tiers - global, data-derived cutoffs.
/// Secondary review spans [80, 100); below is Priority review.
constexpr double kTierAcceptableMin = 80.0;

/// Triage tier for a pair score:
///   "Preferred"  = 100     (every scored parameter in its Preferred/ProSco>=5 region)
///   "Secondary review" = [80,100) (minor strain)
///   "Priority review"  = <80      (worth a manual look)
/// Global rather than per-LW-class: all 12 LW medians are 100 and per-(class,bp)
/// ProSco/Z' normalization already removes the "non-canonical is harder" effect.
inline std::string pair_tier(double score) {
    if (score >= 100.0) return "Preferred";
    if (score >= kTierAcceptableMin) return "Secondary review";
    return "Priority review";
}

/// Report-only "unusual sugar pucker": a C2'-endo (by delta) sugar in a helical cWW
/// G-C/A-U pair, where A-form helices are C3'-endo. `pperp_confirms` = the Pperp
/// geometry also reads C2'-endo (genuine, atypical-but-real) vs C3'-endo (ring
/// likely mis-modeled -> really a pucker OUTLIER). Does NOT affect the score.
struct UnusualPucker {
    std::string res_id;
    std::optional<double> delta;
    std::optional<double> pperp;
    bool pperp_confirms = false;
};

/// Report-only sugar-pucker OUTLIER: the delta-derived pucker disagrees with the
/// Pperp-derived pucker (a mis-modeled ribose ring). Does NOT affect the score.
struct PuckerOutlier {
    double pperp = 0.0;
    std::optional<double> delta;
    std::string pucker_by_pperp;   ///< "C3'-endo" / "C2'-endo"
    std::string pucker_by_delta;
};

/// Report-only UNCLASSIFIED pucker: delta lies outside both Suitename sieve
/// bands (C3' [60,105], C2' [125,165]), so no pucker is assigned from delta.
/// Not an outlier and does NOT affect the score (the suite is already !!).
struct PuckerUnclassified {
    double delta = 0.0;
    std::optional<double> pperp;   ///< absent when there is no following P
    std::string pucker_by_pperp;   ///< "C3'-endo" / "C2'-endo"; empty without Pperp
};

/// Per-pair scoring result.
struct PairScore {
    std::string res_id1, res_id2;
    std::string bp_type, lw_class;
    double score = 0.0;
    double penalty = 0.0;
    std::vector<IssuePenalty> issues;  ///< triggered issues, sorted by descending penalty
    std::vector<UnusualPucker> unusual_pucker;  ///< report-only (helical C2'-endo)
    /// Worst scored heavy-atom overlap (A), signed; negative when separated.
    /// Raw and ungated - the feature behind the `clash` severity, so weight
    /// fitting can recompute it. nullopt if a residue has no heavy atoms.
    std::optional<double> clash_overlap;

    /// Triage tier derived from `score` (Preferred/Secondary review/Priority review).
    std::string tier() const { return pair_tier(score); }
};

/// One suite torsion flagged as genuinely off its target conformer. Reported
/// only when it *fires*: ProSco below the Preferred cutoff (5) against the
/// empirical per-conformer distribution. No Z'/penalty - suiteness is the score.
struct AngleDeviation {
    std::string angle;              ///< suite dihedral name (alpha, ..., delta_prev)
    double value = 0.0;             ///< measured angle, degrees
    double target = 0.0;            ///< recommended angle, degrees (damped)
    double gap = 0.0;               ///< signed measured - conformer centre, degrees
    std::optional<double> prosco;   ///< density percentile (lower = more anomalous)
    std::optional<double> conformer_center;  ///< the undamped cluster centre
};

/// Refinement-facing backbone recommendation for one residue (see
/// BackboneRecommendation). tier: "fixable" (>=1 firing torsion, correct these) /
/// "review" (low suiteness but no firing torsion - jointly distorted, rebuild) /
/// "flag_only" (outlier, nearest conformer is only a hint).
struct BackboneRecommendation {
    std::string tier;
    std::string target_conformer;
    std::vector<AngleDeviation> deviations;
};

/// Per-residue backbone score (suiteness * 100).
struct ResidueScore {
    std::string res_id;
    double suiteness = 0.0;
    double score = 0.0;
    std::optional<BackboneRecommendation> recommendation;  ///< none for ok residues
    std::optional<double> delta;    ///< delta torsion (ring-based pucker signal)
    std::optional<double> pperp;    ///< base-phosphate perpendicular distance (A)
    std::optional<PuckerOutlier> pucker_outlier;  ///< set only when delta/Pperp disagree
    std::optional<PuckerUnclassified> pucker_unclassified;  ///< set only when delta is outside both bands
};

/// Full per-residue backbone torsion row: every residue with a valid 7D suite, including outliers,
/// captured before the RNA-only / L-RNA scoring filters. Populated only when
/// score_structure is called with collect_backbone=true.
struct BackboneTorsion {
    std::string res_id;
    std::string base_type;  ///< normalized (PSU->U, ...)
    std::optional<double> alpha, beta, gamma, delta, epsilon, zeta;
    std::string conformer;  ///< "!!" for outliers
    double suiteness = 0.0;
    std::string suite_bin;  ///< empty -> null on emit
    bool is_outlier = true;
    bool scored = false;    ///< row enters the backbone score (RNA suite, not L-RNA)
};

/// Per-axis structure scores. `pairs_score` and `residues_score` are reported
/// separately: they are on different scales and measure different quality axes.
/// There is deliberately no combined "overall" - such a label overclaims
/// whole-structure quality (clashscore, bond geometry, etc. are not considered).
struct StructureScore {
    double pairs_score = 0.0;
    double residues_score = 0.0;
    int skipped_pairs = 0;  ///< selected pairs that were not scorable
    std::vector<PairScore> pair_scores;
    std::vector<ResidueScore> residue_scores;
    std::vector<BackboneTorsion> backbone_torsions;  ///< populated only when collect_backbone
    /// Every RNA-RNA steric clash in the structure, worst first - including
    /// clashes between residues that are not a base pair. Report-only at the
    /// structure level; the per-pair `clash` parameter is what affects scores.
    /// Protein contacts are deliberately excluded (see clashes.hpp).
    std::vector<ClashContact> clashes;

    /// Count of scored pairs per triage tier.
    struct TierSummary { int preferred = 0, acceptable = 0, review = 0; };
    TierSummary tier_summary() const {
        TierSummary t;
        for (const auto& ps : pair_scores) {
            const std::string tier = ps.tier();
            if (tier == "Preferred") ++t.preferred;
            else if (tier == "Secondary review") ++t.acceptable;
            else ++t.review;
        }
        return t;
    }
};

class Scorer {
public:
    Scorer(const std::filesystem::path& distributions_path,
           const std::filesystem::path& prosco_path,
           const std::filesystem::path& weights_path,
           const std::filesystem::path& richardson_suites_path);
    ~Scorer();
    Scorer(Scorer&&) noexcept;

    /// Scorable: base-base, not ambiguous, single (no '|') lw_class, and
    /// neither residue DNA (the scoring tables are calibrated on RNA).
    static bool is_scorable(const classification::ScoredCandidate& pair,
                            const core::BaseTyping& typing);

    /// Score one scorable pair. res1/res2 in the pair's residue order.
    PairScore score_pair(const classification::ScoredCandidate& pair,
                         const core::Residue& res1, const core::Residue& res2,
                         hbond::HBondFinder& finder, const hbond::HBondChemistry& chem,
                         const core::BaseTyping& typing) const;

    /// Whole-structure score: pair scores + residue suiteness.
    /// collect_backbone=true also fills StructureScore::backbone_torsions with a
    /// full per-residue torsion row (the --details backbone_torsions table).
    StructureScore score_structure(const std::vector<classification::ScoredCandidate>& selected,
                                   const core::Structure& structure, hbond::HBondFinder& finder,
                                   const hbond::HBondChemistry& chem,
                                   const core::BaseTyping& typing,
                                   bool collect_backbone = false) const;

    /// Write the scoring tables in native form (scoring_tables_cache.bin next to
    /// prosco_distributions.json unless @p path is given), stamped with the hash
    /// of the JSON they came from. Run by `rubric dump-caches` after each build.
    void save_tables_cache(const std::filesystem::path& path = {}) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_SCORER_HPP
