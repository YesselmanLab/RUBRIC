/**
 * @file richardson.hpp
 * @brief Richardson RNA backbone suite classifier.
 *
 * Loads richardson_suites.json (54 conformers + widths + triage/sieve +
 * satellite overrides) and classifies a 7D suite
 * [delta_prev, eps_prev, zeta_prev, alpha, beta, gamma, delta] into a suiteness
 * score in [0,1] (0 for triage/sieve/4D/7D outliers).
 */
#ifndef PAIRFINDER_SCORING_RICHARDSON_HPP
#define PAIRFINDER_SCORING_RICHARDSON_HPP

#include <array>
#include <filesystem>
#include <memory>
#include <string>

namespace pairfinder::scoring {

/// Names of the seven suite dihedrals, in suite-array order:
/// [delta(i-1), epsilon(i-1), zeta(i-1), alpha(i), beta(i), gamma(i), delta(i)].
inline constexpr std::array<const char*, 7> kSuiteAngleNames = {
    "delta_prev", "epsilon_prev", "zeta_prev",
    "alpha", "beta", "gamma", "delta"};

/// Result of classifying one 7D backbone suite.
struct SuiteResult {
    std::string conformer;   ///< matched conformer name, or "!!" for outliers
    double suiteness = 0.0;  ///< in [0,1], 0 for outliers
    bool is_outlier = true;
    double distance = -1.0;  ///< 7D L3 distance to the matched center (-1 if none)
    std::string bin;         ///< pucker-pucker-gamma bin ("33p"); empty if triage/sieve fail
};

/// Signed minimal angular gap value-center in (-180, 180]. Positive => rotate
/// the angle down toward center; negative => rotate up.
double signed_angle_gap(double value, double center);

class RichardsonClassifier {
public:
    explicit RichardsonClassifier(const std::filesystem::path& suites_path);
    ~RichardsonClassifier();
    RichardsonClassifier(RichardsonClassifier&&) noexcept;

    /// Suiteness in [0,1] (0 for outliers). suite = 7 angles, each already %360.
    /// Lean hot path (no allocation) - used to compute the residue score.
    double suiteness(const std::array<double, 7>& suite) const;

    /// Full classification: conformer name + suiteness + outlier flag. Used when
    /// the conformer identity is needed (backbone recommendation).
    SuiteResult classify(const std::array<double, 7>& suite) const;

    /// 7D cluster center for a conformer name, or nullptr if unknown.
    const std::array<double, 7>* center_for(const std::string& conformer) const;

    /// Nearest conformer by summed absolute angular gap (L1) - fix target for an
    /// outlier suite that fits no cluster. Empty string if the table is empty.
    std::string nearest_conformer(const std::array<double, 7>& suite) const;

    /// Baseline per-axis tolerance widths (fallback flagging for sparse cells).
    const std::array<double, 7>& normal_widths() const;

    /// Sugar pucker from a residue's delta torsion: 3 = C3'-endo, 2 = C2'-endo,
    /// 0 = between the two sieve bands (C3' [60,105], C2' [125,165]).
    int pucker(double delta) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace pairfinder::scoring

#endif  // PAIRFINDER_SCORING_RICHARDSON_HPP
