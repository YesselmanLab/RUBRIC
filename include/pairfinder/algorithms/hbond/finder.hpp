/**
 * @file finder.hpp
 * @brief Slot-based H-bond finder with global (per-pair) optimal selection.
 *
 * Slot-based HBondFinder. The pipeline
 * entry point is find_between(res1, res2): collect donor->acceptor candidates
 * within max_distance, score each against predicted donor-H / acceptor-LP slot
 * directions, then pick the maximum non-conflicting set (two candidates conflict
 * if they share a donor- or acceptor-slot). Backtracking MIS for <=50
 * candidates, greedy for larger sets.
 */
#ifndef PAIRFINDER_ALGORITHMS_HBOND_FINDER_HPP
#define PAIRFINDER_ALGORITHMS_HBOND_FINDER_HPP

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <pairfinder/algorithms/hbond/chemistry.hpp>
#include <pairfinder/algorithms/hbond/hbond.hpp>
#include <pairfinder/core/residue.hpp>
#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::algorithms::hbond {

/// A potential H-bond before selection (carries positions + chosen slot idx).
struct HBondCandidate {
    std::string donor_res_id;
    std::string acceptor_res_id;
    std::string donor_atom;
    std::string acceptor_atom;
    double distance = 0.0;
    Vector3d donor_pos;
    Vector3d acceptor_pos;
    int h_slot_idx = -1;
    int lp_slot_idx = -1;
    double alignment_score = 0.0;
};

/// Cached atom coordinates + predicted slots for one residue.
class ResidueSlotCache {
public:
    ResidueSlotCache(std::string res_id, std::string base_type, AtomMap atoms,
                     const HBondChemistry& chem);

    const std::vector<HbSlot>& get_h_slots(const std::string& atom_name);
    const std::vector<HbSlot>& get_lp_slots(const std::string& atom_name);

private:
    std::string res_id_;
    std::string base_type_;
    AtomMap atoms_;
    const HBondChemistry* chem_;
    std::optional<Vector3d> base_normal_;
    std::unordered_map<std::string, std::vector<HbSlot>> h_slots_;
    std::unordered_map<std::string, std::vector<HbSlot>> lp_slots_;
};

/// Finds H-bonds between residues using slot-based optimal selection.
class HBondFinder {
public:
    HBondFinder(const HBondChemistry& chemistry, double max_distance = 4.0,
                double min_alignment = 1.4, double short_distance_threshold = 3.2);

    /// Selected H-bonds between two residues (the pipeline entry point).
    std::vector<HBond> find_between(const core::Residue& res1, const core::Residue& res2);

    /// Selected H-bonds PLUS dropped candidates within distance_threshold, for
    /// the LW classifier's face voting only.
    std::vector<HBond> find_face_voters(const core::Residue& res1, const core::Residue& res2,
                                        double distance_threshold = 4.0);

    /// All candidates (both directions) within max_distance (no selection).
    std::vector<HBondCandidate> find_candidates(const core::Residue& res1,
                                                const core::Residue& res2);

private:
    using ResidueMap = std::unordered_map<std::string, const core::Residue*>;

    ResidueSlotCache& get_cache(const core::Residue& residue);
    void add_directional_candidates(const core::Residue& donor,
                                    const core::Residue& acceptor,
                                    std::vector<HBondCandidate>& out);
    void compute_alignments(std::vector<HBondCandidate>& candidates,
                            const ResidueMap& residues);
    std::vector<HBond> select_optimal(std::vector<HBondCandidate>& candidates,
                                      const ResidueMap& residues);
    std::vector<int> find_best_selection(const std::vector<HBondCandidate>& valid,
                                         const ResidueMap& residues);

    const HBondChemistry* chem_;
    double max_distance_;
    double min_alignment_;
    double short_distance_threshold_;
    std::unordered_map<std::string, ResidueSlotCache> residue_cache_;
};

}  // namespace pairfinder::algorithms::hbond

#endif  // PAIRFINDER_ALGORITHMS_HBOND_FINDER_HPP
