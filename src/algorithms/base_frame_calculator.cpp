/**
 * @file base_frame_calculator.cpp
 * @brief Per-base-template least-squares base frame.
 *
 * Each residue's observed ring atoms onto the ring
 * atoms are fitted onto those of that base's standard template (Atomic_<X>.pdb), whose coordinates
 * live in the standard reference frame. Purine (A,G,I) and pyrimidine (C,U,T)
 * rings - and even A vs G - have distinct standard geometries, so the template
 * is selected per base. The fit's (R, t, rms) is the residue's base frame.
 */
#include <pairfinder/algorithms/base_frame_calculator.hpp>

#include <algorithm>
#include <array>
#include <optional>

#include <pairfinder/core/atom_type.hpp>
#include <pairfinder/geometry/superposition.hpp>
#include <pairfinder/io/pdb_parser.hpp>

namespace pairfinder::algorithms {

namespace {

using core::AtomType;

// Ring atoms in 3DNA order. Pyrimidine templates carry only the first
// six; the remaining purine-only atoms are simply absent there.
constexpr std::array<AtomType, 9> kRingTypes = {
    AtomType::C4, AtomType::N3, AtomType::C2, AtomType::N1, AtomType::C6,
    AtomType::C5, AtomType::N7, AtomType::C8, AtomType::N9,
};

// Shared standard ring geometry (from 3DNA), used ONLY for the nucleotide
// type-check gate (not the per-base frame fit). Same order as kRingTypes.
constexpr std::array<std::array<double, 3>, 9> kStandardRing = {{
    {{-1.265, 3.177, 0.000}},  // C4
    {{-2.342, 2.364, 0.001}},  // N3
    {{-1.999, 1.087, 0.000}},  // C2
    {{-0.700, 0.641, 0.000}},  // N1
    {{0.424, 1.460, 0.000}},   // C6
    {{0.071, 2.833, 0.000}},   // C5
    {{0.870, 3.969, 0.000}},   // N7
    {{0.023, 4.962, 0.000}},   // C8
    {{-1.289, 4.551, 0.000}},  // N9
}};

// Ring-atom RMSD cutoff for accepting a residue as a nucleotide (3DNA's value).
constexpr double kNtRmsdCutoff = 0.2618;

// A frame is built only from primary-conformer atoms.
bool frame_altloc(const core::Atom* a) {
    if (a == nullptr) return false;
    return a->alt_loc == ' ' || a->alt_loc == 'A' || a->alt_loc == '1';
}

// Ring atom of a given type, but only if it is a primary-conformer atom.
const core::Atom* ring_atom(const core::Residue& r, AtomType t) {
    const core::Atom* a = r.get_atom_type(t);
    return frame_altloc(a) ? a : nullptr;
}

// Buffers and other non-nucleotide molecules, excluded outright (as in 3DNA).
bool is_excluded_molecule(const std::string& res_name) {
    static const std::array<const char*, 11> kExcluded = {
        "MES", "HEPES", "TRIS", "EDO", "GOL", "SO4",
        "PO4", "ACT",   "FMT",  "EFZ", "LYA"};
    return std::any_of(kExcluded.begin(), kExcluded.end(),
                       [&](const char* m) { return res_name == m; });
}

bool has_c1_prime(const core::Residue& r) {
    return frame_altloc(r.get_atom_type(AtomType::C1_PRIME)) ||
           frame_altloc(r.get_atom("C1R"));
}

// Type-check rms: fit present ring atoms onto the shared standard ring.
struct TypeCheck {
    std::optional<double> rms;
    bool found_purine = false;
};

TypeCheck nt_type_rms(const core::Residue& r) {
    std::vector<geometry::Vector3d> standard, experimental;
    int nN = 0, purine_count = 0;
    for (std::size_t i = 0; i < kRingTypes.size(); ++i) {
        const core::Atom* a = ring_atom(r, kRingTypes[i]);
        if (a == nullptr) continue;
        standard.push_back({kStandardRing[i][0], kStandardRing[i][1], kStandardRing[i][2]});
        experimental.push_back(a->coords);
        if (i == 1 || i == 3 || i == 6 || i == 8) ++nN;  // N3,N1,N7,N9
        if (i >= 6) ++purine_count;
    }
    TypeCheck tc;
    tc.found_purine = purine_count > 0;
    if ((nN == 0 && !has_c1_prime(r)) || experimental.size() < 3) return tc;
    const auto fit = geometry::superpose(standard, experimental);
    if (fit.valid) tc.rms = fit.rms;
    return tc;
}

// Pyrimidine fallback: fit only the first six ring atoms onto the shared ring.
std::optional<double> try_pyrimidine_rms(const core::Residue& r) {
    std::vector<geometry::Vector3d> standard, experimental;
    for (std::size_t i = 0; i < 6; ++i) {
        const core::Atom* a = ring_atom(r, kRingTypes[i]);
        if (a == nullptr) continue;
        standard.push_back({kStandardRing[i][0], kStandardRing[i][1], kStandardRing[i][2]});
        experimental.push_back(a->coords);
    }
    if (experimental.size() < 3) return std::nullopt;
    const auto fit = geometry::superpose(standard, experimental);
    if (!fit.valid) return std::nullopt;
    return fit.rms;
}

// The nucleotide gate: excluded molecules, then the ring-atom RMSD
// type-check with the purine-atom pyrimidine fallback.
bool passes_type_check(const core::Residue& r) {
    if (is_excluded_molecule(r.base_type())) return false;
    const TypeCheck tc = nt_type_rms(r);
    if (tc.rms.has_value() && *tc.rms <= kNtRmsdCutoff) return true;
    if (!tc.found_purine) return false;
    const auto pyr = try_pyrimidine_rms(r);
    return pyr.has_value() && *pyr <= kNtRmsdCutoff;
}

}  // namespace

const std::unordered_map<AtomType, geometry::Vector3d>*
BaseTemplateLibrary::ring_atoms(const std::string& res_name) const {
    const std::string fname = registry_->template_filename(res_name);
    if (fname.empty()) return nullptr;
    if (auto it = cache_.find(fname); it != cache_.end()) return &it->second;
    if (missing_.count(fname)) return nullptr;

    const auto path = templates_dir_ / fname;
    core::Structure tmpl;
    try {
        tmpl = io::parse_pdb(path.string());
    } catch (const std::exception&) {
        missing_[fname] = true;
        return nullptr;
    }

    std::unordered_map<AtomType, geometry::Vector3d> ring;
    for (const auto& chain : tmpl.chains()) {
        for (const auto& res : chain.residues()) {
            for (AtomType t : kRingTypes) {
                if (const core::Atom* a = res.get_atom_type(t)) ring[t] = a->coords;
            }
            if (!ring.empty()) break;
        }
        if (!ring.empty()) break;
    }
    if (ring.empty()) {
        missing_[fname] = true;
        return nullptr;
    }
    auto [it, _] = cache_.emplace(fname, std::move(ring));
    return &it->second;
}

core::ReferenceFrame calculate_base_frame(const core::Residue& residue,
                                          const BaseTemplateLibrary& lib) {
    core::ReferenceFrame frame;
    if (!passes_type_check(residue)) return frame;
    const auto* tmpl_ring = lib.ring_atoms(residue.base_type());
    if (tmpl_ring == nullptr) return frame;

    std::vector<geometry::Vector3d> standard;
    std::vector<geometry::Vector3d> experimental;
    standard.reserve(kRingTypes.size());
    experimental.reserve(kRingTypes.size());
    for (AtomType t : kRingTypes) {
        const core::Atom* atom = ring_atom(residue, t);
        if (atom == nullptr) continue;
        auto it = tmpl_ring->find(t);
        if (it == tmpl_ring->end()) continue;
        standard.push_back(it->second);
        experimental.push_back(atom->coords);
    }
    if (experimental.size() < 3) return frame;

    const auto fit = geometry::superpose(standard, experimental);
    if (!fit.valid) return frame;

    frame.origin = fit.translation;
    frame.rotation = fit.rotation;
    frame.rms_fit = fit.rms;
    frame.valid = true;
    return frame;
}

std::vector<std::pair<std::string, core::ReferenceFrame>>
calculate_all_frames(const core::Structure& structure,
                     const BaseTemplateLibrary& lib) {
    std::vector<std::pair<std::string, core::ReferenceFrame>> frames;
    for (const auto& chain : structure.chains()) {
        for (const auto& res : chain.residues()) {
            core::ReferenceFrame f = calculate_base_frame(res, lib);
            if (f.valid) frames.emplace_back(res.res_id(), f);
        }
    }
    return frames;
}

}  // namespace pairfinder::algorithms
