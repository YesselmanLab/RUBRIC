/**
 * @file scorer.cpp
 * @brief Empirical pair scorer.
 *
 * Severity uses the Cerny method (Cerny et al. NAR 2026 gkaf1335):
 *   severity = 0              if ProSco >= 5  (Preferred tier)
 *   severity = min(1, |Z'|/5) if ProSco <  5
 * Z' is an asymmetric non-parametric standard score loaded from z_tables.json.
 * The 10 individual parameters are scored independently (shear, stretch, stagger,
 * buckle, propeller, opening, distance, hbond_angles, incorrect_hbond_count,
 * clash).
 */
#include <pairfinder/scoring/scorer.hpp>

#include <pairfinder/scoring/clashes.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

#include <pairfinder/algorithms/rna_chains.hpp>
#include <pairfinder/algorithms/selection.hpp>
#include <pairfinder/scoring/hbond_angles.hpp>
#include <pairfinder/scoring/pair_parameters.hpp>
#include <pairfinder/scoring/richardson.hpp>
#include <pairfinder/scoring/torsions.hpp>

namespace pairfinder::scoring {

namespace {

using json = nlohmann::json;
using geometry::Vector3d;

// ---- Sugar-pucker checks (Richardson/Kapral): the base-phosphate perpendicular
//      (Pperp) and the two report-only categories built on it. ----
constexpr double kPperpCutoff = 2.9;   // A: >= C3'-endo, < C2'-endo
const std::unordered_set<std::string> kHelicalWcBp = {"G-C", "C-G", "A-U", "U-A"};

/// Base-phosphate perpendicular distance (A): perpendicular from the 3' P (the
/// next residue's phosphorus) to the extended C1'->N glycosidic vector. Long
/// (>=2.9) = C3'-endo, short (<2.9) = C2'-endo. nullopt if any atom is missing.
std::optional<double> compute_pperp(const core::Residue* residue,
                                    const core::Residue* next_res) {
    if (!residue || !next_res) return std::nullopt;
    const core::Atom* c1 = residue->get_atom("C1'");
    const core::Atom* n  = residue->glycosidic_n();
    const core::Atom* p  = next_res->get_atom("P");
    if (!c1 || !n || !p) return std::nullopt;
    const Vector3d a = c1->coords;
    const Vector3d axis = n->coords - a;
    const double norm = axis.norm();
    if (norm < 1e-6) return std::nullopt;
    const Vector3d d = axis * (1.0 / norm);
    const Vector3d v = p->coords - a;
    const Vector3d perp = v - d * v.dot(d);
    return perp.norm();
}

int pucker_from_pperp(double pperp) { return pperp >= kPperpCutoff ? 3 : 2; }
const char* pucker_label(int p) { return p == 3 ? "C3'-endo" : "C2'-endo"; }

// Cerny severity thresholds.
constexpr double kProscoPreferred  = 5.0;   // ProSco >= 5 -> Preferred, severity = 0
constexpr double kZprimeThreshold  = 5.0;   // |Z'| = 5 -> severity = 1

// Backbone recommendation thresholds.
constexpr double kLowSuitenessMax     = 0.5;  // suiteness below -> residue flagged
constexpr double kFallbackWidthFactor = 1.0;  // no ProSco cell -> fires beyond 1 width

const std::unordered_set<std::string> kBaseAtoms = {
    "N1", "N2", "N3", "N4", "N6", "N7", "N9", "O2", "O4", "O6"};
const std::unordered_set<std::string> kSelfReciprocal = {
    "cWW", "tWW", "cHH", "tHH", "cSS", "tSS"};
const std::unordered_set<std::string> kSignFlip = {"shear", "stagger", "buckle"};

// H-bond distance policy. Watson-Crick cells judge
// every base-base bond; all other cells judge only the strongest (primary) bond.
const std::unordered_set<std::string> kWcDistanceCells = {"C-G/cWW", "A-U/cWW"};
constexpr double kMinPhysicalDistance = 2.2;  // below this: not a real H-bond
constexpr double kMaxStrongDistance   = 3.5;  // strong (primary) bond cutoff

// 10 individual parameters.
const std::array<const char*, 10> kIssueKeys = {
    "shear", "stretch", "stagger", "buckle", "propeller", "opening",
    "distance", "hbond_angles", "incorrect_hbond_count", "clash"};

// Canonical bp type -> (canonical_bp, swapped).
std::pair<std::string, bool> canonicalize_bp(const std::string& bp, const std::string& lw) {
    if (!kSelfReciprocal.count(lw)) return {bp, false};
    const auto dash = bp.find('-');
    if (dash == std::string::npos) return {bp, false};
    const std::string left = bp.substr(0, dash), right = bp.substr(dash + 1);
    if (left <= right) return {bp, false};
    return {right + "-" + left, true};
}

double maybe_flip(double v, const std::string& source, const std::string& param, bool swapped) {
    if (!swapped) return v;
    if (source == "pairs" && kSignFlip.count(param)) return -v;
    return v;
}

// True if (canonical bp_type, lw_class) is a Watson-Crick distance cell.
bool is_wc_distance_cell(const std::string& bp, const std::string& lw) {
    const auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
    (void)swapped;
    return kWcDistanceCells.count(canon_bp + "/" + lw) > 0;
}

// The base-base distance(s)
// whose severity should be taken (max downstream). WC cells return every
// physical bond; other cells return only the single strongest bond within
// [kMinPhysicalDistance, kMaxStrongDistance], or empty when none qualifies.
std::vector<double> distances_to_score(const std::vector<double>& distances,
                                       const std::string& bp, const std::string& lw) {
    std::vector<double> phys;
    for (const double d : distances)
        if (d >= kMinPhysicalDistance) phys.push_back(d);
    if (phys.empty()) return {};
    if (is_wc_distance_cell(bp, lw)) return phys;
    std::vector<double> strong;
    for (const double d : phys)
        if (d <= kMaxStrongDistance) strong.push_back(d);
    if (strong.empty()) return {};
    return {*std::min_element(strong.begin(), strong.end())};
}

/// severity = min(1, |Z'|/5) * clamp((5 - ProSco)/5, 0, 1)
///
/// The ProSco ramp keeps severity continuous at the boundary, where plain Cerny
/// steps from 0 to ~0.5. Reduces to plain
/// Cerny once ProSco hits 0. No ProSco -> 0. No Z' -> the ramp alone.
/// Severity from a ProSco / Z' pair, per the Cerny three-tier scheme.
double severity_from_prosco_zprime(std::optional<double> prosco,
                                   std::optional<double> zprime) {
    if (!prosco || *prosco >= kProscoPreferred) return 0.0;
    const double ramp =
        std::min(1.0, std::max(0.0, (kProscoPreferred - *prosco) / kProscoPreferred));
    if (!zprime) return ramp;
    return std::min(1.0, std::abs(*zprime) / kZprimeThreshold) * ramp;
}

/// Severity plus the ProSco/Z' that produced it (kept for per-issue reporting).
struct SevDetail {
    double severity = 0.0;
    std::optional<double> prosco;
    std::optional<double> zprime;
};

/// Cerny three-tier label (Cerny et al. NAR 2026 gkaf1335):
///   Preferred   = ProSco >= 5
///   Allowed     = ProSco < 5 and |Z'| <= 5
///   Of Concern  = ProSco < 5 and |Z'| > 5  (or Z' unavailable)
/// Returns "" when the parameter is not ProSco/Z'-scored (e.g. hbond count,
/// where prosco is none) - triggered issues of that kind carry an empty tier.
std::string prosco_tier(std::optional<double> prosco, std::optional<double> zprime) {
    if (!prosco) return "";
    if (*prosco >= kProscoPreferred) return "Preferred";
    if (!zprime) return "Of Concern";
    return std::abs(*zprime) <= kZprimeThreshold ? "Allowed" : "Of Concern";
}

// H-bond angle geometry (angle_at_atom / hbond_angle) lives in the shared
// header scoring/hbond_angles.hpp so the CLI per-pair emitter uses it too.

struct HBF {
    double distance;
    std::optional<double> donor_angle, acceptor_angle;
};

// ---------------------------------------------------------------------------
//  Scoring tables in native form
//
// The scorer reads three things from its JSON tables: ProSco bins
// (prosco_distributions.json), Z' parameters (z_tables.json) and, from
// parameter_distributions.json, only the most frequent base-base H-bond count per
// class. Parsing ~1.6 MB of JSON into a document tree took ~15-20 ms per run, as
// long as scoring a small RNA; these flat tables hold exactly the values the
// lookups read (the same doubles, taken from the parsed JSON), so every score is
// unchanged. They are loaded from scoring_tables_cache.bin (written by `rubric
// dump-caches`) when its stored content hash matches the JSON files, and built
// from the JSON otherwise, so a stale cache can only make a run slower.
// ---------------------------------------------------------------------------

struct ProscoCell {
    double lo = 0.0, hi = 0.0;
    int n_bins = 0;
    std::vector<double> per_bin;
};

struct ZCell {
    double M = 0.0, sl = 0.0, su = 0.0;
};

/// Cell keys are source + '\x01' + param + '\x01' + class key ("GC/cWW", "_ANY/cWW", ...).
struct ScoringTables {
    std::unordered_map<std::string, ProscoCell> prosco;
    std::unordered_map<std::string, ZCell> z;
    std::unordered_map<std::string, int> hbond_count;   ///< class key -> expected count
};

std::string cell_key(const std::string& source, const std::string& param, const std::string& key) {
    std::string k;
    k.reserve(source.size() + param.size() + key.size() + 2);
    k.append(source).append(1, '\x01').append(param).append(1, '\x01').append(key);
    return k;
}

/// The tables as read from the JSON: a ProSco cell only if
/// it has a non-empty prosco_per_bin, a Z' cell only if it has "M", an H-bond
/// count only where "frequencies" is non-empty (argmax in key order, first wins).
ScoringTables tables_from_json(const json& prosco, const json& z_tables, const json& categorical) {
    ScoringTables t;
    if (prosco.is_object())
        for (const auto& [source, sv] : prosco.items()) {
            if (!sv.is_object()) continue;
            for (const auto& [param, pv] : sv.items()) {
                if (!pv.is_object()) continue;
                for (const auto& [key, cell] : pv.items()) {
                    if (!cell.is_object() || !cell.contains("prosco_per_bin") ||
                        cell["prosco_per_bin"].empty())
                        continue;
                    ProscoCell c;
                    c.lo = cell["support_lo"].get<double>();
                    c.hi = cell["support_hi"].get<double>();
                    c.n_bins = cell["n_bins"].get<int>();
                    for (const auto& v : cell["prosco_per_bin"])
                        c.per_bin.push_back(v.get<double>());
                    t.prosco.emplace(cell_key(source, param, key), std::move(c));
                }
            }
        }
    if (z_tables.is_object())
        for (const auto& [source, sv] : z_tables.items()) {
            if (!sv.is_object()) continue;
            for (const auto& [param, pv] : sv.items()) {
                if (!pv.is_object()) continue;
                for (const auto& [key, cell] : pv.items()) {
                    if (!cell.is_object() || !cell.contains("M")) continue;
                    t.z.emplace(cell_key(source, param, key),
                                ZCell{cell["M"].get<double>(), cell["sl"].get<double>(),
                                      cell["su"].get<double>()});
                }
            }
        }
    if (categorical.is_object() && categorical.contains("pairs") &&
        categorical["pairs"].contains("num_base_hbonds"))
        for (const auto& [key, cell] : categorical["pairs"]["num_base_hbonds"].items()) {
            if (!cell.is_object() || !cell.contains("frequencies")) continue;
            const auto& freqs = cell["frequencies"];
            if (freqs.empty()) continue;
            std::string best_k;
            double best_v = -1.0;
            for (auto fit = freqs.begin(); fit != freqs.end(); ++fit)
                if (fit.value().get<double>() > best_v) { best_v = fit.value().get<double>(); best_k = fit.key(); }
            int count = 1;
            try { count = std::stoi(best_k); } catch (...) { count = 1; }
            t.hbond_count.emplace(key, count);
        }
    return t;
}

/// A whole file in one read (empty if it cannot be opened).
std::string read_bytes(const std::filesystem::path& path) {
    std::string out;
    std::FILE* f = std::fopen(path.string().c_str(), "rb");
    if (f == nullptr) return out;
    if (std::fseek(f, 0, SEEK_END) == 0) {
        const long n = std::ftell(f);
        if (n > 0 && std::fseek(f, 0, SEEK_SET) == 0) {
            out.resize(static_cast<std::size_t>(n));
            out.resize(std::fread(out.data(), 1, out.size(), f));
        }
    }
    std::fclose(f);
    return out;
}

/// 64-bit hash of the given files' names, sizes and bytes (8 bytes per step), so
/// any edit to a source table changes it and retires the cache.
std::uint64_t hash_files(const std::vector<std::filesystem::path>& files) {
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](std::uint64_t w) {
        h ^= w + 0x9E3779B97F4A7C15ULL + (h << 6) + (h >> 2);
        h *= 0xBF58476D1CE4E5B9ULL;
        h ^= h >> 31;
    };
    const auto mix_bytes = [&mix](const char* data, std::size_t n) {
        std::size_t i = 0;
        for (; i + 8 <= n; i += 8) {
            std::uint64_t w;
            std::memcpy(&w, data + i, 8);
            mix(w);
        }
        std::uint64_t tail = 0;
        std::memcpy(&tail, data + i, n - i);
        mix(tail ^ (static_cast<std::uint64_t>(n - i) << 56));
    };
    for (const auto& f : files) {
        const std::string name = f.filename().string();
        mix_bytes(name.data(), name.size());
        const std::string bytes = read_bytes(f);
        mix(static_cast<std::uint64_t>(bytes.size()));
        mix_bytes(bytes.data(), bytes.size());
    }
    return h;
}

constexpr char kTablesMagic[8] = {'P', 'R', 'S', 'T', 'A', 'B', '0', '1'};

/// Binary layout (native byte order; the magic and a probe double reject a
/// mismatched machine): magic, probe, source hash, then each table as a count
/// followed by its entries in sorted key order, so the file is deterministic.
void save_tables(const std::filesystem::path& path, const ScoringTables& t,
                 std::uint64_t source_hash) {
    std::string out;
    const auto raw = [&out](const void* p, std::size_t n) {
        out.append(static_cast<const char*>(p), n);
    };
    const auto str = [&](const std::string& s) {
        const auto n = static_cast<std::uint32_t>(s.size());
        raw(&n, 4);
        raw(s.data(), n);
    };
    const auto u32 = [&](std::uint32_t v) { raw(&v, 4); };
    const auto f64 = [&](double v) { raw(&v, 8); };
    const auto sorted = [](const auto& m) {
        std::vector<std::string> keys;
        for (const auto& kv : m) keys.push_back(kv.first);
        std::sort(keys.begin(), keys.end());
        return keys;
    };
    raw(kTablesMagic, 8);
    f64(1.5);
    raw(&source_hash, 8);
    u32(static_cast<std::uint32_t>(t.prosco.size()));
    for (const auto& k : sorted(t.prosco)) {
        const ProscoCell& c = t.prosco.at(k);
        str(k); f64(c.lo); f64(c.hi); u32(static_cast<std::uint32_t>(c.n_bins));
        u32(static_cast<std::uint32_t>(c.per_bin.size()));
        raw(c.per_bin.data(), c.per_bin.size() * sizeof(double));
    }
    u32(static_cast<std::uint32_t>(t.z.size()));
    for (const auto& k : sorted(t.z)) {
        const ZCell& c = t.z.at(k);
        str(k); f64(c.M); f64(c.sl); f64(c.su);
    }
    u32(static_cast<std::uint32_t>(t.hbond_count.size()));
    for (const auto& k : sorted(t.hbond_count)) {
        str(k);
        u32(static_cast<std::uint32_t>(t.hbond_count.at(k)));
    }
    std::ofstream(path, std::ios::binary) << out;
}

/// The cached tables, or nullopt if the file is absent, malformed, from another
/// byte order, or built from different source tables.
std::optional<ScoringTables> load_tables(const std::filesystem::path& path,
                                         std::uint64_t source_hash) {
    const std::string buf = read_bytes(path);
    if (buf.empty()) return std::nullopt;
    std::size_t pos = 0;
    bool ok = true;
    const auto take = [&](void* dst, std::size_t n) {
        if (!ok || pos + n > buf.size()) { ok = false; return; }
        std::memcpy(dst, buf.data() + pos, n);
        pos += n;
    };
    const auto u32 = [&]() { std::uint32_t v = 0; take(&v, 4); return v; };
    const auto f64 = [&]() { double v = 0; take(&v, 8); return v; };
    const auto str = [&]() {
        const std::uint32_t n = u32();
        if (!ok || pos + n > buf.size()) { ok = false; return std::string(); }
        std::string s(buf.data() + pos, n);
        pos += n;
        return s;
    };
    char magic[8];
    take(magic, 8);
    if (!ok || std::memcmp(magic, kTablesMagic, 8) != 0 || f64() != 1.5) return std::nullopt;
    std::uint64_t stored = 0;
    take(&stored, 8);
    if (!ok || stored != source_hash) return std::nullopt;
    ScoringTables t;
    for (std::uint32_t i = 0, n = u32(); ok && i < n; ++i) {
        if (i == 0) t.prosco.reserve(n);
        std::string k = str();
        ProscoCell c;
        c.lo = f64(); c.hi = f64(); c.n_bins = static_cast<int>(u32());
        c.per_bin.resize(u32());
        take(c.per_bin.data(), c.per_bin.size() * sizeof(double));
        t.prosco.emplace(std::move(k), std::move(c));
    }
    for (std::uint32_t i = 0, n = u32(); ok && i < n; ++i) {
        if (i == 0) t.z.reserve(n);
        std::string k = str();
        ZCell c;
        c.M = f64(); c.sl = f64(); c.su = f64();
        t.z.emplace(std::move(k), c);
    }
    for (std::uint32_t i = 0, n = u32(); ok && i < n; ++i) {
        std::string k = str();
        t.hbond_count.emplace(std::move(k), static_cast<int>(u32()));
    }
    if (!ok || pos != buf.size()) return std::nullopt;
    return t;
}

}  // namespace

struct Scorer::Impl {
    ScoringTables tables;   // ProSco bins, Z' parameters, expected H-bond counts
    std::vector<std::filesystem::path> table_sources;   // the JSON the tables come from
    std::filesystem::path tables_cache;                  // scoring_tables_cache.bin
    std::uint64_t sources_hash = 0;
    // backbone prosco_distributions["backbone"] (per-conformer torsion). Only a
    // flagged residue's recommendation reads it, and many structures have none,
    // so it is parsed on first use rather than for every run.
    mutable json bb_prosco;
    mutable bool bb_loaded = false;
    std::filesystem::path bb_path;
    const json& backbone_prosco() const {
        if (!bb_loaded) {
            bb_loaded = true;
            if (std::filesystem::exists(bb_path)) {
                std::ifstream in(bb_path);
                json d = json::parse(in);
                if (d.contains("prosco_distributions") &&
                    d["prosco_distributions"].contains("backbone"))
                    bb_prosco = std::move(d["prosco_distributions"]["backbone"]);
            }
        }
        return bb_prosco;
    }
    std::unordered_map<std::string, double> pair_weights;
    std::unique_ptr<RichardsonClassifier> richardson;

    /// ProSco bin lookup with bp/lw -> _ANY/lw -> _ANY/_ANY fallback.
    std::optional<double> prosco_lookup(std::optional<double> value, const std::string& source,
                                        const std::string& param, const std::string& bp,
                                        const std::string& lw) const {
        if (!value || std::isnan(*value)) return std::nullopt;
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        const double v = maybe_flip(*value, source, param, swapped);
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = tables.prosco.find(cell_key(source, param, key));
            if (it == tables.prosco.end()) continue;
            const ProscoCell& c = it->second;
            int idx = static_cast<int>((v - c.lo) / (c.hi - c.lo) * c.n_bins);
            idx = std::max(0, std::min(c.n_bins - 1, idx));
            return c.per_bin.at(static_cast<std::size_t>(idx));
        }
        return std::nullopt;
    }

    /// Z' asymmetric standard score with bp/lw -> _ANY/lw -> _ANY/_ANY fallback.
    /// Z' = (value - M) / sigma, where sigma = sl (lower half) or su (upper half).
    std::optional<double> zprime_lookup(std::optional<double> value, const std::string& source,
                                        const std::string& param, const std::string& bp,
                                        const std::string& lw) const {
        if (!value || std::isnan(*value)) return std::nullopt;
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        const double v = maybe_flip(*value, source, param, swapped);
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = tables.z.find(cell_key(source, param, key));
            if (it == tables.z.end()) continue;
            const double M  = it->second.M;
            const double sl = it->second.sl;
            const double su = it->second.su;
            if (v >= M) return su > 1e-9 ? std::optional<double>((v - M) / su) : std::nullopt;
            return sl > 1e-9 ? std::optional<double>((v - M) / sl) : std::nullopt;
        }
        return std::nullopt;
    }

    int expected_hbond_count(const std::string& bp, const std::string& lw) const {
        auto [canon_bp, swapped] = canonicalize_bp(bp, lw);
        (void)swapped;
        for (const std::string& key : {canon_bp + "/" + lw, "_ANY/" + lw, std::string("_ANY/_ANY")}) {
            const auto it = tables.hbond_count.find(key);
            if (it != tables.hbond_count.end()) return it->second;
        }
        return 1;
    }

    double hbond_count_severity(const std::string& bp, const std::string& lw, int actual) const {
        const int canonical = expected_hbond_count(bp, lw);
        if (canonical <= 0 || actual >= canonical) return 0.0;
        return static_cast<double>(canonical - actual) / canonical;
    }

    // --- backbone recommendation (per-conformer ProSco) ---

    /// ProSco of a suite torsion against its empirical per-conformer cell, with
    /// circular unwrap onto the branch used at build time (unwrap_center).
    std::optional<double> bb_prosco_lookup(const std::string& angle,
                                           const std::string& conformer,
                                           double value) const {
        if (conformer.empty()) return std::nullopt;
        const json& bb = backbone_prosco();
        if (bb.is_null()) return std::nullopt;
        const auto ait = bb.find(angle);
        if (ait == bb.end()) return std::nullopt;
        const auto cit = ait->find(conformer);
        if (cit == ait->end()) return std::nullopt;
        const json& cell = *cit;
        if (!cell.contains("prosco_per_bin") || cell["prosco_per_bin"].empty())
            return std::nullopt;
        const double uc = cell["unwrap_center"].get<double>();
        const double vu = uc + signed_angle_gap(value, uc);
        const double lo = cell["support_lo"].get<double>();
        const double hi = cell["support_hi"].get<double>();
        const int n = cell["n_bins"].get<int>();
        int idx = static_cast<int>((vu - lo) / (hi - lo) * n);
        idx = std::max(0, std::min(n - 1, idx));
        return cell["prosco_per_bin"][idx].get<double>();
    }

    /// How far along the way to the conformer centre a recommendation should go.
    ///
    /// A full snap (1.0) would put the recommendation at the cluster centre
    /// itself. That overshoots: real RNA sits *near* conformer
    /// centres with genuine spread, so a distorted-but-real coordinate is usually
    /// closer to its own true value than the idealised centre is.
    ///
    /// Fitted against independently determined better models of the same molecule
    /// by leave-one-cluster-out cross-validation, with separate values per tier.
    ///
    /// Above kFullSnapGap the deposited value carries no usable information and
    /// going all the way to the centre is right.
    ///
    /// Applied by damped_target() when building a recommendation.
    static constexpr double kDampingFixable = 0.45;
    static constexpr double kDampingOther   = 0.70;
    static constexpr double kDampingDefault = 0.55;
    static constexpr double kFullSnapGap    = 90.0;

    static double damped_target(double value, double gap, const std::string& tier) {
        double alpha = kDampingDefault;
        if (std::abs(gap) > kFullSnapGap) {
            alpha = 1.0;
        } else if (tier == "fixable") {
            alpha = kDampingFixable;
        } else if (tier == "flag_only") {
            alpha = kDampingOther;
        }
        double t = std::fmod(value - alpha * gap, 360.0);
        if (t < 0.0) t += 360.0;
        return t;
    }

    /// Per-residue backbone recommendation. deviations = only torsions that fire
    /// (ProSco < 5), ranked most-anomalous first. tier fixable / review / flag_only.
    std::optional<BackboneRecommendation> build_recommendation(
            const std::array<double, 7>& suite, const SuiteResult& sr) const {
        std::string tier, target;
        if (sr.is_outlier) {
            tier = "flag_only";
            target = richardson->nearest_conformer(suite);
        } else if (sr.suiteness < kLowSuitenessMax) {
            tier = "fixable";
            target = sr.conformer;
        } else {
            return std::nullopt;
        }
        BackboneRecommendation rec;
        rec.target_conformer = target;
        const std::array<double, 7>* center =
            target.empty() ? nullptr : richardson->center_for(target);
        if (center) {
            const auto& widths = richardson->normal_widths();
            for (int k = 0; k < 7; ++k) {
                const double val = suite[k];
                const double gap = signed_angle_gap(val, (*center)[k]);
                const std::string angle = kSuiteAngleNames[k];
                const std::optional<double> ps = bb_prosco_lookup(angle, target, val);
                const bool fires = ps ? (*ps < kProscoPreferred)
                                      : (std::abs(gap) > kFallbackWidthFactor * widths[k]);
                if (fires) {
                    AngleDeviation d;
                    d.angle = angle; d.value = val;
                    d.target = damped_target(val, gap, tier);
                    d.gap = gap; d.prosco = ps;
                    d.conformer_center = (*center)[k];
                    rec.deviations.push_back(std::move(d));
                }
            }
            std::sort(rec.deviations.begin(), rec.deviations.end(),
                      [](const AngleDeviation& a, const AngleDeviation& b) {
                          const double ka = a.prosco ? *a.prosco : 1000.0 - std::abs(a.gap);
                          const double kb = b.prosco ? *b.prosco : 1000.0 - std::abs(b.gap);
                          return ka < kb;
                      });
        }
        if (tier == "fixable" && rec.deviations.empty()) tier = "review";
        rec.tier = tier;
        return rec;
    }
};

Scorer::Scorer(const std::filesystem::path& distributions_path,
               const std::filesystem::path& prosco_path,
               const std::filesystem::path& weights_path,
               const std::filesystem::path& richardson_suites_path)
    : impl_(std::make_unique<Impl>()) {
    impl_->richardson = std::make_unique<RichardsonClassifier>(richardson_suites_path);
    // z_tables.json is auto-discovered alongside prosco_distributions.json.
    const auto z_path = prosco_path.parent_path() / "z_tables.json";
    impl_->table_sources = {distributions_path, prosco_path, z_path};
    impl_->tables_cache = prosco_path.parent_path() / "scoring_tables_cache.bin";
    impl_->sources_hash = hash_files(impl_->table_sources);
    if (auto cached = load_tables(impl_->tables_cache, impl_->sources_hash)) {
        impl_->tables = std::move(*cached);
    } else {
        json categorical, prosco, z_tables;
        {
            std::ifstream in(distributions_path);
            json d = json::parse(in);
            if (d.contains("categorical_distributions"))
                categorical = std::move(d["categorical_distributions"]);
        }
        {
            std::ifstream in(prosco_path);
            json d = json::parse(in);
            if (d.contains("prosco_distributions")) prosco = std::move(d["prosco_distributions"]);
        }
        if (std::filesystem::exists(z_path)) {
            std::ifstream in(z_path);
            z_tables = json::parse(in);
        }
        impl_->tables = tables_from_json(prosco, z_tables, categorical);
    }
    // Auto-discover backbone_prosco_distributions.json (per-conformer torsion
    // ProSco used for backbone recommendations).
    impl_->bb_path = prosco_path.parent_path() / "backbone_prosco_distributions.json";
    {
        std::ifstream in(weights_path);
        json d = json::parse(in);
        const json& raw = d["penalty_weights"];
        // Weights already sum to 100 as generated, but renormalise
        // defensively in case the file is from a different generation.
        double total = 0.0;
        for (const char* k : kIssueKeys) total += raw.value(k, 0.0);
        for (const char* k : kIssueKeys)
            impl_->pair_weights[k] = total <= 0 ? 100.0 / kIssueKeys.size()
                                                : raw.value(k, 0.0) * 100.0 / total;
    }
}

void Scorer::save_tables_cache(const std::filesystem::path& path) const {
    save_tables(path.empty() ? impl_->tables_cache : path, impl_->tables, impl_->sources_hash);
}

Scorer::~Scorer() = default;
Scorer::Scorer(Scorer&&) noexcept = default;

namespace {
/// Residue name out of a "<chain>-<name>-<seq>" res_id.
std::string res_name_of(const std::string& res_id) {
    const auto a = res_id.find('-');
    if (a == std::string::npos) return {};
    const auto b = res_id.find('-', a + 1);
    return res_id.substr(a + 1, (b == std::string::npos ? res_id.size() : b) - a - 1);
}
}  // namespace

bool Scorer::is_scorable(const classification::ScoredCandidate& p,
                         const core::BaseTyping& typing) {
    if (p.pair_category != "base-base") return false;
    if (p.is_ambiguous) return false;
    if (p.lw_class.empty() || p.lw_class.find('|') != std::string::npos) return false;
    // The scoring tables are calibrated on RNA, so DNA is not scored - either
    // side DNA excludes the pair, hybrids included. Filtered here rather than
    // at candidate generation so the pair is still emitted with a null score,
    // as ambiguous pairs are. The backbone side already skips DNA.
    if (typing.is_dna(res_name_of(p.res_id1))) return false;
    if (typing.is_dna(res_name_of(p.res_id2))) return false;
    return true;
}

PairScore Scorer::score_pair(const classification::ScoredCandidate& pair, const core::Residue& res1,
                             const core::Residue& res2, hbond::HBondFinder& finder,
                             const hbond::HBondChemistry& chem, const core::BaseTyping& typing) const {
    const PairParameters params = compute_pair_parameters(*pair.frame1, *pair.frame2);
    const std::string bp = typing.normalize(pair.res_name1) + "-" + typing.normalize(pair.res_name2);
    const std::string& lw = pair.lw_class;

    // H-bonds (base-base) with angles.
    std::vector<HBF> hbonds;
    int num_base = 0;
    for (const auto& hb : finder.find_between(res1, res2)) {
        if (!kBaseAtoms.count(hb.donor_atom) || !kBaseAtoms.count(hb.acceptor_atom)) continue;
        const core::Residue* dres = (hb.donor_res_id == pair.res_id1) ? &res1 : &res2;
        const core::Residue* ares = (hb.donor_res_id == pair.res_id1) ? &res2 : &res1;
        const core::Atom* da = dres->get_atom(hb.donor_atom);
        const core::Atom* aa = ares->get_atom(hb.acceptor_atom);
        HBF f{hb.distance, std::nullopt, std::nullopt};
        if (da && aa) {
            f.donor_angle    = hbond_angle(*dres, hb.donor_atom,    da->coords, aa->coords, chem);
            f.acceptor_angle = hbond_angle(*ares, hb.acceptor_atom, aa->coords, da->coords, chem);
        }
        hbonds.push_back(f);
        ++num_base;
    }

    // Helper: Cerny severity for a continuous parameter, plus the ProSco/Z' that
    // produced it (kept so each triggered issue can report its tier).
    auto cerny_sev = [&](std::optional<double> val, const std::string& source,
                         const char* param) -> SevDetail {
        const auto ps = impl_->prosco_lookup(val, source, param, bp, lw);
        const auto zp = impl_->zprime_lookup(val, source, param, bp, lw);
        return {severity_from_prosco_zprime(ps, zp), ps, zp};
    };

    std::unordered_map<std::string, SevDetail> issues;

    // ---- 6 pair geometry parameters, each individually ----
    const std::array<std::pair<const char*, double>, 6> pair_params = {{
        {"shear",     params.shear},
        {"stretch",   params.stretch},
        {"stagger",   params.stagger},
        {"buckle",    params.buckle},
        {"propeller", params.propeller},
        {"opening",   params.opening},
    }};
    for (const auto& [param, val] : pair_params) {
        const auto d = cerny_sev(val, "pairs", param);
        if (d.severity > 0.0) issues[param] = d;
    }

    // ---- H-bond distance: WC cells judge all base-base bonds; other cells
    //      judge only the strong (primary) bond.
    //      Keep the ProSco/Z' of the worst (max-severity) scored bond. ----
    std::vector<double> bond_distances;
    bond_distances.reserve(hbonds.size());
    for (const auto& hb : hbonds) bond_distances.push_back(hb.distance);
    SevDetail dist;
    for (const double d : distances_to_score(bond_distances, bp, lw)) {
        const auto sd = cerny_sev(d, "hbonds", "distance");
        if (sd.severity > dist.severity) dist = sd;
    }
    if (dist.severity > 0.0) issues["distance"] = dist;

    // ---- H-bond angles: max over donor + acceptor angles across all bonds;
    //      keep the ProSco/Z' of the worst (argmax) angle. ----
    SevDetail ang;
    for (const auto& hb : hbonds) {
        const auto sd_d = cerny_sev(hb.donor_angle,    "hbonds", "donor_angle");
        if (sd_d.severity > ang.severity) ang = sd_d;
        const auto sd_a = cerny_sev(hb.acceptor_angle, "hbonds", "acceptor_angle");
        if (sd_a.severity > ang.severity) ang = sd_a;
    }
    if (ang.severity > 0.0) issues["hbond_angles"] = ang;

    // ---- H-bond count: graded (canonical - actual) / canonical (no ProSco/Z') ----
    {
        const double s = impl_->hbond_count_severity(bp, lw, num_base);
        if (s > 0.0) issues["incorrect_hbond_count"] = SevDetail{s, std::nullopt, std::nullopt};
    }

    // ---- Steric clash: graded on heavy-atom overlap (no ProSco/Z') ----
    // Gate on kClashThreshold: the ungated overlap charges score for ordinary
    // vdW proximity. worst_overlap() returns the contact pair_clash() gates on,
    // so this is equivalent and computes it once. Raw value kept on PairScore.
    const std::optional<double> clash_overlap = worst_overlap(res1, res2, chem);
    if (clash_overlap && *clash_overlap > kClashThreshold) {
        const double s = clash_severity(*clash_overlap);
        if (s > 0.0) issues["clash"] = SevDetail{s, std::nullopt, std::nullopt};
    }

    // Build weighted penalty list in kIssueKeys order, then sort by descending penalty.
    std::vector<IssuePenalty> weighted;
    double penalty = 0.0;
    for (const char* k : kIssueKeys) {
        const auto it = issues.find(k);
        if (it == issues.end()) continue;
        const SevDetail& d = it->second;
        // A weights file with no entry for a parameter: treat it as weight 0
        // (detected and reported, but costs no score) rather than throwing.
        const auto wit = impl_->pair_weights.find(k);
        const double w = (wit == impl_->pair_weights.end() ? 0.0 : wit->second) * d.severity;
        IssuePenalty ip;
        ip.issue = k;
        ip.weight = w;
        ip.severity = d.severity;
        ip.prosco = d.prosco;
        ip.zprime = d.zprime;
        ip.tier = prosco_tier(d.prosco, d.zprime);
        weighted.push_back(std::move(ip));
        penalty += w;
    }
    std::stable_sort(weighted.begin(), weighted.end(),
                     [](const IssuePenalty& a, const IssuePenalty& b) {
                         return a.weight > b.weight;
                     });

    PairScore out;
    out.res_id1 = pair.res_id1;
    out.res_id2 = pair.res_id2;
    out.bp_type = bp;
    out.lw_class = lw;
    out.score   = std::max(0.0, std::min(100.0, 100.0 - penalty));
    out.penalty = penalty;
    out.issues  = std::move(weighted);
    out.clash_overlap = clash_overlap;
    return out;
}

namespace {

// Longest O3'(i-1) - P(i) accepted as the phosphodiester bond joining two
// residues into one suite. Real links are extremely tight (~1.60 A); anything
// above 2.0 A is not a stretched bond but residues that are not joined at all.
constexpr double kMaxPhosphodiester = 2.0;

// Map each res_id to its backbone predecessor's res_id: the residue
// immediately preceding it within the same chain, ordered by (seq, icode),
// when their integer sequence numbers differ by at most one (equal when only
// the insertion code advances, e.g. 20 -> 20A -> 21) AND the two are actually
// bonded. Numbering alone is not enough: some entries label separate strands
// with the same numbers, which would build suites between residues that are
// far apart.
std::unordered_map<std::string, std::string> build_predecessor_index(
    const std::unordered_map<std::string, ResidueTorsions>& tors,
    const std::unordered_map<std::string, const core::Residue*>& res_by_id) {
    std::map<std::string, std::vector<std::tuple<long, std::string, std::string>>> by_chain;
    for (const auto& [rid, t] : tors) {
        (void)t;
        const auto p = parse_res_seq(rid);
        if (!p) continue;
        by_chain[p->chain].emplace_back(p->num, p->icode, rid);
    }
    std::unordered_map<std::string, std::string> predecessor;
    for (auto& [chain, residues] : by_chain) {
        (void)chain;
        std::sort(residues.begin(), residues.end());
        for (std::size_t i = 1; i < residues.size(); ++i) {
            const long prev_num = std::get<0>(residues[i - 1]);
            const long cur_num  = std::get<0>(residues[i]);
            if (prev_num != cur_num && prev_num != cur_num - 1) continue;
            const std::string& prev_id = std::get<2>(residues[i - 1]);
            const std::string& cur_id  = std::get<2>(residues[i]);
            const auto pit = res_by_id.find(prev_id);
            const auto cit = res_by_id.find(cur_id);
            if (pit == res_by_id.end() || cit == res_by_id.end()) continue;
            const auto* o3 = pit->second->get_atom("O3'");
            const auto* ph = cit->second->get_atom("P");
            if (!o3 || !ph) continue;
            if ((o3->coords - ph->coords).norm() > kMaxPhosphodiester) continue;
            predecessor[cur_id] = prev_id;
        }
    }
    return predecessor;
}

double mod360(double v) {
    double a = std::fmod(v, 360.0);
    if (a < 0) a += 360.0;
    return a;
}

// build_suite_for: [delta_prev, eps_prev, zeta_prev, alpha, beta, gamma, delta].
std::optional<std::array<double, 7>> build_suite(
    const std::string& res_id, const std::unordered_map<std::string, ResidueTorsions>& tors,
    const std::unordered_map<std::string, std::string>& pred_index) {
    const auto pit = pred_index.find(res_id);
    if (pit == pred_index.end()) return std::nullopt;
    const auto cit  = tors.find(res_id);
    const auto prit = tors.find(pit->second);
    if (cit == tors.end() || prit == tors.end()) return std::nullopt;
    const ResidueTorsions& cur  = cit->second;
    const ResidueTorsions& pred = prit->second;
    const std::array<std::optional<double>, 7> n = {pred.delta, pred.epsilon, pred.zeta,
                                                     cur.alpha, cur.beta, cur.gamma, cur.delta};
    std::array<double, 7> suite{};
    for (int k = 0; k < 7; ++k) {
        if (!n[k]) return std::nullopt;
        suite[k] = mod360(*n[k]);
    }
    return suite;
}

double round2(double v) { return std::round(v * 100.0) / 100.0; }

// Longest C2'-O2' distance accepted as a bond (the same cutoff as the backbone
// measurement in MolProbity's dangle, `rnabb`).
constexpr double kMaxC2O2 = 2.0;

// A residue carries an RNA backbone suite only if it is a polymer nucleotide
// whose ribose is modelled with its 2' oxygen: its C2' must be bonded to an O2'.
// DNA (no 2' oxygen), RNA residues modelled without their O2', and 2'-fluoro or
// 2'-seleno nucleotides therefore have no suite; the wwPDB validation reports
// give no suite for such residues either. A free nucleotide or cofactor (ATP,
// GTP, ... typed LIGAND in the ligand DB) is not a polymer residue at all: it may
// still pair, which is why it is parsed, but it has no place in an RNA backbone
// score.
bool has_rna_suite(const std::string& res_id, const core::Residue* res,
                   const core::BaseTyping& typing) {
    if (res == nullptr) return false;
    const auto a = res_id.find('-');
    if (a != std::string::npos) {
        const auto b = res_id.find('-', a + 1);
        const std::string name =
            res_id.substr(a + 1, (b == std::string::npos ? res_id.size() : b) - a - 1);
        if (typing.is_ligand(name)) return false;
    }
    const core::Atom* c2 = res->get_atom("C2'");
    if (c2 == nullptr) return false;
    for (const char* name : {"O2'", "O2*"}) {
        const core::Atom* o2 = res->get_atom(name);
        if (o2 != nullptr && c2->distance_to(*o2) < kMaxC2O2) return true;
    }
    return false;
}

}  // namespace

StructureScore Scorer::score_structure(
    const std::vector<classification::ScoredCandidate>& selected, const core::Structure& structure,
    hbond::HBondFinder& finder, const hbond::HBondChemistry& chem,
    const core::BaseTyping& typing, bool collect_backbone) const {
    std::unordered_map<std::string, const core::Residue*> res_by_id;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) res_by_id[res.res_id()] = &res;

    const auto tors       = compute_all_torsions(structure);
    const auto pred_index = build_predecessor_index(tors, res_by_id);
    // 3'-successor (invert predecessor) for the Pperp phosphate, and helix
    // membership over the selected pairs - both drive the report-only sugar
    // pucker flags (do not affect the score).
    std::unordered_map<std::string, std::string> succ_index;
    succ_index.reserve(pred_index.size());
    for (const auto& [res, pred] : pred_index) succ_index[pred] = res;
    const auto chains = pairfinder::algorithms::RNAChains::from_structure(structure);
    const auto in_helix = pairfinder::algorithms::selection::helix_membership(selected, chains);

    StructureScore ss;
    for (std::size_t i = 0; i < selected.size(); ++i) {
        const auto& sc = selected[i];
        if (!is_scorable(sc, typing)) { ++ss.skipped_pairs; continue; }
        const auto r1 = res_by_id.find(sc.res_id1);
        const auto r2 = res_by_id.find(sc.res_id2);
        if (r1 == res_by_id.end() || r2 == res_by_id.end()) continue;
        PairScore ps = score_pair(sc, *r1->second, *r2->second, finder, chem, typing);
        // "Unusual sugar pucker": C2'-endo (by delta) in a helical cWW G-C/A-U pair.
        if (in_helix[i] && ps.lw_class == "cWW" && kHelicalWcBp.count(ps.bp_type)) {
            for (const core::Residue* res : {r1->second, r2->second}) {
                // DNA is natively C2'-endo; a ribose without its O2' is not judged.
                if (!has_rna_suite(res->res_id(), res, typing)) continue;
                const auto tit = tors.find(res->res_id());
                if (tit == tors.end() || !tit->second.delta) continue;
                const double delta = *tit->second.delta;
                if (impl_->richardson->pucker(delta) != 2) continue;  // only C2'-endo
                const core::Residue* next_res = nullptr;
                const auto sit = succ_index.find(res->res_id());
                if (sit != succ_index.end()) {
                    const auto nit = res_by_id.find(sit->second);
                    if (nit != res_by_id.end()) next_res = nit->second;
                }
                const auto pp = compute_pperp(res, next_res);
                UnusualPucker up;
                up.res_id = res->res_id();
                up.delta = delta;
                up.pperp = pp;
                up.pperp_confirms = (pp && *pp < kPperpCutoff);
                ps.unusual_pucker.push_back(std::move(up));
            }
        }
        ss.pair_scores.push_back(std::move(ps));
    }
    for (const auto& [rid, t] : tors) {
        const auto suite = build_suite(rid, tors, pred_index);
        if (!suite) continue;
        const SuiteResult sr = impl_->richardson->classify(*suite);
        const auto residue_of = [&](const std::string& id) -> const core::Residue* {
            const auto it = res_by_id.find(id);
            return it == res_by_id.end() ? nullptr : it->second;
        };
        // RNA-only backbone score (see has_rna_suite). The suite's delta, epsilon
        // and zeta come from the PREDECESSOR, so a residue whose 5' neighbour is
        // not an RNA nucleotide with its O2' (a DNA sugar, a ribose modelled
        // without O2', or an initiating free nucleotide such as ATP) has
        // no RNA suite either.
        bool scored = has_rna_suite(rid, residue_of(rid), typing);
        if (scored) {
            if (const auto pit = pred_index.find(rid);
                pit != pred_index.end() &&
                !has_rna_suite(pit->second, residue_of(pit->second), typing))
                scored = false;
        }
        if (scored && sr.suiteness <= 0.0) {
            // L-RNA enantiomer check: a rejected suite may be mirror-image RNA
            // (e.g. a Spiegelmer aptamer). Negate every torsion (360 - t) and
            // re-classify; if the mirror fits a D-RNA cluster it is L-RNA  -
            // exclude it (D-RNA suiteness doesn't apply; the backbone is
            // genuinely high quality). A real D-RNA outlier won't fit even when
            // mirrored, so it still scores 0.
            std::array<double, 7> mirror;
            for (int k = 0; k < 7; ++k)
                mirror[k] = std::fmod(360.0 - (*suite)[k], 360.0);
            if (impl_->richardson->suiteness(mirror) > 0.0) scored = false;
        }
        // Full torsion row for every classifiable residue (incl. DNA, outliers,
        // L-RNA), with `scored` marking the rows that enter the backbone score,
        // so the emitted table covers every residue the parser could type.
        if (collect_backbone) {
            BackboneTorsion bt;
            bt.res_id = rid;
            if (const core::Residue* r = residue_of(rid)) bt.base_type = typing.normalize(r->base_type());
            bt.alpha = t.alpha; bt.beta = t.beta; bt.gamma = t.gamma;
            bt.delta = t.delta; bt.epsilon = t.epsilon; bt.zeta = t.zeta;
            bt.conformer = sr.conformer;
            bt.suiteness = sr.suiteness;
            bt.suite_bin = sr.bin;
            bt.is_outlier = sr.is_outlier;
            bt.scored = scored;
            ss.backbone_torsions.push_back(std::move(bt));
        }
        if (!scored) continue;
        ResidueScore rs;
        rs.res_id   = rid;
        rs.suiteness = sr.suiteness;
        rs.score    = std::max(0.0, std::min(100.0, sr.suiteness * 100.0));
        rs.recommendation = impl_->build_recommendation(*suite, sr);
        // Sugar pucker: delta, Pperp, and the report-only OUTLIER (delta/Pperp disagree).
        if (t.delta) {
            const double delta = *t.delta;
            rs.delta = delta;
            const int pd = impl_->richardson->pucker(delta);
            const core::Residue* residue = nullptr;
            if (const auto rit = res_by_id.find(rid); rit != res_by_id.end())
                residue = rit->second;
            const core::Residue* next_res = nullptr;
            if (const auto sit = succ_index.find(rid); sit != succ_index.end())
                if (const auto nit = res_by_id.find(sit->second); nit != res_by_id.end())
                    next_res = nit->second;
            const auto pp = compute_pperp(residue, next_res);
            if (pd == 0) {
                PuckerUnclassified pu;
                pu.delta = delta;
                if (pp) {
                    pu.pperp = *pp;
                    pu.pucker_by_pperp = pucker_label(pucker_from_pperp(*pp));
                }
                rs.pucker_unclassified = std::move(pu);
            }
            if (pp) {
                rs.pperp = pp;
                if (pd == 2 || pd == 3) {
                    const int pp_pucker = pucker_from_pperp(*pp);
                    if (pp_pucker != pd) {
                        PuckerOutlier po;
                        po.pperp = *pp;
                        po.delta = delta;
                        po.pucker_by_pperp = pucker_label(pp_pucker);
                        po.pucker_by_delta = pucker_label(pd);
                        rs.pucker_outlier = std::move(po);
                    }
                }
            }
        }
        ss.residue_scores.push_back(std::move(rs));
    }

    auto mean = [](const auto& v, auto getter) {
        if (v.empty()) return 0.0;
        double s = 0.0;
        for (const auto& x : v) s += getter(x);
        return s / v.size();
    };
    const double pm = mean(ss.pair_scores,    [](const PairScore& p)    { return p.score; });
    const double rm = mean(ss.residue_scores, [](const ResidueScore& r) { return r.score; });
    ss.pairs_score    = round2(pm);
    ss.residues_score = round2(rm);

    // Whole-structure RNA-RNA clashes (report-only; includes residues that are
    // not a base pair), read straight from the parsed chains.
    {
        // A map of pointers (the residues are not copied); which atom of a clash
        // is listed first follows its iteration order.
        std::unordered_map<std::string, const core::Residue*> by_id;
        by_id.reserve(res_by_id.size());
        for (const auto& [rid, res] : res_by_id) by_id.emplace(rid, res);
        ss.clashes = structure_clashes(by_id, typing, chem);
    }
    return ss;
}

}  // namespace pairfinder::scoring
