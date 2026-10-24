/**
 * @file richardson.cpp
 * @brief Richardson suite classifier.
 */
#include <pairfinder/scoring/richardson.hpp>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>

namespace pairfinder::scoring {

namespace {

constexpr double kPi = 3.14159265358979323846;
using Widths = std::array<double, 7>;

double angle_diff(double a, double b) {
    const double d = std::abs(a - b);
    return std::min(d, 360.0 - d);
}

template <typename Range>
double l3_distance(const std::array<double, 7>& suite, const std::array<double, 7>& center,
                   const Widths& widths, Range begin, Range end) {
    double total = 0.0;
    for (Range k = begin; k < end; ++k) {
        const double d = angle_diff(suite[k], center[k]);
        const double r = d / widths[k];
        total += r * r * r;
    }
    return std::cbrt(total);
}

struct Conformer {
    std::string name;
    std::string bin;
    std::array<double, 7> angles{};
    std::string dominance;
};

}  // namespace

struct RichardsonClassifier::Impl {
    std::vector<Conformer> conformers;
    std::unordered_map<std::string, std::vector<int>> by_bin;
    std::unordered_map<std::string, int> by_name;
    Widths normal{}, satellite{};
    std::unordered_map<std::string, Widths> sat_widths_ov;
    std::unordered_map<std::string, Widths> dom_widths_ov;
    std::array<double, 2> tri_eps{}, tri_alpha{}, tri_beta{}, tri_zeta{};
    std::array<double, 2> d_c3{}, d_c2{};
    std::array<double, 2> g_p{}, g_t{}, g_m{};

    int pucker(double delta) const {
        if (d_c3[0] <= delta && delta <= d_c3[1]) return 3;
        if (d_c2[0] <= delta && delta <= d_c2[1]) return 2;
        return 0;
    }
    std::string gamma(double g) const {
        if (g_p[0] <= g && g <= g_p[1]) return "p";
        if (g_t[0] <= g && g <= g_t[1]) return "t";
        if (g_m[0] <= g && g <= g_m[1]) return "m";
        return "";
    }
    Widths compose(const Widths& ov, const Widths& base) const {
        Widths r = base;
        for (int k = 0; k < 7; ++k)
            if (ov[k] > 0) r[k] = ov[k];
        return r;
    }
    Widths widths_for_satellite(const std::string& name) const {
        const auto it = sat_widths_ov.find(name);
        return it == sat_widths_ov.end() ? satellite : compose(it->second, satellite);
    }
    Widths widths_for_dominant(const std::array<double, 7>& suite,
                               const std::array<double, 7>& dom_center,
                               const std::vector<int>& bin_confs) const {
        double best_d = std::numeric_limits<double>::infinity();
        Widths best = normal;
        bool found = false;
        for (int ci : bin_confs) {
            if (conformers[ci].dominance != "sat") continue;
            const auto it = dom_widths_ov.find(conformers[ci].name);
            if (it == dom_widths_ov.end()) continue;
            const Widths w = compose(it->second, normal);
            const double d = l3_distance(suite, dom_center, w, 1, 5);
            if (d < best_d) { best_d = d; best = w; found = true; }
        }
        return found ? best : normal;
    }
    Widths resolve_widths(const Conformer& conf, const std::vector<int>& bin_confs,
                          const std::array<double, 7>& suite) const {
        if (conf.dominance == "sat") return widths_for_satellite(conf.name);
        if (conf.dominance == "dom") return widths_for_dominant(suite, conf.angles, bin_confs);
        return normal;
    }
};

namespace {
Widths to_widths(const nlohmann::json& a) {
    Widths w{};
    for (int k = 0; k < 7 && k < static_cast<int>(a.size()); ++k) w[k] = a[k].get<double>();
    return w;
}
std::array<double, 2> to_range(const nlohmann::json& a) {
    return {a[0].get<double>(), a[1].get<double>()};
}
}  // namespace

RichardsonClassifier::RichardsonClassifier(const std::filesystem::path& path)
    : impl_(std::make_unique<Impl>()) {
    std::ifstream in(path);
    const nlohmann::json d = nlohmann::json::parse(in);
    for (const auto& c : d["conformers"]) {
        Conformer cf;
        cf.name = c["name"].get<std::string>();
        cf.bin = c["bin"].get<std::string>();
        const auto& ang = c["angles"];
        for (int k = 0; k < 7; ++k) cf.angles[k] = ang[k].get<double>();
        cf.dominance = c.value("dominance", std::string("ord"));
        const int idx = static_cast<int>(impl_->conformers.size());
        impl_->by_bin[cf.bin].push_back(idx);
        impl_->by_name[cf.name] = idx;
        impl_->conformers.push_back(std::move(cf));
    }
    impl_->normal = to_widths(d["widths"]["normal"]);
    impl_->satellite = d["widths"].contains("satellite") ? to_widths(d["widths"]["satellite"])
                                                         : impl_->normal;
    if (d.contains("satellite_overrides")) {
        // NB: iterate the lvalue, not a temporary - .items() on a temporary json
        // dangles and silently yields nothing.
        const nlohmann::json& overrides = d["satellite_overrides"];
        for (const auto& [name, ov] : overrides.items()) {
            if (ov.contains("sat_widths")) impl_->sat_widths_ov[name] = to_widths(ov["sat_widths"]);
            if (ov.contains("dom_widths")) impl_->dom_widths_ov[name] = to_widths(ov["dom_widths"]);
        }
    }
    impl_->tri_eps = to_range(d["triage"]["epsilon"]);
    impl_->tri_alpha = to_range(d["triage"]["alpha"]);
    impl_->tri_beta = to_range(d["triage"]["beta"]);
    impl_->tri_zeta = to_range(d["triage"]["zeta"]);
    impl_->d_c3 = to_range(d["sieve"]["delta"]["C3_endo"]);
    impl_->d_c2 = to_range(d["sieve"]["delta"]["C2_endo"]);
    impl_->g_p = to_range(d["sieve"]["gamma"]["p"]);
    impl_->g_t = to_range(d["sieve"]["gamma"]["t"]);
    impl_->g_m = to_range(d["sieve"]["gamma"]["m"]);
}

RichardsonClassifier::~RichardsonClassifier() = default;
RichardsonClassifier::RichardsonClassifier(RichardsonClassifier&&) noexcept = default;

double signed_angle_gap(double value, double center) {
    double d = std::fmod(value - center + 180.0, 360.0);
    if (d < 0.0) d += 360.0;
    return d - 180.0;
}

SuiteResult RichardsonClassifier::classify(const std::array<double, 7>& suite) const {
    const Impl& m = *impl_;
    SuiteResult out;  // defaults: conformer "" -> set to "!!" on any rejection
    out.conformer = "!!";
    out.is_outlier = true;

    // Triage.
    auto chk = [](double v, const std::array<double, 2>& r) { return r[0] <= v && v <= r[1]; };
    if (!chk(suite[1], m.tri_eps) || !chk(suite[3], m.tri_alpha) ||
        !chk(suite[4], m.tri_beta) || !chk(suite[2], m.tri_zeta))
        return out;
    // Sieve.
    const int pdm = m.pucker(suite[0]);
    const int pd = m.pucker(suite[6]);
    const std::string g = m.gamma(suite[5]);
    if (pdm == 0 || pd == 0 || g.empty()) return out;
    const std::string bin = std::to_string(pdm) + std::to_string(pd) + g;
    out.bin = bin;  // recorded once triage+sieve pass, even for outliers below
    const auto it = m.by_bin.find(bin);
    if (it == m.by_bin.end() || it->second.empty()) return out;

    // 4D screen.
    double best_d4 = std::numeric_limits<double>::infinity();
    int best = -1;
    Widths best_w = m.normal;
    for (int ci : it->second) {
        const Widths w = m.resolve_widths(m.conformers[ci], it->second, suite);
        const double d4 = l3_distance(suite, m.conformers[ci].angles, w, 1, 5);
        if (d4 < best_d4) { best_d4 = d4; best = ci; best_w = w; }
    }
    if (best < 0 || best_d4 >= 1.0) return out;

    // 7D test.
    const double d7 = l3_distance(suite, m.conformers[best].angles, best_w, 0, 7);
    if (d7 > 1.0) return out;

    double s = (std::cos(kPi * d7) + 1.0) / 2.0;
    if (s < 0.01) s = 0.01;
    out.conformer = m.conformers[best].name;
    out.suiteness = std::round(s * 1000.0) / 1000.0;
    out.is_outlier = false;
    out.distance = d7;
    return out;
}

double RichardsonClassifier::suiteness(const std::array<double, 7>& suite) const {
    return classify(suite).suiteness;
}

const std::array<double, 7>* RichardsonClassifier::center_for(
    const std::string& conformer) const {
    const auto it = impl_->by_name.find(conformer);
    return it == impl_->by_name.end() ? nullptr : &impl_->conformers[it->second].angles;
}

std::string RichardsonClassifier::nearest_conformer(
    const std::array<double, 7>& suite) const {
    std::string best_name;
    double best_d = std::numeric_limits<double>::infinity();
    for (const auto& c : impl_->conformers) {
        double d = 0.0;
        for (int k = 0; k < 7; ++k) d += angle_diff(suite[k], c.angles[k]);
        if (d < best_d) { best_d = d; best_name = c.name; }
    }
    return best_name;
}

const std::array<double, 7>& RichardsonClassifier::normal_widths() const {
    return impl_->normal;
}

int RichardsonClassifier::pucker(double delta) const {
    return impl_->pucker(std::fmod(delta, 360.0));
}

}  // namespace pairfinder::scoring
