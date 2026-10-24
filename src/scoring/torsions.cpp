/**
 * @file torsions.cpp
 * @brief Backbone torsions.
 */
#include <pairfinder/scoring/torsions.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <tuple>
#include <vector>

#include <pairfinder/geometry/vector3d.hpp>

namespace pairfinder::scoring {

namespace {

using core::Residue;
using geometry::Vector3d;

constexpr double kPi = 3.14159265358979323846;

const Vector3d* pos(const Residue* r, const std::string& a) {
    if (r == nullptr) return nullptr;
    const core::Atom* at = r->get_atom(a);
    return at ? &at->coords : nullptr;
}

// IUPAC dihedral (degrees) for four points.
double dihedral(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3, const Vector3d& p4) {
    const Vector3d b1 = p2 - p1, b2 = p3 - p2, b3 = p4 - p3;
    Vector3d n1 = b1.cross(b2), n2 = b2.cross(b3);
    const double n1l = n1.norm(), n2l = n2.norm();
    if (n1l < 1e-10 || n2l < 1e-10) return 0.0;
    n1 = n1 / n1l;
    n2 = n2 / n2l;
    const Vector3d b2n = b2 / b2.norm();
    const Vector3d m1 = n1.cross(b2n);
    const double x = n1.dot(n2), y = m1.dot(n2);
    return std::atan2(-y, x) * 180.0 / kPi;
}

std::optional<double> dih(const Vector3d* a, const Vector3d* b, const Vector3d* c,
                          const Vector3d* d) {
    if (!a || !b || !c || !d) return std::nullopt;
    return dihedral(*a, *b, *c, *d);
}

ResidueTorsions compute_one(const Residue* res, const Residue* prev, const Residue* next) {
    const Vector3d* p = pos(res, "P");
    const Vector3d* o5 = pos(res, "O5'");
    const Vector3d* c5 = pos(res, "C5'");
    const Vector3d* c4 = pos(res, "C4'");
    const Vector3d* c3 = pos(res, "C3'");
    const Vector3d* o3 = pos(res, "O3'");
    const Vector3d* prev_o3 = pos(prev, "O3'");
    const Vector3d* next_p = pos(next, "P");
    const Vector3d* next_o5 = pos(next, "O5'");

    ResidueTorsions t;
    t.alpha = dih(prev_o3, p, o5, c5);
    t.beta = dih(p, o5, c5, c4);
    t.gamma = dih(o5, c5, c4, c3);
    t.delta = dih(c5, c4, c3, o3);
    t.epsilon = dih(c4, c3, o3, next_p);
    t.zeta = dih(c3, o3, next_p, next_o5);
    return t;
}

// Sort key: (chain, seq, icode) via
// parse_res_seq, falling back to (chain-before-first-'-', 0, "") on failure.
std::tuple<std::string, long, std::string> sort_key(const std::string& res_id) {
    if (const auto p = parse_res_seq(res_id)) return {p->chain, p->num, p->icode};
    const auto d = res_id.find('-');
    return {d == std::string::npos ? res_id : res_id.substr(0, d), 0L, std::string()};
}

}  // namespace

std::optional<ResSeq> parse_res_seq(const std::string& res_id) {
    const auto p1 = res_id.find('-');
    if (p1 == std::string::npos) return std::nullopt;
    const auto p2 = res_id.find('-', p1 + 1);
    if (p2 == std::string::npos) return std::nullopt;
    const std::string seqnum = res_id.substr(p2 + 1);  // everything after the second '-'
    if (seqnum.empty()) return std::nullopt;
    std::size_t i = (seqnum[0] == '-') ? 1 : 0;
    const std::size_t digit_start = i;
    while (i < seqnum.size() && std::isdigit(static_cast<unsigned char>(seqnum[i]))) ++i;
    if (i == digit_start) return std::nullopt;  // no digits
    const std::size_t num_end = i;
    for (std::size_t j = i; j < seqnum.size(); ++j)
        if (!std::isalpha(static_cast<unsigned char>(seqnum[j]))) return std::nullopt;
    try {
        return ResSeq{res_id.substr(0, p1), std::stol(seqnum.substr(0, num_end)),
                      seqnum.substr(num_end)};
    } catch (...) {
        return std::nullopt;
    }
}

std::unordered_map<std::string, ResidueTorsions> compute_all_torsions(
    const core::Structure& structure) {
    std::unordered_map<std::string, const Residue*> by_id;
    std::vector<std::string> ids;
    for (const auto& chain : structure.chains())
        for (const auto& res : chain.residues()) {
            by_id[res.res_id()] = &res;
            ids.push_back(res.res_id());
        }
    // Sort by (chain, seq, icode); group by chain; prev/next = adjacent. Each
    // id's key is parsed once (not on every comparison); a stable sort on the
    // precomputed keys gives the same order.
    std::vector<std::tuple<std::string, long, std::string>> keys;
    keys.reserve(ids.size());
    for (const auto& id : ids) keys.push_back(sort_key(id));
    std::vector<std::size_t> order(ids.size());
    for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
    std::stable_sort(order.begin(), order.end(),
                     [&keys](std::size_t a, std::size_t b) { return keys[a] < keys[b]; });
    std::vector<std::string> sorted_ids;
    std::vector<const std::string*> chain_of;   // each id's chain, std::get<0>(key)
    sorted_ids.reserve(ids.size());
    chain_of.reserve(ids.size());
    for (const std::size_t k : order) {
        sorted_ids.push_back(std::move(ids[k]));
        chain_of.push_back(&std::get<0>(keys[k]));
    }
    ids = std::move(sorted_ids);

    std::unordered_map<std::string, ResidueTorsions> out;
    std::size_t i = 0;
    while (i < ids.size()) {
        const std::string& chain = *chain_of[i];
        std::size_t j = i;
        while (j < ids.size() && *chain_of[j] == chain) ++j;
        for (std::size_t k = i; k < j; ++k) {
            const Residue* prev = (k > i) ? by_id[ids[k - 1]] : nullptr;
            const Residue* next = (k + 1 < j) ? by_id[ids[k + 1]] : nullptr;
            out[ids[k]] = compute_one(by_id[ids[k]], prev, next);
        }
        i = j;
    }
    return out;
}

}  // namespace pairfinder::scoring
