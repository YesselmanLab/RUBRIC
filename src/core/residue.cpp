/**
 * @file residue.cpp
 * @brief Residue methods + base-type helpers.
 *
 */
#include <pairfinder/core/residue.hpp>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>

namespace pairfinder::core {

namespace {

const std::unordered_set<std::string> kPurines = {"A", "G", "I"};
const std::unordered_set<std::string> kPyrimidines = {"C", "U", "T"};

// DNA prefix + a small set of common modified codes. The full data-driven map
// lives in the bundled base-typing cache; this is the compiled-in fallback.
const std::unordered_map<std::string, std::string>& base_map() {
    static const std::unordered_map<std::string, std::string> m = {
        {"DA", "A"}, {"DG", "G"}, {"DC", "C"}, {"DT", "T"}, {"DU", "U"},
        {"PSU", "U"}, {"5MU", "U"}, {"H2U", "U"}, {"4SU", "U"},
        {"5MC", "C"}, {"OMC", "C"},
        {"1MA", "A"}, {"2MA", "A"}, {"6MA", "A"},
        {"2MG", "G"}, {"7MG", "G"}, {"OMG", "G"}, {"1MG", "G"},
    };
    return m;
}

}  // namespace

std::string normalize_base_type(std::string_view residue_name) {
    const std::string key(residue_name);
    if (kPurines.count(key) || kPyrimidines.count(key)) return key;
    const auto it = base_map().find(key);
    if (it != base_map().end()) return it->second;
    return key;  // unknown modified base: keep original
}

bool is_purine(std::string_view base_type) {
    return kPurines.count(normalize_base_type(base_type)) > 0;
}

bool is_pyrimidine(std::string_view base_type) {
    return kPyrimidines.count(normalize_base_type(base_type)) > 0;
}

std::string pipeline_base_type(std::string_view res_id) {
    const auto p1 = res_id.find('-');
    if (p1 == std::string_view::npos) return std::string(res_id);
    const auto p2 = res_id.find('-', p1 + 1);
    const std::string_view name = res_id.substr(
        p1 + 1, (p2 == std::string_view::npos ? res_id.size() : p2) - p1 - 1);
    if (name.size() == 2 && name[0] == 'D') return std::string(name.substr(1));
    return std::string(name);
}

void apply_pipeline_view(Residue& residue) {
    residue.set_base_type(pipeline_base_type(residue.res_id()));
    residue.rename_atom("OP1", "O1P");
    residue.rename_atom("OP2", "O2P");
    residue.rename_atom("OP3", "O3P");
}

std::string make_res_id(std::string_view chain, std::string_view res_name,
                        std::string_view seq) {
    std::string id;
    id.reserve(chain.size() + res_name.size() + seq.size() + 2);
    id.append(chain).push_back('-');
    id.append(res_name).push_back('-');
    id.append(seq);
    return id;
}

std::string Residue::chain() const {
    const auto pos = res_id_.find('-');
    return pos == std::string::npos ? res_id_ : res_id_.substr(0, pos);
}

std::string Residue::seq_num() const {
    const auto pos = res_id_.rfind('-');
    return pos == std::string::npos ? std::string{} : res_id_.substr(pos + 1);
}

void Residue::add_atom(const std::string& name, const geometry::Vector3d& coords,
                       const std::string& element, char alt_loc) {
    const auto it = index_.find(name);
    if (it != index_.end()) {
        atoms_[it->second] = Atom(name, coords, element, alt_loc);  // overwrite
        return;
    }
    index_[name] = atoms_.size();
    atoms_.emplace_back(name, coords, element, alt_loc);
}

void Residue::rename_atom(const std::string& old_name, const std::string& new_name) {
    const auto it = index_.find(old_name);
    if (it == index_.end()) return;
    const std::size_t i = it->second;
    atoms_[i] = Atom(new_name, atoms_[i].coords, atoms_[i].element, atoms_[i].alt_loc);
    index_.erase(it);
    index_[new_name] = i;
}

const Atom* Residue::get_atom(std::string_view name) const {
    const auto it = index_.find(std::string(name));
    return it == index_.end() ? nullptr : &atoms_[it->second];
}

bool Residue::has_atom_type(AtomType type) const {
    return std::any_of(atoms_.begin(), atoms_.end(),
                       [type](const Atom& a) { return a.type == type; });
}

const Atom* Residue::get_atom_type(AtomType type) const {
    for (const auto& a : atoms_) {
        if (a.type == type) return &a;
    }
    return nullptr;
}

const Atom* Residue::glycosidic_n() const {
    if (is_purine(base_type_)) return get_atom_type(AtomType::N9);
    if (is_pyrimidine(base_type_)) return get_atom_type(AtomType::N1);
    return nullptr;
}

}  // namespace pairfinder::core
