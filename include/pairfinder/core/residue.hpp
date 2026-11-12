/**
 * @file residue.hpp
 * @brief A nucleotide residue: identity + its atoms.
 *
 * Fields:
 *   - constructed from (res_id, base_type), res_id == "<chain>-<resname>-<seq>"
 *   - chain()/seq_num() derived by splitting res_id
 *   - base_type is the parse-time name (upper/stripped; DNA & modified codes are
 *     KEPT, never converted) - conversion to a parent base happens only in
 *     is_purine/is_pyrimidine/normalize_base_type below.
 *   - add_atom overwrites an existing atom of the same name.
 */
#ifndef PAIRFINDER_CORE_RESIDUE_HPP
#define PAIRFINDER_CORE_RESIDUE_HPP

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <pairfinder/core/atom.hpp>

namespace pairfinder::core {

/// True for purine base types (A, G, I), normalizing DNA/modified codes.
bool is_purine(std::string_view base_type);
/// True for pyrimidine base types (C, U, T), normalizing DNA/modified codes.
bool is_pyrimidine(std::string_view base_type);
/// Classification normalizer: standard, DNA prefix (DA->A), modified->parent,
/// else the name unchanged. (Distinct from the parse-time keep-original rule.)
std::string normalize_base_type(std::string_view residue_name);
/// Build a residue id "<chain>-<resname>-<seq>".
std::string make_res_id(std::string_view chain, std::string_view res_name,
                        std::string_view seq);

/// Pipeline base type from a res_id: strips a
/// leading 'D' from 2-char DNA codes (DA->A, DC->C), keeps everything else.
std::string pipeline_base_type(std::string_view res_id);

class Residue;  // fwd

/// Convert a parsed residue to its pipeline form: pipeline base type + phosphate
/// atom renaming OP1->O1P, OP2->O2P, OP3->O3P. Apply AFTER frame computation.
/// Frames/candidates are unaffected; modified-base phosphate H-bond capacity
/// lookups then find the O-form names.
void apply_pipeline_view(Residue& residue);

class Residue {
public:
    Residue() = default;
    Residue(std::string res_id, std::string base_type)
        : res_id_(std::move(res_id)), base_type_(std::move(base_type)) {}

    /// Overwrites any existing atom with the same name.
    void add_atom(const std::string& name, const geometry::Vector3d& coords,
                  const std::string& element = "", char alt_loc = ' ');

    /// Rename an atom in place (no-op if absent; type recomputed from new name).
    void rename_atom(const std::string& old_name, const std::string& new_name);

    const std::string& res_id() const { return res_id_; }
    /// Chain id = first '-'-delimited field of res_id.
    std::string chain() const;
    /// Sequence (+ insertion code) = last '-'-delimited field of res_id.
    std::string seq_num() const;
    /// Base name as parsed (upper/stripped; modified/DNA kept).
    const std::string& base_type() const { return base_type_; }
    /// Override the base type (used to switch a parsed residue to its pipeline
    /// base type).
    void set_base_type(std::string base_type) { base_type_ = std::move(base_type); }

    const std::vector<Atom>& atoms() const { return atoms_; }

    bool has_atom(std::string_view name) const {
        return index_.count(std::string(name)) > 0;
    }
    const Atom* get_atom(std::string_view name) const;

    bool has_atom_type(AtomType type) const;
    const Atom* get_atom_type(AtomType type) const;

    /// Glycosidic nitrogen (N9 purine / N1 pyrimidine), or nullptr.
    const Atom* glycosidic_n() const;

private:
    std::string res_id_;
    std::string base_type_;
    std::vector<Atom> atoms_;
    std::unordered_map<std::string, std::size_t> index_;  ///< name -> atoms_ idx
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_RESIDUE_HPP
