/**
 * @file pdb_parser.cpp
 * @brief GEMMI-backed structure loader (PDB / mmCIF / .gz).
 *
 * Parser conventions:
 *   - first model only
 *   - chain = gemmi chain name; seq = auth seqid.num; insertion code appended
 *   - alt-loc filter: keep ' ' / 'A' / '1' only (drops B-only atoms)
 *   - HETATM kept only for modified nucleotides (registry); waters/ligands dropped
 *   - atom-name normalization: OP1->O1P, OP2->O2P, OP3->O3P, O1'->O4', OL->O1P,
 *     OR->O2P, C5A->C5M, O5T->O5', O3T->O3', and '*'->'\'' (legacy spelling)
 *   - res_id = "<chain>-<resname>-<seq+icode>", base_type = upper(resname)
 * GEMMI auto-detects format by extension, so .pdb, .cif, .pdb.gz, .cif.gz work.
 */
#include <pairfinder/io/pdb_parser.hpp>

#include <cctype>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <gemmi/gz.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/model.hpp>

#include <pairfinder/core/base_typing.hpp>
#include <pairfinder/core/nucleotide_registry.hpp>
#include <pairfinder/core/resource_locator.hpp>
#include <pairfinder/core/residue.hpp>

namespace pairfinder::io {

namespace {

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// The ligand H-bond DB's nucleotide typing, loaded once. The registry below
// covers frame templates; this covers everything the rest of the pipeline can
// type, and the two disagree on ~64 codes.
const core::BaseTyping& base_typing() {
    static const core::BaseTyping typing(resources::config("ligand_hbond_db.json").string());
    return typing;
}

// modified_nucleotides.json: env PARSE_NT_CONFIG, else the bundled resource.
// Loaded once.
const core::NucleotideRegistry& nucleotide_registry() {
    static const core::NucleotideRegistry registry([] {
        if (const char* env = std::getenv("PARSE_NT_CONFIG")) return std::string(env);
        return resources::config("modified_nucleotides.json").string();
    }());
    return registry;
}

bool is_water(const std::string& name) {
    static const std::unordered_set<std::string> kWaters = {"HOH", "WAT", "H2O", "OH2", "SOL"};
    return kWaters.count(upper(name)) != 0;
}

// keep alt-loc ' ', '\0', 'A', '1'
bool keep_altloc(char a) { return a == ' ' || a == '\0' || a == 'A' || a == '1'; }

// Atom-name exact matches, applied to the trimmed GEMMI name after
// '*'->'\'' replacement. Returns the normalized (trimmed) atom name.
std::string normalize_atom(std::string name) {
    for (char& c : name)
        if (c == '*') c = '\'';
    static const std::unordered_map<std::string, std::string> kMap = {
        {"OP1", "O1P"}, {"OP2", "O2P"}, {"OP3", "O3P"}, {"OL", "O1P"}, {"OR", "O2P"},
        {"O1'", "O4'"}, {"C5A", "C5M"}, {"O5T", "O5'"}, {"O3T", "O3'"}};
    const auto it = kMap.find(name);
    return it == kMap.end() ? name : it->second;
}

}  // namespace

core::Structure parse_pdb(const std::filesystem::path& path) {
    gemmi::Structure gs;
    try {
        gs = gemmi::read_structure(gemmi::MaybeGzipped(path.string()));
    } catch (const std::exception& e) {
        throw std::runtime_error("cannot parse structure file " + path.string() + ": " + e.what());
    }

    core::Structure structure(path.stem().string());
    if (gs.models.empty()) return structure;
    const gemmi::Model& model = gs.models.front();

    std::unordered_map<std::string, std::pair<std::size_t, std::size_t>> res_loc;
    std::unordered_map<std::string, std::size_t> chain_idx;

    for (const gemmi::Chain& chain : model.chains) {
        const std::string chain_id = chain.name.empty() ? "A" : chain.name;
        for (const gemmi::Residue& residue : chain.residues) {
            const std::string res_name = residue.name;
            const bool is_hetatm = (residue.het_flag == 'H');
            // HETATM kept only for modified nucleotides (standard ATOM always kept);
            // waters/ligands dropped (HETATM included only for modified bases).
            if (is_hetatm) {
                if (is_water(res_name)) continue;
                // Keep if EITHER registry knows it. modified_nucleotides.json
                // drives frame-template selection and misses ~64 codes the
                // ligand DB types as RNA/DNA (U37, U8U, M3X, SRA, ...). Dropping
                // those does not just lose one residue: it breaks the chain, so
                // the neighbours lose their predecessor and no suite can be
                // built.
                if (!nucleotide_registry().lookup(res_name).has_value() &&
                    !base_typing().is_known_nucleotide(res_name))
                    continue;
            }

            const long seq = residue.seqid.num.value;
            std::string ins;
            const char icode = residue.seqid.icode;
            if (icode != ' ' && icode != '\0') ins = std::string(1, icode);
            const std::string seq_with_ins = std::to_string(seq) + ins;
            const std::string res_id = core::make_res_id(chain_id, res_name, seq_with_ins);
            const std::string base_type = upper(res_name);

            // The residue is resolved at its first kept atom (a residue whose atoms
            // are all filtered out is never created), then reused for the rest:
            // one lookup per residue rather than one per atom.
            core::Residue* target = nullptr;
            for (const gemmi::Atom& atom : residue.atoms) {
                if (!keep_altloc(atom.altloc)) continue;
                const std::string atom_name = normalize_atom(atom.name);
                if (target != nullptr) {
                    target->add_atom(atom_name, {atom.pos.x, atom.pos.y, atom.pos.z}, "",
                                     atom.altloc);
                    continue;
                }

                auto loc_it = res_loc.find(res_id);
                if (loc_it == res_loc.end()) {
                    auto ci_it = chain_idx.find(chain_id);
                    std::size_t ci;
                    if (ci_it == chain_idx.end()) {
                        ci = structure.chains().size();
                        structure.add_chain(core::Chain(chain_id));
                        chain_idx.emplace(chain_id, ci);
                    } else {
                        ci = ci_it->second;
                    }
                    core::Chain& c = structure.chains()[ci];
                    const std::size_t ri = c.residues().size();
                    c.add_residue(core::Residue(res_id, base_type));
                    loc_it = res_loc.emplace(res_id, std::make_pair(ci, ri)).first;
                }
                auto [ci, ri] = loc_it->second;
                target = &structure.chains()[ci].residues()[ri];
                target->add_atom(atom_name, {atom.pos.x, atom.pos.y, atom.pos.z}, "", atom.altloc);
            }
        }
    }

    return structure;
}

}  // namespace pairfinder::io
