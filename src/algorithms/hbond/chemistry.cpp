/**
 * @file chemistry.cpp
 * @brief H-bond capacities, connectivity, and slot prediction.
 */
#include <pairfinder/algorithms/hbond/chemistry.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>

#include <nlohmann/json.hpp>

namespace pairfinder::algorithms::hbond {

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string upper(std::string s) {
    for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

}  // namespace

Vector3d hb_normalize(const Vector3d& v) {
    const double n = v.norm();
    if (n < 1e-10) return v;
    return v / n;
}

double angle_between(const Vector3d& a, const Vector3d& b) {
    const Vector3d ua = hb_normalize(a);
    const Vector3d ub = hb_normalize(b);
    double dot = ua.dot(ub);
    if (dot > 1.0) dot = 1.0;
    if (dot < -1.0) dot = -1.0;
    return std::acos(dot) * 180.0 / kPi;
}

Vector3d rotate_vector(const Vector3d& v, const Vector3d& axis_in, double angle_deg) {
    const double a = angle_deg * kPi / 180.0;
    const Vector3d axis = hb_normalize(axis_in);
    const double cos_a = std::cos(a);
    const double sin_a = std::sin(a);
    const Vector3d term1 = v * cos_a;
    const Vector3d term2 = axis.cross(v) * sin_a;
    const Vector3d term3 = axis * (axis.dot(v) * (1.0 - cos_a));
    return term1 + term2 + term3;
}

Vector3d compute_base_normal(const AtomMap& atoms) {
    static const char* kRing[] = {"C2", "C4", "C6", "N1", "N3"};
    std::vector<Vector3d> pos;
    for (const char* name : kRing) {
        const auto it = atoms.find(name);
        if (it != atoms.end()) pos.push_back(it->second);
    }
    if (pos.size() < 3) return {0.0, 0.0, 1.0};
    const Vector3d v1 = pos[1] - pos[0];
    const Vector3d v2 = pos[2] - pos[0];
    return hb_normalize(v1.cross(v2));
}

std::string HBondChemistry::key(const std::string& base, const std::string& atom) {
    return upper(base) + '\x01' + atom;
}

namespace {
/// key(base, atom) written into a reused buffer. The lookups below run for every
/// atom of every candidate pair, and building a fresh string each time (copy,
/// upper-case, concatenate) was most of their cost; the key itself is unchanged.
const std::string& key_into(std::string& buf, const std::string& base, const std::string& atom) {
    buf.clear();
    for (const char c : base) buf += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    buf += '\x01';
    buf += atom;
    return buf;
}
}  // namespace

const std::vector<std::string>&
HBondChemistry::connectivity(const std::string& base, const std::string& atom) const {
    thread_local std::string buf;
    const auto it = connectivity_.find(key_into(buf, base, atom));
    return it == connectivity_.end() ? empty_ : it->second;
}

std::optional<int> HBondChemistry::donor_capacity(const std::string& base,
                                                  const std::string& atom) const {
    thread_local std::string buf;
    const auto it = donor_cap_.find(key_into(buf, base, atom));
    if (it == donor_cap_.end()) return std::nullopt;
    return it->second;
}

std::optional<int> HBondChemistry::acceptor_capacity(const std::string& base,
                                                     const std::string& atom) const {
    thread_local std::string buf;
    const auto it = acceptor_cap_.find(key_into(buf, base, atom));
    if (it == acceptor_cap_.end()) return std::nullopt;
    return it->second;
}

// ---------------------------------------------------------------------------
//  Slot prediction
// ---------------------------------------------------------------------------

std::vector<HbSlot> HBondChemistry::predict_h_slots(const std::string& base,
                                                    const std::string& atom,
                                                    const AtomMap& atoms,
                                                    const Vector3d& base_normal) const {
    const auto cap = donor_capacity(base, atom);
    if (!cap) return {};
    const auto a_it = atoms.find(atom);
    if (a_it == atoms.end()) return {};
    const Vector3d donor_pos = a_it->second;

    std::vector<Vector3d> ant;
    for (const std::string& name : connectivity(base, atom)) {
        const auto it = atoms.find(name);
        if (it != atoms.end()) ant.push_back(it->second);
    }

    // Fallback isotropic slots for modified nucleotides without connectivity.
    if (ant.empty()) {
        std::vector<HbSlot> slots;
        slots.push_back({{1.0, 0.0, 0.0}, 1, {}});
        if (*cap >= 2) slots.push_back({{0.0, 1.0, 0.0}, 1, {}});
        if (*cap >= 3) slots.push_back({{0.0, 0.0, 1.0}, 1, {}});
        return slots;
    }

    std::vector<HbSlot> slots;
    if (*cap == 2 && ant.size() == 1) {  // sp2 NH2: two H at +/-120 deg
        const Vector3d ant_to_donor = hb_normalize(donor_pos - ant[0]);
        slots.push_back({rotate_vector(ant_to_donor, base_normal, 60.0), 1, {}});
        slots.push_back({rotate_vector(ant_to_donor, base_normal, -60.0), 1, {}});
    } else if (*cap == 1 && ant.size() == 2) {  // sp2 imino NH
        const Vector3d avg = (ant[0] + ant[1]) / 2.0;
        slots.push_back({hb_normalize(donor_pos - avg), 1, {}});
    } else if (*cap == 1 && ant.size() == 1) {  // H opposite single antecedent
        slots.push_back({hb_normalize(donor_pos - ant[0]), 1, {}});
    }
    return slots;
}

std::vector<HbSlot> HBondChemistry::predict_lp_slots(const std::string& base,
                                                     const std::string& atom,
                                                     const AtomMap& atoms,
                                                     const Vector3d& base_normal) const {
    const auto cap = acceptor_capacity(base, atom);
    if (!cap) return {};
    const auto a_it = atoms.find(atom);
    if (a_it == atoms.end()) return {};
    const Vector3d acceptor_pos = a_it->second;

    // Phosphate oxygens: three isotropic LP slots.
    if (atom == "OP1" || atom == "OP2" || atom == "O1P" || atom == "O2P") {
        return {{{1.0, 0.0, 0.0}, 2, {}},
                {{0.0, 1.0, 0.0}, 2, {}},
                {{0.0, 0.0, 1.0}, 2, {}}};
    }
    // Ribose oxygens: two perpendicular LP slots.
    if (atom == "O2'" || atom == "O4'" || atom == "O3'" || atom == "O5'") {
        Vector3d perp1 = base_normal.cross({1.0, 0.0, 0.0});
        if (perp1.norm() < 0.1) perp1 = base_normal.cross({0.0, 1.0, 0.0});
        perp1 = hb_normalize(perp1);
        const Vector3d perp2 = base_normal.cross(perp1);
        return {{perp1, 2, {}}, {perp2, 2, {}}};
    }

    std::vector<Vector3d> ant;
    for (const std::string& name : connectivity(base, atom)) {
        const auto it = atoms.find(name);
        if (it != atoms.end()) ant.push_back(it->second);
    }
    if (ant.empty()) {  // fallback isotropic
        std::vector<HbSlot> slots;
        slots.push_back({{1.0, 0.0, 0.0}, 2, {}});
        if (*cap >= 2) slots.push_back({{0.0, 1.0, 0.0}, 2, {}});
        return slots;
    }

    std::vector<HbSlot> slots;
    if (*cap == 2 && ant.size() == 1) {  // sp2 carbonyl: two LP at +/-120 deg
        const Vector3d ant_to_acc = hb_normalize(acceptor_pos - ant[0]);
        slots.push_back({rotate_vector(ant_to_acc, base_normal, 60.0), 1, {}});
        slots.push_back({rotate_vector(ant_to_acc, base_normal, -60.0), 1, {}});
    } else if (*cap == 1 && ant.size() == 2) {  // sp2 ring nitrogen
        const Vector3d avg = (ant[0] + ant[1]) / 2.0;
        slots.push_back({hb_normalize(acceptor_pos - avg), 1, {}});
    } else if (*cap == 1 && ant.size() == 1) {
        slots.push_back({hb_normalize(acceptor_pos - ant[0]), 1, {}});
    }
    return slots;
}

// ---------------------------------------------------------------------------
//  Table construction
// ---------------------------------------------------------------------------

namespace {

struct Conn { const char* base; const char* atom; std::vector<const char*> ant; };
struct Cap  { const char* base; const char* atom; int cap; };

}  // namespace

HBondChemistry::HBondChemistry(const std::filesystem::path& ligand_db_path, bool use_cache) {
    // --- Base connectivity ---
    static const Conn kConn[] = {
        {"A","N6",{"C6"}}, {"A","N1",{"C2","C6"}}, {"A","N3",{"C2","C4"}},
        {"A","N7",{"C5","C8"}}, {"A","C2",{"N1","N3"}}, {"A","C8",{"N7","N9"}},
        {"G","N1",{"C2","C6"}}, {"G","N2",{"C2"}}, {"G","O6",{"C6"}},
        {"G","N3",{"C2","C4"}}, {"G","N7",{"C5","C8"}}, {"G","C8",{"N7","N9"}},
        {"C","N4",{"C4"}}, {"C","N3",{"C2","C4"}}, {"C","O2",{"C2"}},
        {"C","C5",{"C4","C6"}}, {"C","C6",{"C5","N1"}},
        {"U","N3",{"C2","C4"}}, {"U","O2",{"C2"}}, {"U","O4",{"C4"}},
        {"U","C5",{"C4","C6"}}, {"U","C6",{"C5","N1"}},
        {"T","N3",{"C2","C4"}}, {"T","O2",{"C2"}}, {"T","O4",{"C4"}}, {"T","C6",{"C5","N1"}},
        {"A","O2'",{"C2'"}}, {"G","O2'",{"C2'"}}, {"C","O2'",{"C2'"}}, {"U","O2'",{"C2'"}}, {"T","O2'",{"C2'"}},
        {"A","O4'",{"C1'","C4'"}}, {"G","O4'",{"C1'","C4'"}}, {"C","O4'",{"C1'","C4'"}},
        {"U","O4'",{"C1'","C4'"}}, {"T","O4'",{"C1'","C4'"}},
        {"A","O3'",{"C3'"}}, {"G","O3'",{"C3'"}}, {"C","O3'",{"C3'"}}, {"U","O3'",{"C3'"}}, {"T","O3'",{"C3'"}},
        {"A","O5'",{"C5'"}}, {"G","O5'",{"C5'"}}, {"C","O5'",{"C5'"}}, {"U","O5'",{"C5'"}}, {"T","O5'",{"C5'"}},
        {"P","N1",{"C2","C6"}}, {"P","N3",{"C2","C4"}}, {"P","O2",{"C2"}}, {"P","O4",{"C4"}},
        {"P","O2'",{"C2'"}}, {"P","O4'",{"C1'","C4'"}}, {"P","O3'",{"C3'"}}, {"P","O5'",{"C5'"}},
        {"I","N1",{"C2","C6"}}, {"I","O6",{"C6"}}, {"I","N3",{"C2","C4"}}, {"I","N7",{"C5","C8"}},
        {"I","O2'",{"C2'"}}, {"I","O4'",{"C1'","C4'"}}, {"I","O3'",{"C3'"}}, {"I","O5'",{"C5'"}},
        {"DA","N6",{"C6"}}, {"DA","N1",{"C2","C6"}}, {"DA","N3",{"C2","C4"}}, {"DA","N7",{"C5","C8"}},
        {"DA","O4'",{"C1'","C4'"}}, {"DA","O3'",{"C3'"}}, {"DA","O5'",{"C5'"}},
        {"DG","N1",{"C2","C6"}}, {"DG","N2",{"C2"}}, {"DG","O6",{"C6"}}, {"DG","N3",{"C2","C4"}},
        {"DG","N7",{"C5","C8"}}, {"DG","O4'",{"C1'","C4'"}}, {"DG","O3'",{"C3'"}}, {"DG","O5'",{"C5'"}},
        {"DC","N4",{"C4"}}, {"DC","N3",{"C2","C4"}}, {"DC","O2",{"C2"}},
        {"DC","O4'",{"C1'","C4'"}}, {"DC","O3'",{"C3'"}}, {"DC","O5'",{"C5'"}},
        {"DT","N3",{"C2","C4"}}, {"DT","O2",{"C2"}}, {"DT","O4",{"C4"}},
        {"DT","O4'",{"C1'","C4'"}}, {"DT","O3'",{"C3'"}}, {"DT","O5'",{"C5'"}},
        // Amino-acid backbone + side-chain antecedents
        {"ALA","N",{"CA"}},{"ALA","O",{"C"}}, {"ARG","N",{"CA"}},{"ARG","O",{"C"}},
        {"ARG","NE",{"CD","CZ"}},{"ARG","NH1",{"CZ"}},{"ARG","NH2",{"CZ"}},
        {"ASN","N",{"CA"}},{"ASN","O",{"C"}},{"ASN","OD1",{"CG"}},{"ASN","ND2",{"CG"}},
        {"ASP","N",{"CA"}},{"ASP","O",{"C"}},{"ASP","OD1",{"CG"}},{"ASP","OD2",{"CG"}},
        {"CYS","N",{"CA"}},{"CYS","O",{"C"}},{"CYS","SG",{"CB"}},
        {"GLN","N",{"CA"}},{"GLN","O",{"C"}},{"GLN","OE1",{"CD"}},{"GLN","NE2",{"CD"}},
        {"GLU","N",{"CA"}},{"GLU","O",{"C"}},{"GLU","OE1",{"CD"}},{"GLU","OE2",{"CD"}},
        {"GLY","N",{"CA"}},{"GLY","O",{"C"}},
        {"HIS","N",{"CA"}},{"HIS","O",{"C"}},{"HIS","ND1",{"CG","CE1"}},{"HIS","NE2",{"CD2","CE1"}},
        {"ILE","N",{"CA"}},{"ILE","O",{"C"}}, {"LEU","N",{"CA"}},{"LEU","O",{"C"}},
        {"LYS","N",{"CA"}},{"LYS","O",{"C"}},{"LYS","NZ",{"CE"}},
        {"MET","N",{"CA"}},{"MET","O",{"C"}}, {"PHE","N",{"CA"}},{"PHE","O",{"C"}},
        {"PRO","N",{"CA","CD"}},{"PRO","O",{"C"}}, {"SER","N",{"CA"}},{"SER","O",{"C"}},{"SER","OG",{"CB"}},
        {"THR","N",{"CA"}},{"THR","O",{"C"}},{"THR","OG1",{"CB"}},
        {"TRP","N",{"CA"}},{"TRP","O",{"C"}},{"TRP","NE1",{"CD1","CE2"}},
        {"TYR","N",{"CA"}},{"TYR","O",{"C"}},{"TYR","OH",{"CZ"}}, {"VAL","N",{"CA"}},{"VAL","O",{"C"}},
    };
    for (const auto& c : kConn) {
        std::vector<std::string> ant(c.ant.begin(), c.ant.end());
        connectivity_[key(c.base, c.atom)] = std::move(ant);
    }

    // Fast path: load the precomputed capacity tables (an exact serialization of
    // the built-in + DB-extended tables below) instead of parsing the 1.5 MB
    // ligand DB. connectivity_ is hardcoded above and unaffected by the DB.
    if (use_cache && load_cache(ligand_db_path.parent_path() / "hbond_chemistry_cache.json"))
        return;

    // --- Donor capacity ---
    static const Cap kDonor[] = {
        {"A","N6",2},{"C","N4",2},{"G","N2",2},
        {"G","N1",1},{"U","N3",1},{"T","N3",1},
        {"A","O2'",1},{"G","O2'",1},{"C","O2'",1},{"U","O2'",1},{"T","O2'",1},{"P","O2'",1},{"I","O2'",1},
        {"P","N1",1},{"P","N3",1},{"I","N1",1},
        {"DA","N6",2},{"DG","N1",1},{"DG","N2",2},{"DC","N4",2},{"DT","N3",1},
        {"ALA","N",1},{"ARG","N",1},{"ASN","N",1},{"ASP","N",1},{"CYS","N",1},{"GLN","N",1},
        {"GLU","N",1},{"GLY","N",1},{"HIS","N",1},{"ILE","N",1},{"LEU","N",1},{"LYS","N",1},
        {"MET","N",1},{"PHE","N",1},{"SER","N",1},{"THR","N",1},{"TRP","N",1},{"TYR","N",1},{"VAL","N",1},
        {"ARG","NE",1},{"ARG","NH1",2},{"ARG","NH2",2},{"ASN","ND2",2},{"CYS","SG",1},
        {"GLN","NE2",2},{"HIS","ND1",1},{"HIS","NE2",1},{"LYS","NZ",3},{"SER","OG",1},
        {"THR","OG1",1},{"TRP","NE1",1},{"TYR","OH",1},
    };
    for (const auto& c : kDonor) donor_cap_[key(c.base, c.atom)] = c.cap;

    // --- Acceptor capacity ---
    static const Cap kAcceptor[] = {
        {"G","O6",2},{"U","O2",2},{"U","O4",2},{"C","O2",2},{"T","O2",2},{"T","O4",2},
        {"A","N1",1},{"A","N3",1},{"A","N7",1},{"G","N3",1},{"G","N7",1},{"C","N3",1},
        {"A","O2'",2},{"G","O2'",2},{"C","O2'",2},{"U","O2'",2},{"T","O2'",2},
        {"A","O4'",1},{"G","O4'",1},{"C","O4'",1},{"U","O4'",1},{"T","O4'",1},
        {"A","OP1",3},{"G","OP1",3},{"C","OP1",3},{"U","OP1",3},{"T","OP1",3},
        {"A","OP2",3},{"G","OP2",3},{"C","OP2",3},{"U","OP2",3},{"T","OP2",3},
        {"A","O1P",3},{"G","O1P",3},{"C","O1P",3},{"U","O1P",3},
        {"A","O2P",3},{"G","O2P",3},{"C","O2P",3},{"U","O2P",3},
        {"P","O2",2},{"P","O4",2},{"P","O2'",2},{"P","O4'",1},{"P","OP1",3},{"P","OP2",3},
        {"I","O6",2},{"I","N3",1},{"I","N7",1},{"I","O2'",2},{"I","O4'",1},{"I","OP1",3},{"I","OP2",3},
        {"DA","N1",1},{"DA","N3",1},{"DA","N7",1},{"DA","O4'",1},
        {"DG","O6",2},{"DG","N3",1},{"DG","N7",1},{"DG","O4'",1},
        {"DC","O2",2},{"DC","N3",1},{"DC","O4'",1},{"DT","O2",2},{"DT","O4",2},{"DT","O4'",1},
        {"ALA","O",2},{"ARG","O",2},{"ASN","O",2},{"ASP","O",2},{"CYS","O",2},{"GLN","O",2},
        {"GLU","O",2},{"GLY","O",2},{"HIS","O",2},{"ILE","O",2},{"LEU","O",2},{"LYS","O",2},
        {"MET","O",2},{"PHE","O",2},{"PRO","O",2},{"SER","O",2},{"THR","O",2},{"TRP","O",2},
        {"TYR","O",2},{"VAL","O",2},
        {"ASN","OD1",2},{"ASP","OD1",3},{"ASP","OD2",3},{"GLN","OE1",2},{"GLU","OE1",3},{"GLU","OE2",3},
        {"HIS","ND1",1},{"HIS","NE2",1},{"MET","SD",2},{"SER","OG",2},{"THR","OG1",2},{"TYR","OH",2},
    };
    for (const auto& c : kAcceptor) acceptor_cap_[key(c.base, c.atom)] = c.cap;

    // --- Extend from ligand_hbond_db.json (only keys not already present) ---
    std::ifstream in(ligand_db_path);
    if (!in.is_open()) return;
    nlohmann::json db;
    try {
        db = nlohmann::json::parse(in);
    } catch (const std::exception&) {
        return;
    }
    // Phosphate-oxygen spelling aliases (OP1<->O1P, ...): the DB lists OP-form,
    // but the parser normalizes to O-form. Standard maps already carry both; mirror
    // that for DB-derived modified residues so phosphate H-bonds are found
    // regardless of spelling.
    static const std::unordered_map<std::string, std::string> phos_alias = {
        {"OP1", "O1P"}, {"OP2", "O2P"}, {"OP3", "O3P"},
        {"O1P", "OP1"}, {"O2P", "OP2"}, {"O3P", "OP3"}};
    auto add = [&](std::unordered_map<std::string, int>& table, const std::string& res_name,
                   const std::string& atom, int cap) {
        table.emplace(key(res_name, atom), cap);  // no-op if present
        if (auto it = phos_alias.find(atom); it != phos_alias.end())
            table.emplace(key(res_name, it->second), cap);
    };
    for (const auto& [res_name, res_data] : db.items()) {
        if (res_name == "_metadata" || !res_data.is_object()) continue;
        if (res_data.contains("donors")) {
            for (const auto& d : res_data["donors"])
                add(donor_cap_, res_name, d.value("atom", std::string{}), d.value("capacity", 1));
        }
        if (res_data.contains("acceptors")) {
            for (const auto& a : res_data["acceptors"])
                add(acceptor_cap_, res_name, a.value("atom", std::string{}),
                    a.value("capacity", 1));
        }
    }
}

void HBondChemistry::save_cache(const std::filesystem::path& path) const {
    nlohmann::json root = nlohmann::json::object();
    nlohmann::json donors = nlohmann::json::object();
    nlohmann::json acceptors = nlohmann::json::object();
    for (const auto& [k, cap] : donor_cap_) donors[k] = cap;
    for (const auto& [k, cap] : acceptor_cap_) acceptors[k] = cap;
    root["donor"] = std::move(donors);
    root["acceptor"] = std::move(acceptors);
    std::ofstream out(path);
    out << root.dump(2) << "\n";
}

bool HBondChemistry::load_cache(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in.is_open()) return false;
    nlohmann::json root;
    try {
        in >> root;
    } catch (const nlohmann::json::exception&) {
        return false;
    }
    if (!root.is_object() || !root.contains("donor") || !root.contains("acceptor"))
        return false;
    std::unordered_map<std::string, int> donors, acceptors;
    for (auto it = root["donor"].begin(); it != root["donor"].end(); ++it) {
        if (!it->is_number_integer()) return false;
        donors[it.key()] = it->get<int>();
    }
    for (auto it = root["acceptor"].begin(); it != root["acceptor"].end(); ++it) {
        if (!it->is_number_integer()) return false;
        acceptors[it.key()] = it->get<int>();
    }
    donor_cap_ = std::move(donors);
    acceptor_cap_ = std::move(acceptors);
    return true;
}

}  // namespace pairfinder::algorithms::hbond
