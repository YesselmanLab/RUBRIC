/**
 * @file saenger.cpp
 * @brief Saenger 28-type table + matching.
 */
#include <pairfinder/algorithms/saenger.hpp>

#include <cctype>
#include <set>
#include <string>
#include <utility>

namespace pairfinder::algorithms::classification {

namespace {

struct SaengerType {
    const char* base1;
    const char* base2;
    std::vector<std::pair<const char*, const char*>> hbonds;  // (base1_atom, base2_atom)
    const char* lw_class;
};

// The 28 Saenger types, order preserved.
const std::vector<SaengerType>& saenger_types() {
    static const std::vector<SaengerType> t = {
        {"A", "A", {{"N6", "N1"}, {"N1", "N6"}}, "tWW"},
        {"A", "A", {{"N7", "N6"}, {"N6", "N7"}}, "tHH"},
        {"G", "G", {{"O6", "N1"}, {"N1", "O6"}}, "tWW"},
        {"G", "G", {{"N2", "N3"}, {"N3", "N2"}}, "tSS"},
        {"A", "A", {{"N7", "N6"}, {"N6", "N1"}}, "tHW"},
        {"G", "G", {{"N1", "O6"}, {"N2", "N7"}}, "cWH"},
        {"G", "G", {{"N7", "N1"}, {"O6", "N2"}}, "tWH"},
        {"A", "G", {{"N6", "O6"}, {"N1", "N1"}}, "cWW"},
        {"A", "G", {{"N6", "O6"}, {"N7", "N1"}}, "cHW"},
        {"A", "G", {{"N6", "N3"}, {"N1", "N2"}}, "tWS"},
        {"A", "G", {{"N6", "N3"}, {"N7", "N2"}}, "tHS"},
        {"U", "U", {{"N3", "O4"}, {"O4", "N3"}}, "tWW"},
        {"U", "U", {{"O2", "N3"}, {"N3", "O2"}}, "tWW"},
        {"C", "C", {{"O2", "N4"}, {"N4", "O2"}}, "tWW"},
        {"C", "C", {{"O2", "N4"}, {"N3", "N3"}, {"N4", "O2"}}, "tWW"},
        {"U", "U", {{"O2", "N3"}, {"N3", "O4"}}, "cWW"},
        {"C", "U", {{"O2", "O4"}, {"N3", "N3"}, {"N4", "O2"}}, "tWW"},
        {"C", "U", {{"N3", "N3"}, {"N4", "O4"}}, "cWW"},
        {"G", "C", {{"O6", "N4"}, {"N1", "N3"}, {"N2", "O2"}}, "cWW"},
        {"A", "U", {{"N6", "O4"}, {"N1", "N3"}}, "cWW"},
        {"A", "U", {{"N6", "O2"}, {"N1", "N3"}}, "tWW"},
        {"G", "C", {{"N1", "O2"}, {"N2", "N3"}}, "tWW"},
        {"A", "U", {{"N6", "O4"}, {"N7", "N3"}}, "cWH"},
        {"A", "U", {{"N6", "O2"}, {"N7", "N3"}}, "tWH"},
        {"A", "C", {{"N6", "N3"}, {"N7", "N4"}}, "tWH"},
        {"A", "C", {{"N6", "N3"}, {"N1", "N4"}}, "tWW"},
        {"G", "U", {{"N1", "O4"}, {"O6", "N3"}}, "tWW"},
        {"G", "U", {{"N1", "O2"}, {"O6", "N3"}}, "cWW"},
    };
    return t;
}

std::string norm_upper(const std::string& s, const core::BaseTyping& typing) {
    std::string n = typing.normalize(s);
    for (char& c : n) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return n;
}

}  // namespace

SaengerMatch match_saenger(const core::Residue& res1, const core::Residue& res2,
                           const std::vector<hbond::HBond>& hbonds,
                           const core::BaseTyping& typing) {
    SaengerMatch out;
    if (hbonds.empty()) return out;
    const std::string base1 = norm_upper(res1.base_type(), typing);
    const std::string base2 = norm_upper(res2.base_type(), typing);

    std::set<std::pair<std::string, std::string>> fwd, rev;
    for (const auto& hb : hbonds) {
        if (hb.donor_res_id == res1.res_id()) {
            fwd.insert({hb.donor_atom, hb.acceptor_atom});
            rev.insert({hb.acceptor_atom, hb.donor_atom});
        } else if (hb.donor_res_id == res2.res_id()) {
            fwd.insert({hb.acceptor_atom, hb.donor_atom});
            rev.insert({hb.donor_atom, hb.acceptor_atom});
        }
    }

    auto subset = [](const std::vector<std::pair<const char*, const char*>>& pat,
                     const std::set<std::pair<std::string, std::string>>& s) {
        for (const auto& [a, b] : pat)
            if (s.find({a, b}) == s.end()) return false;
        return true;
    };
    auto faces = [](const SaengerType& st, bool swapped) -> std::pair<char, char> {
        char fb1 = st.lw_class[1], fb2 = st.lw_class[2];
        return swapped ? std::pair{fb2, fb1} : std::pair{fb1, fb2};
    };

    for (const auto& st : saenger_types()) {
        if (base1 == st.base1 && base2 == st.base2 && subset(st.hbonds, fwd)) {
            auto [f1, f2] = faces(st, false);
            return {true, f1, f2};
        }
        if (base1 == st.base2 && base2 == st.base1 && subset(st.hbonds, rev)) {
            auto [f1, f2] = faces(st, true);
            return {true, f1, f2};
        }
    }
    return out;
}

}  // namespace pairfinder::algorithms::classification
