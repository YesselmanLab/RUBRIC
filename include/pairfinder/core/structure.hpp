/**
 * @file structure.hpp
 * @brief A parsed structure: chains of residues, with res_id lookup.
 */
#ifndef PAIRFINDER_CORE_STRUCTURE_HPP
#define PAIRFINDER_CORE_STRUCTURE_HPP

#include <string>
#include <vector>

#include <pairfinder/core/chain.hpp>

namespace pairfinder::core {

class Structure {
public:
    Structure() = default;
    explicit Structure(std::string id) : id_(std::move(id)) {}

    const std::string& id() const { return id_; }

    std::vector<Chain>& chains() { return chains_; }
    const std::vector<Chain>& chains() const { return chains_; }

    Chain& add_chain(Chain chain) {
        chains_.push_back(std::move(chain));
        return chains_.back();
    }

    /// Total residue count across all chains.
    std::size_t residue_count() const {
        std::size_t n = 0;
        for (const auto& c : chains_) n += c.residues().size();
        return n;
    }

private:
    std::string id_;
    std::vector<Chain> chains_;
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_STRUCTURE_HPP
