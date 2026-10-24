/**
 * @file chain.hpp
 * @brief A chain: ordered residues sharing a chain id.
 */
#ifndef PAIRFINDER_CORE_CHAIN_HPP
#define PAIRFINDER_CORE_CHAIN_HPP

#include <string>
#include <vector>

#include <pairfinder/core/residue.hpp>

namespace pairfinder::core {

class Chain {
public:
    Chain() = default;
    explicit Chain(std::string id) : id_(std::move(id)) {}

    const std::string& id() const { return id_; }

    std::vector<Residue>& residues() { return residues_; }
    const std::vector<Residue>& residues() const { return residues_; }

    Residue& add_residue(Residue residue) {
        residues_.push_back(std::move(residue));
        return residues_.back();
    }

private:
    std::string id_;
    std::vector<Residue> residues_;
};

}  // namespace pairfinder::core

#endif  // PAIRFINDER_CORE_CHAIN_HPP
