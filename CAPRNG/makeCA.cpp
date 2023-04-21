#include "makeCA.h"

#include <bitset>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void makeRulePattern(const unsigned long& rule,
                     std::vector<bool>& rule_bits_vec) {
    std::bitset<8> rule_bits(rule % 256);
    for (size_t i = 0; i < 8; ++i) {
        rule_bits_vec[i] = rule_bits[i];
    }
}

inline void toBinaryList(const unsigned long& seed,
                         const unsigned long& size,
                         std::vector<bool>& seed_bin) {
    unsigned long i = 0;
    if (size > sizeof(unsigned long) * CHAR_BIT) {
        i = sizeof(unsigned long) * CHAR_BIT - 1;
    } else {
        i = size - 1;
    }

    for (; i > 0; --i) {
        seed_bin[size - 1 - i] = (seed >> (i - 1)) & 1;
    }
}

inline void genNextRow(const std::vector<bool>& prevRow,
                       const std::vector<bool>& rulePattern,
                       std::vector<bool>& nextRow,
                       const unsigned long& size) {
    // #pragma omp parallel for default(none) shared(prevRow, rulePattern, nextRow, length)
    for (unsigned long idx = 0; idx < size; ++idx) {
        bool left_cell = prevRow[(idx - 1 + size) % size];
        bool center_cell = prevRow[idx];
        bool right_cell = prevRow[(idx + 1) % size];

        int decimal_idx = left_cell * 4 + center_cell * 2 + right_cell;
        // #pragma omp critical
        nextRow[idx] = rulePattern[decimal_idx];
    }
}

inline void makeCA(const unsigned long& size,
                   const std::vector<bool>& rulePattern,
                   const std::vector<bool>& seedBinPattern,
                   std::vector<std::vector<bool>>& caMat) {
    genNextRow(seedBinPattern, rulePattern, caMat[0], size);
    for (unsigned long i = 1; i < size; ++i) {
        genNextRow(caMat[i - 1], rulePattern, caMat[i], size);
    }
}

std::vector<std::vector<bool>> makeCAandReturn(const unsigned long& rule,
                                               const unsigned long& seed,
                                               const unsigned long& size) {
    std::vector<bool> rulePattern(8);
    makeRulePattern(rule, rulePattern);

    std::vector<bool> seedBinPattern(size, false);
    toBinaryList(seed, size, seedBinPattern);

    std::vector<std::vector<bool>>
        caMat(size, std::vector<bool>(size));
    makeCA(size, rulePattern, seedBinPattern, caMat);

    return caMat;
}