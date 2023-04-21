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
                       const size_t& start_col, const size_t& end_col) {
    // #pragma omp parallel for default(none) shared(prevRow, rulePattern, nextRow, length)
    for (size_t idx = start_col; idx < end_col; ++idx) {
        bool left_cell = prevRow[(idx + end_col - 1) % end_col];
        bool center_cell = prevRow[idx];
        bool right_cell = prevRow[(idx + 1) % end_col];

        int decimal_idx = left_cell * 4 + center_cell * 2 + right_cell;
        // #pragma omp critical
        nextRow[idx] = rulePattern[decimal_idx];
    }
}

inline void initCA(const std::vector<bool>& rulePattern,
                   const std::vector<bool>& seedBinPattern,
                   std::vector<std::vector<bool>>& caMat,
                   const size_t& start_row, const size_t& end_row,
                   const size_t& start_col, const size_t& end_col) {
    for (size_t i = start_row + 1; i < end_row; ++i) {
        genNextRow(caMat[i - 1], rulePattern,
                   caMat[i], start_col,
                   end_col);
    }
}

inline void makeCA(const std::vector<bool>& rulePattern,
                   std::vector<std::vector<bool>>& caMat,
                   const size_t& start_row, const size_t& end_row,
                   const size_t& start_col, const size_t& end_col) {
    for (size_t i = start_row; i < end_row; ++i) {
        for (size_t j = start_col; j < end_col; ++j) {
            bool N = caMat[(i + end_row - 1) % end_row][j];
            bool S = caMat[(i + 1) % end_row][j];
            bool E = caMat[i][(j + 1) % end_col];
            bool W = caMat[i][(j + end_col - 1) % end_col];
            bool C = caMat[i][j];

            caMat[i][j] = rulePattern[0];
            if (rulePattern[1]) {
                caMat[i][j] = N ^ caMat[i][j];
            }
            if (rulePattern[2]) {
                caMat[i][j] = S ^ caMat[i][j];
            }
            if (rulePattern[3]) {
                caMat[i][j] = E ^ caMat[i][j];
            }
            if (rulePattern[4]) {
                caMat[i][j] = W ^ caMat[i][j];
            }
            if (rulePattern[5]) {
                caMat[i][j] = C ^ caMat[i][j];
            }
        }
    }
}

void makeCAandReturn(const unsigned long& rule,
                     const unsigned long& seed,
                     const unsigned long& size,
                     std::vector<std::vector<bool>>& caMat,
                     const size_t& start_row, const size_t& end_row,
                     const size_t& start_col, const size_t& end_col,
                     const bool init) {
    std::vector<bool> rulePattern(8);
    makeRulePattern(rule, rulePattern);

    std::vector<bool> seedBinPattern(size, false);
    toBinaryList(seed, size, seedBinPattern);
    // std::vector<std::vector<bool>>
    //     caMat(size, std::vector<bool>(size));
    // initialize it here
    if (init) {
        genNextRow(seedBinPattern, rulePattern,
                   caMat[start_row],
                   start_col, end_col);

        initCA(rulePattern, seedBinPattern,
               caMat, start_row, end_row,
               start_col, end_col);
    }
    makeCA(rulePattern,
           caMat, start_row, end_row,
           start_col, end_col);
}