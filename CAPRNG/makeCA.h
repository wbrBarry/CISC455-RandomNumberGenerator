#ifndef makeCA_h
#define makeCA_h

#include <bitset>
#include <chrono>
#include <climits>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void makeRulePattern(const unsigned long& rule,
                     std::vector<bool>& rule_bits_vec);

inline void toBinaryList(const unsigned long& seed,
                         const unsigned long& size,
                         std::vector<bool>& seed_bin);

inline void genNextRow(const std::vector<bool>& prevRow,
                       const std::vector<bool>& rulePattern,
                       std::vector<bool>& nextRow,
                       const size_t& start_col, const size_t& end_col);

inline void initCA(const std::vector<bool>& rulePattern,
                   const std::vector<bool>& seedBinPattern,
                   std::vector<std::vector<bool>>& caMat,
                   const size_t& start_row, const size_t& end_row,
                   const size_t& start_col, const size_t& end_col);
inline void makeCA(const std::vector<bool>& rulePattern,
                   std::vector<std::vector<bool>>& caMat,
                   const size_t& start_row, const size_t& end_row,
                   const size_t& start_col, const size_t& end_col);

void makeCAandReturn(const unsigned long& rule,
                     const unsigned long& seed,
                     const unsigned long& size,
                     std::vector<std::vector<bool>>& caMat,
                     const size_t& start_row, const size_t& end_row,
                     const size_t& start_col, const size_t& end_col,
                     const bool init = false);

#endif