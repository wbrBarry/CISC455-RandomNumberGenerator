#include <iostream>
#include <sstream>
#include <vector>

#include "makeCA.h"

inline void binary_to_double(const std::vector<bool>& bits) {
    std::ofstream fileD("rdDouble.bin", std::ios::app);

    const size_t max_bits = 52;  // Maximum bits for the double representation

    unsigned itr = 0;
    while (itr < bits.size()) {
        unsigned i = 0;
        double result = 0.0;
        for (; i < bits.size() && i < max_bits; ++i) {
            result += bits[i] * (1.0 / (1 << (i + 1)));
        }
        fileD << result
              << std::endl;
        itr += i;
    }

    fileD.close();
}

inline void binary_to_uint(const std::vector<bool>& bits, unsigned& count) {
    std::ofstream fileU("rdUint.bin", std::ios::app);

    const size_t max_bits = 8 * sizeof(unsigned);  // Maximum bits for the unsigned integer representation

    unsigned itr = 0;
    while (itr < bits.size()) {
        unsigned i = 0;
        unsigned result = 0;
        for (; i < bits.size() && i < max_bits; ++i) {
            if (bits[i]) {
                result |= 1 << (bits.size() - 1 - i);
            }
        }
        fileU << result
              << std::endl;
        ++count;
        itr += i;
    }

    fileU.close();
}

inline void transpose(const std::vector<std::vector<bool>>& matrix, std::vector<std::vector<bool>> result) {
    for (int i = 0; i < matrix[0].size(); i++) {
        for (int j = 0; j < matrix.size(); j++) {
            result[i][j] = matrix[j][i];
        }
    }
}

inline void gen_randomNumber(std::vector<std::vector<bool>>& combined_ca,
                             const unsigned& num_random, unsigned& count) {
    // rd in row
    for (const auto& row : combined_ca) {
        binary_to_double(row);
        binary_to_uint(row, count);
    }

    // transpose
    int m = combined_ca.size();
    int n = combined_ca[0].size();
    std::vector<std::vector<bool>> result(n, std::vector<bool>(m));

    transpose(combined_ca, result);

    // rd in col
    for (const auto& row : result) {
        binary_to_double(row);
        binary_to_uint(row, count);
    }
}

inline void update_rules_grid(std::vector<std::vector<int>>& rules_grid) {
    std::vector<std::vector<int>> updated_rules_grid = rules_grid;

    int grid_size = rules_grid.size();
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            int north = (i - 1 + grid_size) % grid_size;
            int south = (i + 1) % grid_size;
            int west = (j - 1 + grid_size) % grid_size;
            int east = (j + 1) % grid_size;

            unsigned N = rules_grid[north][j];
            unsigned S = rules_grid[south][j];
            unsigned E = rules_grid[i][east];
            unsigned W = rules_grid[i][west];
            unsigned C = rules_grid[i][j];

            std::vector<bool> rule_bits_vec(8);
            makeRulePattern(C, rule_bits_vec);

            int operation = rule_bits_vec[5] ? 1 : 0;  // 0 for XOR, 1 for XNOR

            int new_rule = 0;

            if (operation == 0) {  // XOR
                new_rule = rule_bits_vec[1] * N ^
                           rule_bits_vec[2] * S ^
                           rule_bits_vec[3] * E ^
                           rule_bits_vec[4] * W ^
                           rule_bits_vec[5] * C;
            } else {  // XNOR
                new_rule = rule_bits_vec[1] * ~N ^
                           rule_bits_vec[2] * ~S ^
                           rule_bits_vec[3] * ~E ^
                           rule_bits_vec[4] * ~W ^
                           rule_bits_vec[5] * ~C;
            }

            if (rule_bits_vec[0]) {
                new_rule += 1;
            }

            if (!rule_bits_vec[6]) {
                new_rule += 1;
            }

            // if (!rule_bits_vec[7]) {
            //     new_rule += 1;
            //     new_rule += 1;
            // }

            updated_rules_grid[i][j] = new_rule & 0xFF;  // Keep only the last 8 bits
        }
    }

    rules_grid = updated_rules_grid;
}

inline void parse_rules_grid(const std::string& rules_str,
                             int grid_size,
                             std::vector<std::vector<int>>& rules_grid) {
    std::stringstream ss(rules_str);
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            ss >> rules_grid[i][j];
        }
    }
}

inline void writeCSV(const std::string& filename,
                     const std::vector<std::vector<bool>>& CA_mat) {
    std::ofstream file(filename);

    for (const auto& row : CA_mat) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i < row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " seed num_random [grid_size ca_size rules_grid]" << std::endl;
        return 1;
    }

    unsigned long seed = std::stoul(argv[1]);
    unsigned long num_random = std::stoul(argv[2]);

    int grid_size = 3;
    int ca_size = 64;

    if (argc >= 4) {
        grid_size = std::stoi(argv[3]);
    }

    if (argc >= 5) {
        ca_size = std::stoi(argv[4]);
    }

    std::vector<std::vector<int>> rules_grid(grid_size, std::vector<int>(grid_size));

    if (argc >= 6) {
        parse_rules_grid(argv[5], grid_size, rules_grid);
    } else {
        rules_grid = {{26, 28, 25},
                      {22, 11, 18},
                      {25, 26, 28}};
    }

    unsigned count_random = 0;
    unsigned period = 0;
    std::vector<std::vector<bool>> combined_ca(grid_size * ca_size, std::vector<bool>(grid_size * ca_size));
    while (count_random <= num_random) {
        // write to the combined_ca
        // for (int i = 0; i < grid_size; ++i) {
        //     for (int j = 0; j < grid_size; ++j) {
        //         auto ca = makeCAandReturn(rules_grid[i][j], seed, ca_size);
        //         for (int x = 0; x < ca_size; ++x) {
        //             for (int y = 0; y < ca_size; ++y) {
        //                 combined_ca[i * ca_size + x][j * ca_size + y] = ca[x][y];
        //             }
        //         }
        //     }
        // }
        // gen_randomNumber(combined_ca, num_random, count_random);

        for (const auto& rules : rules_grid) {
            for (const auto& rule : rules) {
                std::cout << rule << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        update_rules_grid(rules_grid);

        std::ostringstream oss;
        oss << "output" << period++ << ".csv";
        // writeCSV(oss.str(), combined_ca);

        // // test
        if (period > 500000) {
            break;
        }
    }

    // Continue with the rest of the algorithm

    return 0;
}