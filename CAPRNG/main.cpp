#include <omp.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "makeCA.h"


inline void transpose(const std::vector<std::vector<bool>>& matrix,
                      std::vector<std::vector<bool>> &result) {
    const size_t m = matrix.size();
    const size_t n = matrix[0].size();

#pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = matrix[j][i];
        }
    }
}


inline void parse_rules_grid(const std::string& rules_str,
                             int grid_size,
                             std::vector<std::vector<size_t>>& rules_grid) {
    std::stringstream ss(rules_str);
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            ss >> rules_grid[i][j];
        }
    }
}

inline void writeCSV(const std::string&& filename,
                     const std::vector<std::vector<bool>>& CA_mat) {
    std::ofstream fileH("H_" + filename);

    for (const auto& row : CA_mat) {
        for (size_t i = 0; i < row.size(); ++i) {
            fileH << row[i];
            // if (i < row.size() - 1) {
            //     file << ",";
            // }
        }
        fileH << "\n";
    }

    const size_t m = CA_mat.size();
    const size_t n = CA_mat[0].size();
    fileH.close();

    std::vector<std::vector<bool>> t_mat (n, std::vector<bool>(m));
    transpose(CA_mat, t_mat);
    std::ofstream fileV("V_" +  filename);
    for (const auto& row : t_mat) {
        for (size_t i = 0; i < row.size(); ++i) {
            fileV << row[i];
            // if (i < row.size() - 1) {
            //     file << ",";
            // }
        }
        fileV << "\n";
    }

    fileV.close();
}

inline void testALL(const int& grid_size,
                    const int& ca_size,
                    unsigned long& seed,
                    unsigned long& cycle_length,
                    const std::string& unique_id,
                    std::vector<std::vector<size_t>>& rules_grid) {
    unsigned count_random = 0;
    unsigned period = 0;
    std::vector<std::vector<bool>>
        combined_ca(grid_size * ca_size,
                    std::vector<bool>(grid_size * ca_size));

    // #pragma omp parallel for collapse(2)
    // first round
    // initialize the grid, this round will use the seed
    for (size_t i = 0; i < grid_size; ++i) {
        for (size_t j = 0; j < grid_size; ++j) {
            const size_t start_row = j * ca_size;
            const size_t end_row = start_row + ca_size;

            const size_t start_col = i * ca_size;
            const size_t end_col = start_col + ca_size;

            makeCAandReturn(rules_grid[i][j],
                            seed, ca_size,
                            combined_ca,
                            start_row, end_row,
                            start_col, end_col, true);

            std::ostringstream oss;
            oss << "Origin" << ".txt";
            writeCSV(oss.str(), combined_ca);

        }
    }

    while (period++ < cycle_length) {
        

        for (size_t i = 0; i < grid_size; ++i) {
            for (size_t j = 0; j < grid_size; ++j) {

                makeCAandReturn(rules_grid[i][j],
                                seed, ca_size,
                                combined_ca,
                                0, grid_size * ca_size,
                                0, grid_size * ca_size);
            }
        }

        std::ostringstream oss;
        oss << unique_id << period  <<".txt";
        writeCSV(oss.str(), combined_ca);
    }
}

int main(int argc, char* argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();
    if (argc < 3) {
        std::cerr
            << "Usage: "
            << argv[0]
            << " seed cycle_length [grid_size ca_size rules_grid id] "
            << std::endl;
        return 1;
    }

    unsigned long seed = std::stoul(argv[1]);
    unsigned long cycle_length = std::stoul(argv[2]);
    int grid_size = 3;
    int ca_size = 64;

    if (argc >= 4) {
        grid_size = std::stoi(argv[3]);
    }

    if (argc >= 5) {
        ca_size = std::stoi(argv[4]);
    }

    std::vector<std::vector<size_t>>
        rules_grid(grid_size,
                   std::vector<size_t>(grid_size));

    if (argc >= 6) {
        parse_rules_grid(argv[5], grid_size, rules_grid);
    } else {
        rules_grid = {{26, 28, 25},
                      {22, 11, 18},
                      {25, 26, 28}};
    }
    std::string unique_id = "";
    if (argc >= 7) {
        unique_id = argv[6];
    }
    testALL(grid_size, ca_size, seed,
            cycle_length, unique_id,
            rules_grid);

    // double result = test_this_ca_matrix(combined_ca);

    // auto end_time = std::chrono::high_resolution_clock::now();

    // auto duration_s = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    // std::cout << "Time elapsed: "
    //           << duration_s.count()
    //           << " s"
    //           << std::endl;

    // if (unique_id.size() > 0) {
    //     std::ofstream file(unique_id + ".txt");
    //     file << result;
    //     file.close();
    // }

    // std::cout << result << std::endl;

    return 0;
}