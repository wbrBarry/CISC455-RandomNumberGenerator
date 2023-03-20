# Initializing the population
#
# this is initializing a grid where each cell
# in the grid represent a rule
# like the following
# [ 26 27 29
#   1   3  4
#   3  4  1 ]
#

import random as rd
from CA_grid import CA_grid
def permutation(pop_size : int, chrom_size: int)-> list[CA_grid]:


    # return a CA_grid object and each of them contains a rules map
    return [
        CA_grid([
            [rd.randint(0, 255) for _ in range(chrom_size)] for _ in range(chrom_size)
        ]) for _ in range(pop_size)
    ]


if __name__ == "__main__":
    
    ca_grids = permutation(2, 8)
    
    for i, ca_grid in enumerate(ca_grids):
        print(f"CA_grid {i + 1}:")
        ca_grid.display()