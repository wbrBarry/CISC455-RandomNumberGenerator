from CA_class import CA
from evaluation import fitness
import numpy as np
class CA_grid:
    
    def __init__(self, grid_rules, seed = 123, size = 16) -> None:
        self.grid_rules = grid_rules
        self.seed = seed
        self.size = size
        self.fitnessValue = 0
        
    def evalFitness(self):
        length = len(self.grid_rules)
        caGrid = []
        for row in self.grid_rules:
            caRow = []
            for rule in row:  
                matrix = CA(rule, self.seed, self.size).ret_np_matrix()
                
                if len(caRow) == 0:
                    caRow = matrix
                    
                else:
                    caRow = np.hstack((caRow, matrix))
                
            if len(caGrid) == 0:
                caGrid = caRow
            else:
                caGrid = np.vstack((caGrid, caRow))
        # Iterate through rows
        for row in caGrid:
            row_bit_string = "".join(str(bit) for bit in row)
            self.fitnessValue = fitness(row_bit_string)

        # Iterate through columns
        for col in range(caGrid.shape[1]):
            column = caGrid[:, col]
            col_bit_string = "".join(str(bit) for bit in column)
            self.fitnessValue = fitness(col_bit_string)
        
        print(caGrid)
        print(self.fitnessValue)

                
                
            
        
        
            
        
                
    
    def display(self):
        for row in self.grid_rules:
            print(row)
        print()
        

if __name__ == "__main__":
    import random as rd
    import time
    start_time = time.perf_counter()
    CA_grid([
            [rd.randint(0, 255) for _ in range(8)] for _ in range(8)
        ]).evalFitness()
    end_time = time.perf_counter()
    elapsed_time = end_time - start_time
    print(f"Elapsed time: {elapsed_time:.4f} seconds")
