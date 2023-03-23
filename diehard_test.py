import math
import random
import decimal
import itertools
import numpy
from math import isqrt
import scipy.special as spc
import scipy.fftpack as sff
import scipy.stats as sst

import subprocess

def decimal_to_binary(decimal_num):
    bi = bin(decimal_num).replace("0b", "")
    if len(bi) < 4:
        bi = "0" * (4 - len(bi)) + bi
    return bi

def sqrt_digits(n, d):
    ''' 
    input param: 
    n - prime number
    d - population size 

    output param:
    digits - irrational number in population size
    '''
    decimal.getcontext().prec = d + 1 # set precision to d+1 decimal places
    x = decimal.Decimal(n)
    y = decimal.Decimal(0)
    for i in range(d):
        y = (x + decimal.Decimal(n)/x) / decimal.Decimal(2)
        if y == x: # check if the square root has converged
            break
        x = y
    digits = str(y).replace('.', '')
    return digits[:d]

def initialize_binary(input,population_size):
    ''' 
    input param: 
    input - seed prime number
    population_size - population size that use input

    output param:
    population - a list of initialized populations in binary version
    '''
    irrational_number = sqrt_digits(input,population_size)
    population = [None] * population_size
    # convert number to a list of digits
    decimal_list = list(map(int, str(irrational_number)))
    # use each digit to generate binary sequence
    for i in range(population_size):
        converted_binary = decimal_to_binary(decimal_list[i])
        # e.g. split binary 1000 to [1,0,0,0] for mutation
        population[i] = list(str(converted_binary))
    return population

def initialize_decimal(input,population_size,chunk_size):
    ''' 
    input param: 
    input - seed prime number
    population_size - population size that use input

    output param:
    population - a list of initialized populations in decimal version
    '''
    irrational_number = sqrt_digits(input,population_size*8)
    chunks = []
    num_str = str(irrational_number)
    # Loop through the string and slice it into chunks of 8 digits
    for i in range(0, len(num_str), chunk_size):
        chunk = num_str[i:i+chunk_size]
        chunks.append(str(chunk))  # convert chunk back to integer and append to the list
    population = chunks
    return population

def diehard(population):
    # 1. save generated population to a bin file
    file_name = "randombytes.bin"
    store_diehard_test = []
    with open(file_name, "w") as f:
        for each_pop in population:
            f.write("%s\n" % each_pop)
    print("write done")

    # 2. save command
    
    # dieharder -D test_name -D pvalues -D assessment -d diehard_opso -f randombyte.bin -c ',' 
    # print result: diehard_opso,0.11388305, PASSED
    command = "dieharder -D test_name -D pvalues -D assessment"

    # 3. Run Dieharder on the file and capture the output
    # can change the test number here
    test = [0,1,2,3,8,9,10,13,15] # craps 跑很久(16),  RGB Generalized Minimum Distance Test 跑很久 (201)
    for test_index in test:
        result = subprocess.run(command.split() + ["-f", file_name, "-d", str(test_index)], stdout=subprocess.PIPE)
        # 4. Print the result
        separated_list=list(result.stdout.decode().strip().replace(' ', '').split("|"))
        store_diehard_test.append(separated_list)
        print(test_index,"test passed\n")
    
    return store_diehard_test

  # main
if __name__ == "__main__":
    input = 17 # seed prime number
    population_size = 1000
    chunk_size = 8 # number of digits of one random number
    mutation_rate = 0.05
    mating_pool_size = 5
    tournament_size = 10

    ### 1. initialization: initialize population with decimals
    population = initialize_decimal(input, population_size, chunk_size)
    ### 2. fitness function: diehard test
    fitness_result = diehard(population)
    print("fitness_result:",fitness_result)

