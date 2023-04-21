# CISC 455  Evolutionary Optimization W23 / group project: Random Number generation
## Group Member:
- Wanqing Li
- Lu Chen
- Yuzhe He
- Baorong Wei

# Introduction
In this work, we use an evolutionary algorithm to find the best rule set for a random number generator based on cellular automata. We evaluate the fitness of the generator's performance using the NIST statistical test suite, which includes a battery of tests for randomness. Our approach produces better generation results compared to prior works that use entropy as a fitness function. Our study demonstrates the effectiveness of using cellular automata and evolutionary algorithms in the design of high-quality random number generators.

# Directory structure
```
.
|-- CAPRNG                                *CAs are created here
|   |-- CMakeLists.txt
|   |-- main.cpp
|   |-- makeCA.cpp
|   `-- makeCA.h
|-- CA_draw.ipynb                         *virtulization tool
|-- modified_nist                         *CAs are evaluated here
|   |-- data
|   |-- experiments                       *Test data are saved here
|   |   |-- AlgorithmTesting              *Test types and it's results
|   |-- include                           *header files
|   |-- makefile
|   |-- obj                               *obj files for .c in src
|   |-- src                               *main body of the suite
|   `-- templates
|-- utilites                              *Exploratory Work
|-- readme.md
`-- run_EA.ipynb                          *run this one evo algo
```

# How to run?
## Requirements:
- Linux-based environment
- CMake
- GCC, G++ (supporting C++17 standard)
- Jupyter or similar environment that runs .ipynb files

For the run_EA.ipynb, you can execute it directly in Google Colab. However, for the CAPRNG folder and modified_nist, you will need to compile them into executables and manually upload them to the Colab runtime. It is possible to have the source files in Colab and compile from there, which was our initial approach, but this method can be unreliable and challenging to troubleshoot. Therefore, in this README, we will outline the most dependable way to run the program.

## Prepare All Necessary Files!
Our project incorporates various languages. The run_EA.ipynb file should be your starting point, as it serves as the main program that launches other programs, extracts their results, and cleans up any residual data they generate. As such, you must first set up the other programs before executing run_EA.ipynb.

### 1. Set up makeCA
```
cd <insert-your-path-here>/CAPRNG
mkdir build
cd ./build
cmake ..
make
```
If everything works correctly, you should see a file named
**ca_fit** in the build folder. Save it for now, 
and we will return to it later.

### 2. Set up NIST
```
cd <insert-your-path-here>/modified_nist
make
```
If everything works correctly, you should see a file named
**assess** in the modified_nist folder. Save it for now, 
and we will return to it later.

### 3. Optional - Persistent Run
If you use Google Drive, create a folder named **CISC455**. Inside this folder, create two more folders: **CAFIT2** and **NIST_tests**. Move the **ca_fit** file into **CAFIT2** and the **assess** file into **NIST_tests**. Additionally, copy the **experiments** and **templates** folders into **NIST_tests**. If you don't use Google Drive, skip this step.

### 3.1 If you completed step 3
You should be able to run the run_EA.ipynb without any issues.

### 3.2 If you didn't complete step 3
Manually upload the **ca_fit** and **assess** files, as well as the **experiments** and **templates** folders, into the runtime. These files will be lost once your session ends, and you will need to repeat this process the next time you run the program.


# Reference: 
CA1 <br/>
https://sci-hub.se/10.1109/12.888056<br/>
Tomassini, M., Sipper, M., & Perrenoud, M. (2000). On the generation of high-quality random numbers by two-dimensional cellular automata. IEEE Transactions on computers, 49(10), 1146-1151. <br/>

CA2<br/>
https://sci-hub.se/10.1109/melcon.2006.1653219<br/>
Szaban, M., Seredynski, F., & Bouvry, P. (2006, May). Evolving collective behavior of cellular automata for cryptography. In MELECON 2006-2006 IEEE Mediterranean Electrotechnical Conference (pp. 799-802). IEEE.<br/>

Fitness:<br/>
Dieharder github<br/>
https://github.com/GINARTeam/Diehard-statistical-test<br/>

NIST test github<br/>
https://gist.github.com/StuartGordonReid<br/>


