# CISC 455  Evolutionary Optimization W23 / group project: Random Number generation
## Group Member:
- Wanqing Li
- Lu Chen
- Yuzhe He
- Baorong Wei

### how to run diehard_test.py

*Step 1. install Dieharder package*
Method 1: You can install dieharder test suite directly from Linux Package Managesment. (Ubuntu)
$ sudo apt-get install dieharder

Method 2: Or you can clone repository["https://github.com/GINARTeam/Diehard-statistical-test"] and build the Diehard Test by yourself using their C source code (using Makefile) (* This method doesn't work on my laptop)

*Step 2. Run python file*
you can edit the test you want by change the *test* list in function *diehard*.
Currently, the test in used is [0,1,2,3,8,9,10,13,15], to keep consistency with CA paper [Tomassini, M., Sipper, M., & Perrenoud, M. (2000). On the generation of high-quality random numbers by two-dimensional cellular automata. IEEE Transactions on computers, 49(10), 1146-1151.]

Installed dieharder tests:
-  Test Number                         Test Name                Test Reliability
- ===============================================================================
-   -d 0                            Diehard Birthdays Test              Good
-   -d 1                               Diehard OPERM5 Test              Good
-   -d 2                    Diehard 32x32 Binary Rank Test              Good
  -d 3                      Diehard 6x8 Binary Rank Test              Good
  -d 4                            Diehard Bitstream Test              Good
  -d 5                                      Diehard OPSO           Suspect
  -d 6                                 Diehard OQSO Test           Suspect
  -d 7                                  Diehard DNA Test           Suspect
  -d 8                Diehard Count the 1s (stream) Test              Good
  -d 9                  Diehard Count the 1s Test (byte)              Good
  -d 10                         Diehard Parking Lot Test              Good
  -d 11         Diehard Minimum Distance (2d Circle) Test             Good
  -d 12         Diehard 3d Sphere (Minimum Distance) Test             Good
  -d 13                             Diehard Squeeze Test              Good
  -d 14                                Diehard Sums Test        Do Not Use
  -d 15                                Diehard Runs Test              Good
  -d 16                               Diehard Craps Test              Good
  -d 17                     Marsaglia and Tsang GCD Test              Good
  -d 100                                STS Monobit Test              Good
  -d 101                                   STS Runs Test              Good
  -d 102                   STS Serial Test (Generalized)              Good
  -d 200                       RGB Bit Distribution Test              Good
  -d 201           RGB Generalized Minimum Distance Test              Good
  -d 202                           RGB Permutations Test              Good
  -d 203                             RGB Lagged Sum Test              Good
  -d 204                RGB Kolmogorov-Smirnov Test Test              Good
  -d 205                               Byte Distribution              Good
  -d 206                                         DAB DCT              Good
  -d 207                              DAB Fill Tree Test              Good
  -d 208                            DAB Fill Tree 2 Test              Good
  -d 209                              DAB Monobit 2 Test              Good
  
  
  # Reference: 
  
