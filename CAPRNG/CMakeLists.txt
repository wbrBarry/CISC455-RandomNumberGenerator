cmake_minimum_required(VERSION 3.1)
project(ca_prng)

set(CMAKE_CXX_STANDARD 14)

# Set optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

add_executable(ca_prng main.cpp makeCA.cpp)