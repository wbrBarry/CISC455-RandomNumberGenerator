cmake_minimum_required(VERSION 3.10)

project(ca_fit)

# Use C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set build type to Release and enable highest optimization level
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add your source files
set(SOURCES main.cpp makeCA.cpp)

# Find required libraries
# find_package(PkgConfig REQUIRED)
# pkg_check_modules(FFTW3 REQUIRED fftw3)
# find_package(Boost COMPONENTS math_c99 REQUIRED)
find_package(OpenMP REQUIRED)

# Create the executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link the required libraries
target_link_libraries(${PROJECT_NAME}
    PRIVATE
        # ${FFTW3_LIBRARIES}
        OpenMP::OpenMP_CXX
)

# # Include directories
# target_include_directories(${PROJECT_NAME}
#     PRIVATE
#         ${FFTW3_INCLUDE_DIRS}
# )