cmake_minimum_required(VERSION 3.14)
project(mini-newtonian-nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Release with native optimizations
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall -Wextra")

# Build shared library for Julia ccall
add_library(nbody SHARED
    nbody_accel.cpp
)

target_include_directories(nbody PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Install to lib/ directory
install(TARGETS nbody LIBRARY DESTINATION lib)
