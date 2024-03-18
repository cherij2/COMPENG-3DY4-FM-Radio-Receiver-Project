#!/bin/bash

# Navigate to the build directory
cd build

# Remove any existing build artifacts
rm -rf *

# Run CMake to configure the project
cmake ..

# Build the project using make
make
