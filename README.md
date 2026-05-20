# Charging Effect on Ion Trajectories during Plasma Etching

This project simulates the self-consistent charging of high-aspect-ratio trenches during plasma etching. Currently it is built as a standalone module utilizing the [ViennaPS](https://github.com/ViennaTools/ViennaPS) simulation framework to calculate surface charge accumulation, solve localized electric fields, and trace deflected ion trajectories.

## Prerequisites
* **CMake** (v3.18 or higher)
* A C++17 compatible compiler
* [ParaView](https://www.paraview.org/) (Recommended for visualizing `.vtp` outputs)

## Building and Compilation

This project uses CPM to automatically fetch and link ViennaPS. To generate the build environment and compile the project, run the following from the project root:

```bash
# 1. Generate the build files (downloads dependencies on first run)
cmake -S . -B build

# 2. Compile the executable
cmake --build build

# 3. Run the executable
./build/charging_effect