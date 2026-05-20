# Set the minimum required version of CMake
# FATAL_ERROR stops the build immediately if an older version is used.
cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Project name
project("charging_effect")

# Includes the CPM script from cmake folder.  
include("cmake/cpm.cmake") 

# Use CPM to check if ViennaPS is already installed.
# If not, download version 4.4.0 from github and build it as part of the project.
CPMFindPackage(
  NAME ViennaPS
  VERSION 4.4.0
  GIT_REPOSITORY "https://github.com/ViennaTools/ViennaPS"
)

# Tells CMake to compile C++ source file(s) into an executable program.
add_executable(charging_effect src/charging_effect.cpp)

# Tells CMake that your executable needs the ViennaPS library to work.
target_link_libraries(charging_effect ViennaPS)

# Copies "Parameters.txt" from source folder directly into build folder.
configure_file(Parameters.txt ${CMAKE_CURRENT_BINARY_DIR}/Parameters.txt COPYONLY)