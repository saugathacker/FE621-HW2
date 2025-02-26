cmake_minimum_required(VERSION 3.30)
project(FE621HW2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a header-only library for OptionType
add_library(OptionType INTERFACE)
target_include_directories(OptionType INTERFACE ${CMAKE_SOURCE_DIR})

# Create BinomialTree library
add_library(BinomialTreeLib BinomialTree.cpp)
target_link_libraries(BinomialTreeLib PRIVATE OptionType)

# Create TrinomialTree library (which depends on BinomialTree)
add_library(TrinomialTreeLib TrinomialTree.cpp)
target_link_libraries(TrinomialTreeLib PRIVATE BinomialTreeLib)

# Create main executable (links both libraries)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE BinomialTreeLib TrinomialTreeLib)