cmake_minimum_required(VERSION 3.30)
project(FE621HW2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Create a header-only library for OptionType
add_library(OptionType INTERFACE)
target_include_directories(OptionType INTERFACE ${CMAKE_SOURCE_DIR})

# Create util library (for norm_cdf and other utility functions)
add_library(UtilLib util.cpp)

# Create BSM Library
add_library(BSMLib BSM.cpp)
target_link_libraries(BSMLib PRIVATE OptionType UtilLib)

# Create BarrierBSM Library
add_library(BarrierBSMLib BarrierBSM.cpp)
target_link_libraries(BarrierBSMLib PRIVATE BSMLib UtilLib)

# Create BinomialTree library
add_library(BinomialTreeLib BinomialTree.cpp)
target_link_libraries(BinomialTreeLib PRIVATE OptionType UtilLib)

# Create TrinomialTree library (which depends on BinomialTree)
add_library(TrinomialTreeLib TrinomialTree.cpp)
target_link_libraries(TrinomialTreeLib PRIVATE BinomialTreeLib UtilLib)

# Create main executable
add_executable(main main.cpp)
target_link_libraries(main PRIVATE BinomialTreeLib TrinomialTreeLib BSMLib BarrierBSMLib UtilLib)