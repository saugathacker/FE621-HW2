cmake_minimum_required(VERSION 3.16)
project(FE621HW2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only library for OptionType
add_library(OptionType INTERFACE)
target_include_directories(OptionType INTERFACE ${CMAKE_SOURCE_DIR})

# Utility Library (Common Functions)
add_library(UtilLib STATIC util.cpp)
target_link_libraries(UtilLib PUBLIC OptionType)  # PUBLIC allows all users of UtilLib to access OptionType

# Black-Scholes Model Library
add_library(BSMLib STATIC BSM.cpp)
target_link_libraries(BSMLib PUBLIC UtilLib OptionType)  # PUBLIC allows BSM dependents to access UtilLib

# Option Data Library
add_library(OptionDataLib STATIC OptionData.cpp)
target_link_libraries(OptionDataLib PUBLIC UtilLib BSMLib BinomialTreeLib OptionType)

# Ticker Library
add_library(TickerLib STATIC Ticker.cpp)
target_link_libraries(TickerLib PUBLIC OptionDataLib)

# FileHandler Library
add_library(FileHandlerLib STATIC FileHandler.cpp)
target_link_libraries(FileHandlerLib PUBLIC TickerLib)

# Barrier BSM Library
add_library(BarrierBSMLib STATIC BarrierBSM.cpp)
target_link_libraries(BarrierBSMLib PUBLIC BSMLib UtilLib)

# Binomial Tree Library
add_library(BinomialTreeLib STATIC BinomialTree.cpp)
target_link_libraries(BinomialTreeLib PUBLIC UtilLib OptionType)

# Trinomial Tree Library (which depends on BinomialTree)
add_library(TrinomialTreeLib STATIC TrinomialTree.cpp)
target_link_libraries(TrinomialTreeLib PUBLIC BinomialTreeLib UtilLib)

# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC 
    OptionType 
    UtilLib 
    BSMLib 
    BinomialTreeLib 
    TrinomialTreeLib 
    BarrierBSMLib 
    OptionDataLib 
    TickerLib 
    FileHandlerLib)
