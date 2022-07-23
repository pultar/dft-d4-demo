cmake_minimum_required(VERSION 3.21)

# set the project name
project(dftb4-example)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# include directory
set(MY_NAME dftb4-example)
set(PROGRAM_INCLUDES ${D4}/include)

# linking
link_directories(${D4}/lib)

# add the executable
add_executable(${MY_NAME} main.cc)
target_include_directories(${MY_NAME} PUBLIC ${PROGRAM_INCLUDES})
target_link_libraries(${MY_NAME} dftd4)
