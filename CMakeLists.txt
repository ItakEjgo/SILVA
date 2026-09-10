cmake_minimum_required(VERSION 3.20)
project(mvq LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)

add_compile_options(
    -O3
    -DNDEBUG
    -mcx16
    -march=native
    -DHOMEGROWN
    -Wall
)

include_directories(include)
find_package(Threads REQUIRED)
link_directories(/usr/local/lib)

add_executable(main 
    src/main.cpp
    src/hilbert.c
)

target_link_libraries(main PRIVATE Threads::Threads jemalloc)
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_executable(verify_bench 
    verify_bench.cpp
    src/hilbert.c
)
target_link_libraries(verify_bench PRIVATE Threads::Threads jemalloc)
set_target_properties(verify_bench PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
