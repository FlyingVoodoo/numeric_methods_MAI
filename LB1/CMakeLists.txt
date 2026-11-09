cmake_minimum_required(VERSION 3.16)

project(LW1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(LW1
    main.cpp
)

target_include_directories(LW1 PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_options(LW1 PRIVATE
    -Wall
    -Wextra
    -Wpedantic
)