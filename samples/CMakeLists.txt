project(linalg)
cmake_minimum_required(VERSION 3.0)

add_executable(linalg-test
    linalg-test.cpp
    linalg-test2.cpp
    test-user-defined-conversions.cpp
    test-matrices.cpp
    test-ostream-overloads.cpp
    test-constant-expressions.cpp
    test-higher-order-functions.cpp
    test-quaternions.cpp
)
set_target_properties(linalg-test PROPERTIES 
    CXX_STANDARD 11
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
