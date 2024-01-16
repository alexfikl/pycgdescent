cmake_minimum_required(VERSION 3.15)
project(${SKBUILD_PROJECT_NAME} LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ version selection")
set(PYBIND11_NEWPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

if(DEFINED ENV{CG_DISABLE_BLAS})
    set(BLAS_LIBRARIES)
    add_compile_definitions(CG_DISABLE_BLAS)
else ()
    find_package(BLAS REQUIRED)
endif()

pybind11_add_module(
    _cg_descent
    src/wrapper/cg_descent.c
    src/wrapper/cg_descent_wrap.cpp
    )
target_link_libraries(_cg_descent PRIVATE ${BLAS_LIBRARIES})
target_compile_definitions(
    _cg_descent
    PRIVATE VERSION_INFO=${SKBUILD_PROJECT_VERSION})

install(TARGETS _cg_descent LIBRARY DESTINATION .)
