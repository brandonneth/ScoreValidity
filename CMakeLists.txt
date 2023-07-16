###############################################################################
# Copyright (c) 2016-19, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################
cmake_minimum_required(VERSION 3.8)

project( raja-example )

set (ENABLE_TESTS Off CACHE BOOL "")
set (ENABLE_EXAMPLES Off CACHE BOOL "")
set (ENABLE_REPRODUCERS Off CACHE BOOL "")
set (ENABLE_EXERCISES Off CACHE BOOL "")
set (ENABLE_DOCUMENTATION Off CACHE BOOL "")
set (ENABLE_BENCHMARKS Off CACHE BOOL "")
set(BLT_CXX_STD "c++17" CACHE STRING "")
set(ENABLE_OPENMP On)

set(CMAKE_CXX_STANDARD 17)
set(BLT_CXX_STD c++17)


include(../blt/SetupBLT.cmake)

find_package(RAJA REQUIRED)
blt_print_target_properties(TARGET RAJA)

include(FindPkgConfig)
pkg_check_modules(ISL REQUIRED IMPORTED_TARGET isl)
pkg_check_modules(BARVINOK REQUIRED IMPORTED_TARGET barvinok)

if (ENABLE_CUDA)
  set (CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
add_subdirectory(src)
