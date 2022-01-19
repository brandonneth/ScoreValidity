###############################################################################
# Copyright (c) 2016-19, Lawrence Livermore National Security, LLC
# and RAJA project contributors. See the RAJA/COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
###############################################################################

set(all_depends RAJA)

if (ENABLE_OPENMP)
  set(all_depends 
    ${all_depends}
    openmp)
endif ()

if (ENABLE_CUDA)
  set(all_depends 
    ${all_depends}
    cuda)
endif ()

blt_add_executable(
  NAME benchmark1.exe
  SOURCES benchmark1.cpp
  DEPENDS_ON ${all_depends})

blt_add_executable(
  NAME benchmark2.exe
  SOURCES benchmark2.cpp
  DEPENDS_ON ${all_depends})

blt_add_executable(
  NAME benchmark3.exe
  SOURCES benchmark3.cpp
  DEPENDS_ON ${all_depends})

blt_add_executable(
  NAME intro-example.exe
  SOURCES intro-example.cpp
  DEPENDS_ON ${all_depends})
