cmake_minimum_required(VERSION 2.8)

project(NKA-F95 Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -Ounsafe -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g90 -gline")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "G95")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
endif()

add_executable(nka_example f90_assert.F90 nka_type.F90 nka_example.F90)
