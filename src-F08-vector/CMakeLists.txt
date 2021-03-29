cmake_minimum_required(VERSION 3.12)

project(NKA Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")

# CMake doesn't set appropriate build-type flags for NAG and Intel.
if (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -Ounsafe -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g90 -gline")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g")
endif()

add_executable(nka_example f90_assert.F90 vector_class.F90 nka_type.F90 grid_vector_type.F90 nka_example.F90)
