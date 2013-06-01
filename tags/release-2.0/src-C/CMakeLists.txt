cmake_minimum_required(VERSION 2.8)

project(NKA-C C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(nka_example nonlinear_krylov_accelerator.c nka_example.c)
target_link_libraries(nka_example m)
