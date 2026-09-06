cmake_minimum_required(VERSION 3.16)
project(fastcpd_readme_comparison LANGUAGES CXX)

find_package(fastcpd 1.3 REQUIRED)
set(FPOP_SOURCE_DIR "" CACHE PATH "fpop 2019.8.26 source checkout")
if(NOT EXISTS "${FPOP_SOURCE_DIR}/src/colibri.cc")
  message(FATAL_ERROR "Set FPOP_SOURCE_DIR to the pinned fpop source checkout")
endif()
add_executable(readme_comparison
  comparison.cc
  "${FPOP_SOURCE_DIR}/src/colibri.cc"
  "${FPOP_SOURCE_DIR}/src/liste.cc")
target_compile_features(readme_comparison PRIVATE cxx_std_17)
target_include_directories(readme_comparison PRIVATE "${FPOP_SOURCE_DIR}/src")
target_link_libraries(readme_comparison PRIVATE fastcpd::fastcpd)
