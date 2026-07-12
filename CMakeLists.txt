cmake_minimum_required(VERSION 3.16)

project(fastcpd
  VERSION 1.1.2
  LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(FASTCPD_BUILD_EXAMPLES "Build fastcpd C++ examples" ON)
set(FASTCPD_ARMADILLO_INCLUDE_DIR "" CACHE PATH
  "Optional Armadillo include directory for no-wrapper builds")
if(FASTCPD_ARMADILLO_INCLUDE_DIR)
  set(_fastcpd_default_use_arma_wrapper OFF)
else()
  set(_fastcpd_default_use_arma_wrapper ON)
endif()
option(FASTCPD_USE_ARMA_WRAPPER
  "Use libarmadillo instead of linking BLAS/LAPACK directly"
  ${_fastcpd_default_use_arma_wrapper})
unset(_fastcpd_default_use_arma_wrapper)

if(FASTCPD_USE_ARMA_WRAPPER)
  find_package(Armadillo REQUIRED)
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  if(NOT FASTCPD_ARMADILLO_INCLUDE_DIR)
    message(FATAL_ERROR
      "FASTCPD_ARMADILLO_INCLUDE_DIR is required when "
      "FASTCPD_USE_ARMA_WRAPPER=OFF")
  endif()
endif()
find_package(absl 20260526 CONFIG REQUIRED)

add_library(fastcpd cpp/fastcpd.cc)
add_library(fastcpd::fastcpd ALIAS fastcpd)

target_compile_features(fastcpd PUBLIC cxx_std_17)
target_compile_definitions(fastcpd PRIVATE NO_RCPP)
target_include_directories(fastcpd
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${FASTCPD_ARMADILLO_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src)

if(FASTCPD_USE_ARMA_WRAPPER)
  if(TARGET Armadillo::Armadillo)
    target_link_libraries(fastcpd PUBLIC Armadillo::Armadillo)
  else()
    target_include_directories(fastcpd PUBLIC ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(fastcpd PUBLIC ${ARMADILLO_LIBRARIES})
  endif()
else()
  target_compile_definitions(fastcpd PUBLIC ARMA_DONT_USE_WRAPPER)
  target_link_libraries(fastcpd PUBLIC LAPACK::LAPACK BLAS::BLAS)
endif()

target_link_libraries(fastcpd PRIVATE
  absl::inlined_vector
  absl::prefetch)

if(FASTCPD_BUILD_EXAMPLES)
  add_executable(fastcpd_cpp_example examples/cpp/simple.cc)
  target_link_libraries(fastcpd_cpp_example PRIVATE fastcpd::fastcpd)
endif()

install(TARGETS fastcpd
  EXPORT fastcpdTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/fastcpdConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)
configure_package_config_file(
  cmake/fastcpdConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/fastcpdConfig.cmake"
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastcpd)

install(EXPORT fastcpdTargets
  NAMESPACE fastcpd::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastcpd)
install(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/fastcpdConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/fastcpdConfigVersion.cmake"
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/fastcpd)
