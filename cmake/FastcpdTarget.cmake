include_guard(GLOBAL)

get_filename_component(_fastcpd_root_dir
  "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(_fastcpd_src_dir "${_fastcpd_root_dir}/src")

add_library(fastcpd "${_fastcpd_root_dir}/cpp/fastcpd.cc")
add_library(fastcpd::fastcpd ALIAS fastcpd)

set_target_properties(fastcpd PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_features(fastcpd PUBLIC cxx_std_17)
target_compile_definitions(fastcpd PRIVATE NO_RCPP)
if(MSVC)
  # fastcpd.cc instantiates every family through the shared template core and
  # can exceed the default COFF section limit in optimized Windows builds.
  target_compile_options(fastcpd PRIVATE /bigobj)
endif()
target_include_directories(fastcpd
  PUBLIC
    $<BUILD_INTERFACE:${_fastcpd_root_dir}/include>
    $<BUILD_INTERFACE:${FASTCPD_ARMADILLO_INCLUDE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    ${_fastcpd_src_dir})

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

unset(_fastcpd_src_dir)
unset(_fastcpd_root_dir)
