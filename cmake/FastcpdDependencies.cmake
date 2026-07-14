include_guard(GLOBAL)

set(FASTCPD_ARMADILLO_INCLUDE_DIR "" CACHE PATH
  "Optional Armadillo include directory for no-wrapper builds")

if(FASTCPD_FETCH_DEPENDENCIES)
  include(FetchContent)

  FetchContent_Declare(fastcpd_armadillo
    URL https://sourceforge.net/projects/arma/files/armadillo-14.4.0.tar.xz
    URL_HASH
      SHA256=023242fd59071d98c75fb015fd3293c921132dc39bf46d221d4b059aae8d79f4)
  FetchContent_GetProperties(fastcpd_armadillo)
  if(NOT fastcpd_armadillo_POPULATED)
    FetchContent_Populate(fastcpd_armadillo)
  endif()
  if(NOT FASTCPD_ARMADILLO_INCLUDE_DIR)
    set(FASTCPD_ARMADILLO_INCLUDE_DIR
      "${fastcpd_armadillo_SOURCE_DIR}/include")
  endif()

  set(ABSL_BUILD_TESTING OFF CACHE BOOL "" FORCE)
  set(ABSL_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
  set(ABSL_PROPAGATE_CXX_STD ON CACHE BOOL "" FORCE)
  FetchContent_Declare(fastcpd_abseil
    URL
      https://github.com/abseil/abseil-cpp/releases/download/20260526.0/abseil-cpp-20260526.0.tar.gz
    URL_HASH
      SHA256=6e1aee535473414164bf83e4ebc40240dec71a4701f8a642d906e95bea1aea0c)
  FetchContent_GetProperties(fastcpd_abseil)
  if(NOT fastcpd_abseil_POPULATED)
    FetchContent_Populate(fastcpd_abseil)
    add_subdirectory(
      "${fastcpd_abseil_SOURCE_DIR}"
      "${fastcpd_abseil_BINARY_DIR}"
      EXCLUDE_FROM_ALL)
  endif()
else()
  find_package(absl 20260526 CONFIG REQUIRED)
endif()

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
  if(NOT FASTCPD_ARMADILLO_INCLUDE_DIR)
    message(FATAL_ERROR
      "FASTCPD_ARMADILLO_INCLUDE_DIR is required when "
      "FASTCPD_USE_ARMA_WRAPPER=OFF")
  endif()

  if(WIN32 AND FASTCPD_FETCH_DEPENDENCIES)
    FetchContent_Declare(fastcpd_openblas
      URL
        https://github.com/OpenMathLib/OpenBLAS/releases/download/v0.3.29/OpenBLAS-0.3.29_x64.zip
      URL_HASH
        SHA256=b42a74d1c9c77bdab2cf2688031b9bc4a322ade71c549427f4950d85fd590fca)
    FetchContent_GetProperties(fastcpd_openblas)
    if(NOT fastcpd_openblas_POPULATED)
      FetchContent_Populate(fastcpd_openblas)
    endif()

    add_library(fastcpd_openblas SHARED IMPORTED GLOBAL)
    set_target_properties(fastcpd_openblas PROPERTIES
      IMPORTED_IMPLIB
        "${fastcpd_openblas_SOURCE_DIR}/lib/libopenblas.lib"
      IMPORTED_LOCATION
        "${fastcpd_openblas_SOURCE_DIR}/bin/libopenblas.dll")

    add_library(BLAS::BLAS INTERFACE IMPORTED)
    set_property(TARGET BLAS::BLAS PROPERTY
      INTERFACE_LINK_LIBRARIES fastcpd_openblas)
    add_library(LAPACK::LAPACK INTERFACE IMPORTED)
    set_property(TARGET LAPACK::LAPACK PROPERTY
      INTERFACE_LINK_LIBRARIES fastcpd_openblas)
  else()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    if(NOT TARGET BLAS::BLAS)
      add_library(BLAS::BLAS INTERFACE IMPORTED)
      set_property(TARGET BLAS::BLAS PROPERTY
        INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    endif()
    if(NOT TARGET LAPACK::LAPACK)
      add_library(LAPACK::LAPACK INTERFACE IMPORTED)
      set_property(TARGET LAPACK::LAPACK PROPERTY
        INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    endif()
  endif()
endif()
