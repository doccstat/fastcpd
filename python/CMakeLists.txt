cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(fastcpd_python LANGUAGES CXX)

include(GNUInstallDirs)

option(FASTCPD_FETCH_DEPENDENCIES
  "Fetch the C++ dependencies needed for a self-contained build" ON)

include(../src/standalone/cmake/FastcpdDependencies.cmake)
include(../src/standalone/cmake/FastcpdTarget.cmake)

set(PYBIND11_FINDPYTHON ON)
find_package(Python 3.10 COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(interface MODULE pybind.cc)
target_link_libraries(interface PRIVATE fastcpd::fastcpd)
target_compile_features(interface PRIVATE cxx_std_17)

install(FILES
  __init__.py
  confidence.py
  py.typed
  segmentation.py
  variance_estimation.py
  DESTINATION fastcpd
  COMPONENT python)
install(TARGETS interface
  LIBRARY DESTINATION fastcpd COMPONENT python
  RUNTIME DESTINATION fastcpd COMPONENT python)
if(WIN32 AND FASTCPD_FETCH_DEPENDENCIES)
  install(FILES "${fastcpd_openblas_SOURCE_DIR}/bin/libopenblas.dll"
    DESTINATION fastcpd
    COMPONENT python)
endif()
