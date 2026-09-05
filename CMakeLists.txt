cmake_minimum_required(VERSION 3.16)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project(fastcpd
  VERSION 1.3.0
  LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(FASTCPD_BUILD_EXAMPLES "Build fastcpd C++ examples" ON)
option(FASTCPD_FETCH_DEPENDENCIES
  "Fetch the C++ dependencies needed for a self-contained build" OFF)
option(FASTCPD_INSTALL_CPP "Install the standalone C++ library" ON)

include(cmake/FastcpdDependencies.cmake)
include(cmake/FastcpdTarget.cmake)

if(FASTCPD_BUILD_EXAMPLES)
  add_executable(fastcpd_mean_change examples/cpp/mean_change.cc)
  target_link_libraries(fastcpd_mean_change PRIVATE fastcpd::fastcpd)

  add_executable(fastcpd_interface_contract examples/cpp/interface_contract.cc)
  target_link_libraries(fastcpd_interface_contract PRIVATE fastcpd::fastcpd)
endif()

if(FASTCPD_INSTALL_CPP)
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
endif()
