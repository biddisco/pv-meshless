cmake_minimum_required(VERSION 3.22...3.30)
project(pv_meshless)

find_package(ParaView REQUIRED)

if(NOT PARAVIEW_USE_MPI)
  message(FATAL_ERROR "pv-meshless requires an MPI enabled ParaView")
endif()

find_package(MPI REQUIRED)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fetch the H5hut library (parallel HDF5-based H5Part successor).
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(FetchH5hut)

set("_paraview_plugin_default_${PROJECT_NAME}" ON)

paraview_plugin_scan(
  PLUGIN_FILES
  "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  ENABLE_BY_DEFAULT
  ON
  PROVIDES_PLUGINS
  plugins
  REQUIRES_MODULES
  required_modules)

foreach(module IN LISTS required_modules)
  if(NOT TARGET "${module}")
    message(
      FATAL_ERROR "Skipping ${PROJECT_NAME}: Missing required module: ${module}"
    )
  endif()
endforeach()

paraview_plugin_build(
  RUNTIME_DESTINATION
  "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION
  "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY
  "${PARAVIEW_PLUGIN_SUBDIR}"
  PLUGINS
  ${plugins})

if(BUILD_TESTING)
  enable_testing()
  add_subdirectory(Testing)
endif()
