# ------------------------------------------------------------------------------
# Fetch H5hut (renamed/updated H5Part library) using CMake FetchContent.
# ------------------------------------------------------------------------------
include(FetchContent)

# HDF5 is required by H5hut; the parallel build is preferred when available.
find_package(HDF5 REQUIRED)

set(H5hut_VERSION "git.53666583861f944b1284eb7324337dfe025242d0"
  CACHE STRING "H5hut version/tag/SHA to fetch")

# If a git ref is requested, build from source; otherwise try find-or-fetch.
if("${H5hut_VERSION}" MATCHES "^git\\.(.+)$")
  set(_h5hut_git_tag "${CMAKE_MATCH_1}")
  FetchContent_Declare(
    H5hut
    GIT_REPOSITORY https://github.com/H5hut/H5hut
    GIT_TAG ${_h5hut_git_tag})
else()
  FetchContent_Declare(
    H5hut
    GIT_REPOSITORY https://github.com/H5hut/H5hut
    GIT_TAG ${H5hut_VERSION}
    FIND_PACKAGE_ARGS ${H5hut_VERSION})
endif()

FetchContent_MakeAvailable(H5hut)

if(TARGET H5hut)
  message(STATUS "H5hut available as target H5hut (source in ${h5hut_SOURCE_DIR})")
else()
  message(FATAL_ERROR "H5hut target was not created")
endif()
