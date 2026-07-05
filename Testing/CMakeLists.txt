if(NOT BUILD_TESTING)
  return()
endif()

set(PLUGIN_TEST_DIR "${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary")
file(MAKE_DIRECTORY "${PLUGIN_TEST_DIR}")

find_package(MPI REQUIRED)

# --------------------------------------------------------------------------
# Test helpers
# --------------------------------------------------------------------------
add_library(pv_meshlessTestLib STATIC TestUtils.cxx)
target_include_directories(pv_meshlessTestLib
                           PRIVATE "${PROJECT_SOURCE_DIR}/Filters")
target_link_libraries(
  pv_meshlessTestLib PRIVATE VTK::CommonCore VTK::ParallelMPI
                             VTK::TestingRendering)

# --------------------------------------------------------------------------
# H5Part parallel writer / reader round-trip test
# --------------------------------------------------------------------------
add_executable(TestH5PartParallelWriter TestH5PartParallelWriter.cxx)
target_include_directories(TestH5PartParallelWriter
                           PRIVATE "${PROJECT_SOURCE_DIR}/Filters")
target_link_libraries(
  TestH5PartParallelWriter PRIVATE pvMeshlessFilters pv_meshlessTestLib
                                   VTK::TestingRendering MPI::MPI_CXX)
vtk_module_autoinit(TARGETS TestH5PartParallelWriter MODULES pvMeshlessFilters
                    VTK::TestingRendering)

foreach(N 1 2 3 4)
  set(test_name "TestH5PartParallelWriter-P${N}")
  if(N EQUAL 1)
    add_test(
      NAME ${test_name}
      COMMAND
        $<TARGET_FILE:TestH5PartParallelWriter> -generateParticles 1000000 -F
        temp.h5 -D ${PLUGIN_TEST_DIR} -T "${PLUGIN_TEST_DIR}" -testName
        ${test_name})
  else()
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${N})
      add_test(
        NAME ${test_name}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestH5PartParallelWriter> -generateParticles 1000000 -F
          temp.h5 -D ${PLUGIN_TEST_DIR} -T "${PLUGIN_TEST_DIR}" -testName
          ${test_name})
    endif()
  endif()
endforeach()
