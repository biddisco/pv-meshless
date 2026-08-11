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
set_property(TARGET TestH5PartParallelWriter APPEND PROPERTY BUILD_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartParallelWriter APPEND PROPERTY INSTALL_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartParallelWriter PROPERTY BUILD_RPATH_USE_LINK_PATH TRUE)
set_property(TARGET TestH5PartParallelWriter PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

foreach(N 1 2 3 4)
  set(test_name "TestH5PartParallelWriter-P${N}")
  set(test_added FALSE)
  if(N EQUAL 1)
    add_test(
      NAME ${test_name}
      COMMAND
        $<TARGET_FILE:TestH5PartParallelWriter> -generateParticles 1000000 -F
        temp.h5 -D ${PLUGIN_TEST_DIR} -T "${PLUGIN_TEST_DIR}" -testName
        ${test_name})
    set(test_added TRUE)
  else()
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${N})
      add_test(
        NAME ${test_name}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestH5PartParallelWriter> -generateParticles 1000000 -F
          temp.h5 -D ${PLUGIN_TEST_DIR} -T "${PLUGIN_TEST_DIR}" -testName
          ${test_name})
      set(test_added TRUE)
    endif()
  endif()
  if(test_added)
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:pvMeshlessFilters>:$ENV{LD_LIBRARY_PATH}")
  endif()
endforeach()

# --------------------------------------------------------------------------
# H5Part reader V2 partition-box parallel read test
# --------------------------------------------------------------------------
add_executable(TestH5PartReaderV2Partitioning TestH5PartReaderV2Partitioning.cxx)
target_include_directories(
  TestH5PartReaderV2Partitioning
  PRIVATE "${PROJECT_SOURCE_DIR}/Filters"
          "$<TARGET_PROPERTY:H5hut,INTERFACE_INCLUDE_DIRECTORIES>"
          "${HDF5_INCLUDE_DIRS}")
target_link_libraries(
  TestH5PartReaderV2Partitioning
  PRIVATE pvMeshlessFilters VTK::TestingRendering MPI::MPI_CXX H5hut
          ${HDF5_LIBRARIES})
vtk_module_autoinit(TARGETS TestH5PartReaderV2Partitioning MODULES
                    pvMeshlessFilters VTK::TestingRendering)
set_property(TARGET TestH5PartReaderV2Partitioning APPEND PROPERTY BUILD_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2Partitioning APPEND PROPERTY INSTALL_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2Partitioning PROPERTY BUILD_RPATH_USE_LINK_PATH TRUE)
set_property(TARGET TestH5PartReaderV2Partitioning PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

foreach(N 1 2 4)
  set(test_name "TestH5PartReaderV2Partitioning-P${N}")
  set(test_added FALSE)
  if(N EQUAL 1)
    add_test(
      NAME ${test_name}
      COMMAND $<TARGET_FILE:TestH5PartReaderV2Partitioning> -T
              "${PLUGIN_TEST_DIR}")
    set(test_added TRUE)
  else()
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${N})
      add_test(
        NAME ${test_name}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestH5PartReaderV2Partitioning> -T "${PLUGIN_TEST_DIR}")
      set(test_added TRUE)
    endif()
  endif()
  if(test_added)
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:pvMeshlessFilters>:$ENV{LD_LIBRARY_PATH}")
  endif()
endforeach()

# --------------------------------------------------------------------------
# H5Part reader V2 MaxParticlesPerRank contiguous vs strided sampling test
# --------------------------------------------------------------------------
add_executable(TestH5PartReaderV2StridedSampling
               TestH5PartReaderV2StridedSampling.cxx)
target_include_directories(
  TestH5PartReaderV2StridedSampling
  PRIVATE "${PROJECT_SOURCE_DIR}/Filters"
          "$<TARGET_PROPERTY:H5hut,INTERFACE_INCLUDE_DIRECTORIES>"
          "${HDF5_INCLUDE_DIRS}")
target_link_libraries(
  TestH5PartReaderV2StridedSampling
  PRIVATE pvMeshlessFilters VTK::TestingRendering MPI::MPI_CXX H5hut
          ${HDF5_LIBRARIES})
vtk_module_autoinit(TARGETS TestH5PartReaderV2StridedSampling MODULES
                    pvMeshlessFilters VTK::TestingRendering)
set_property(TARGET TestH5PartReaderV2StridedSampling APPEND PROPERTY BUILD_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2StridedSampling APPEND PROPERTY INSTALL_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2StridedSampling PROPERTY BUILD_RPATH_USE_LINK_PATH TRUE)
set_property(TARGET TestH5PartReaderV2StridedSampling PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

foreach(N 1 4)
  set(test_name "TestH5PartReaderV2StridedSampling-P${N}")
  set(test_added FALSE)
  if(N EQUAL 1)
    add_test(
      NAME ${test_name}
      COMMAND $<TARGET_FILE:TestH5PartReaderV2StridedSampling> -T
              "${PLUGIN_TEST_DIR}")
    set(test_added TRUE)
  else()
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${N})
      add_test(
        NAME ${test_name}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestH5PartReaderV2StridedSampling> -T
          "${PLUGIN_TEST_DIR}")
      set(test_added TRUE)
    endif()
  endif()
  if(test_added)
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:pvMeshlessFilters>:$ENV{LD_LIBRARY_PATH}")
  endif()
endforeach()

# --------------------------------------------------------------------------
# H5Part reader V2 CombineVectorComponents test
# --------------------------------------------------------------------------
add_executable(TestH5PartReaderV2CombineVectors
               TestH5PartReaderV2CombineVectors.cxx)
target_include_directories(
  TestH5PartReaderV2CombineVectors
  PRIVATE "${PROJECT_SOURCE_DIR}/Filters"
          "$<TARGET_PROPERTY:H5hut,INTERFACE_INCLUDE_DIRECTORIES>"
          "${HDF5_INCLUDE_DIRS}")
target_link_libraries(
  TestH5PartReaderV2CombineVectors
  PRIVATE pvMeshlessFilters VTK::TestingRendering MPI::MPI_CXX H5hut
          ${HDF5_LIBRARIES})
vtk_module_autoinit(TARGETS TestH5PartReaderV2CombineVectors MODULES
                    pvMeshlessFilters VTK::TestingRendering)
set_property(TARGET TestH5PartReaderV2CombineVectors APPEND PROPERTY BUILD_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2CombineVectors APPEND PROPERTY INSTALL_RPATH
  "$ORIGIN/../../lib/paraview-6.1/plugins/pv_meshless")
set_property(TARGET TestH5PartReaderV2CombineVectors PROPERTY BUILD_RPATH_USE_LINK_PATH TRUE)
set_property(TARGET TestH5PartReaderV2CombineVectors PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

foreach(N 1 4)
  set(test_name "TestH5PartReaderV2CombineVectors-P${N}")
  set(test_added FALSE)
  if(N EQUAL 1)
    add_test(
      NAME ${test_name}
      COMMAND $<TARGET_FILE:TestH5PartReaderV2CombineVectors> -T
              "${PLUGIN_TEST_DIR}")
    set(test_added TRUE)
  else()
    if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL ${N})
      add_test(
        NAME ${test_name}
        COMMAND
          ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestH5PartReaderV2CombineVectors> -T
          "${PLUGIN_TEST_DIR}")
      set(test_added TRUE)
    endif()
  endif()
  if(test_added)
    set_tests_properties(${test_name} PROPERTIES
      ENVIRONMENT "LD_LIBRARY_PATH=$<TARGET_FILE_DIR:pvMeshlessFilters>:$ENV{LD_LIBRARY_PATH}")
  endif()
endforeach()
