INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/vtkCSCSMeshless/sph)

#--------------------------------------------------
# a macro to add one test
#--------------------------------------------------
MACRO(add_local_test test_name) 
  ADD_EXECUTABLE(${test_name} ${test_name}.cxx)
  TARGET_LINK_LIBRARIES(${test_name}
    pv_meshless 
    )
  ADD_TEST(${test_name} ${PROJECT_BINARY_DIR}/../bin/${test_name})
ENDMACRO(add_local_test)

#--------------------------------------------------
# Cxx Tests 
#--------------------------------------------------
SET(TEST_LIST 
  TestKernel3D
  TestKernel2D 
  TestCenterOfMassFilter
  TestMomentsOfInertiaFilter
)

FOREACH(test_name ${TEST_LIST})
  add_local_test(${test_name})
ENDFOREACH(test_name)

#--------------------------------------------------
# H5Part
#--------------------------------------------------
IF (VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_EXECUTABLE(TestH5PartParallelWriter TestH5PartParallelWriter.cxx)
  TARGET_LINK_LIBRARIES(TestH5PartParallelWriter
    pv_meshless 
    )
  ADD_TEST(TestH5PartParallelWriter.cxx ${PROJECT_BINARY_DIR}/../bin/TestH5PartParallelWriter.cxx)
ENDIF (VTK_USE_PARALLEL AND VTK_USE_MPI)

#--------------------------------------------------
# Trilinos
#--------------------------------------------------
IF(${Trilinos_FOUND})

  ADD_EXECUTABLE(TestParticlePartitionFilter TestParticlePartitionFilter.cxx)
  TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
  ADD_TEST(TestParticlePartitionFilter.cxx ${PROJECT_BINARY_DIR}/../bin/TestParticlePartitionFilter)

  ADD_EXECUTABLE(TestSPHProbeParallel TestSPHProbeParallel.cxx)
  TARGET_LINK_LIBRARIES(TestSPHProbeParallel
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
  ADD_TEST(TestSPHProbeParallel.cxx ${PROJECT_BINARY_DIR}/../bin/TestSPHProbeParallel)
  
ENDIF(${Trilinos_FOUND})

