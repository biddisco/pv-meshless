#--------------------------------------------------
# Tests here are added using the 
# ADD_TEST(NAME ... COMMAND ...)
# form and will only run when ctest is run using
# ctest -C Release / ctest -C Debug / etc
#--------------------------------------------------

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
  ADD_TEST(TestH5PartParallelWriter ${PROJECT_BINARY_DIR}/../bin/TestH5PartParallelWriter)
ENDIF (VTK_USE_PARALLEL AND VTK_USE_MPI)

#--------------------------------------------------
# Trilinos
#--------------------------------------------------
IF (${Trilinos_FOUND})

  ADD_EXECUTABLE(TestParticlePartitionFilter TestParticlePartitionFilter.cxx)
  TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
  
  ADD_TEST(
    NAME TestParticlePartitionFilter 
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestParticlePartitionFilter> 
      -N 10000
      -F temp.h5
      -D ${PLUGIN_TEST_DIR}
      -T ${PLUGIN_TEST_DIR}
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/ParticlePartition.png
      -R
  )

  ADD_EXECUTABLE(TestSPHProbeParallel TestSPHProbeParallel.cxx)
  TARGET_LINK_LIBRARIES(TestSPHProbeParallel
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
    
  ADD_TEST(
    NAME TestSPHProbeSerial
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_subset.h5part
      -T "${PLUGIN_TEST_DIR}"
      -range "0.0310571286827326 372.910675048828"
      -index "693"
      -coord "0.492221683263779 0.506363093852997 0.486311554908752"
    )
  
  ADD_TEST(
    NAME TestSPHProbeParallel
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_densitycheck.h5part
      -T "${PLUGIN_TEST_DIR}"
      -range "0.0219669789075851 177405216"
      -index "625851"
      -coord "0.567801 0.586049 0.559156"
    )
  
ENDIF(${Trilinos_FOUND})

