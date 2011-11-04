#--------------------------------------------------
# Tests here are added using the 
# ADD_TEST(NAME ... COMMAND ...)
# form and will only run on windows when ctest is run using
# ctest -C Release / Debug / etc
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

  ADD_DEFINITIONS(-DPV_MESHLESS_TRILINOS)

  ADD_EXECUTABLE(TestParticlePartitionFilter 
    TestParticlePartitionFilter.cxx
#    R250.cxx
  )
  TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
  
  ADD_TEST(
    NAME TestParticlePartitionFilter 
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestParticlePartitionFilter> 
      -N 1000
      -F temp.h5
      -D ${PLUGIN_TEST_DIR}
      -T ${PLUGIN_TEST_DIR}
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/ParticlePartition.png
      -R
  )
 
ENDIF(${Trilinos_FOUND})

#--------------------------------------------------
# Test Probe Filter
#--------------------------------------------------
ADD_EXECUTABLE(TestSPHProbeParallel TestSPHProbeParallel.cxx)
TARGET_LINK_LIBRARIES(TestSPHProbeParallel
  pv_meshless 
  ${Trilinos_LIBRARIES}
  )
    
ADD_TEST(
  NAME TestSPHProbeKernelSmallSerial
  COMMAND 
    $<TARGET_FILE:TestSPHProbeParallel> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F jet.h5part
    -T "${PLUGIN_TEST_DIR}"
    -particlesize 0.001
    -scalar ShepardCoeff
    -value_range "0.228203 1.1817"
    -peak_position "-0.163415 0.129112 0.00150253"
  )
  
ADD_TEST(
  NAME TestSPHProbeShepardSmallSerial
  COMMAND 
    $<TARGET_FILE:TestSPHProbeParallel> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F AGN_subset.h5part
    -T "${PLUGIN_TEST_DIR}"
    -neighbours 32
    -scalar SmoothedDensity
    -value_range "0.019823620095849 1565.79528808594"
    -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
  )
  
#--------------------------------------------------
# Test Probe Filter : Parallel
#--------------------------------------------------
IF (${Trilinos_FOUND})
  ADD_TEST(
    NAME TestSPHProbeShepardSmallParallel2
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_subset.h5part
      -T "${PLUGIN_TEST_DIR}"
      -neighbours 32
      -ghost_region 0.01
      -scalar SmoothedDensity
      -value_range "0.019823620095849 1565.79528808594"
      -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
    )
  
  ADD_TEST(
    NAME TestSPHProbeShepardSmallParallel4
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_subset.h5part
      -T "${PLUGIN_TEST_DIR}"
      -neighbours 32
      -ghost_region 0.01
      -scalar SmoothedDensity
      -value_range "0.019823620095849 1565.79528808594"
      -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
    )
  
  ADD_TEST(
    NAME TestSPHProbeKernelSmallParallel2
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 2
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F jet.h5part
      -T "${PLUGIN_TEST_DIR}"
      -ghost_region 0.002
      -particlesize 0.001
      -scalar ShepardCoeff
      -value_range "0.228203 1.1817"
      -peak_position "-0.163415 0.129112 0.00150253"
    )
  
  ADD_TEST(
    NAME TestSPHProbeKernelSmallParallel4
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F jet.h5part
      -T "${PLUGIN_TEST_DIR}"
      -ghost_region 0.002
      -particlesize 0.001
      -scalar ShepardCoeff
      -value_range "0.228203 1.1817"
      -peak_position "-0.163415 0.129112 0.00150253"
    )
  
  ADD_TEST(
    NAME TestSPHProbeShepardLarge
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestSPHProbeParallel> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_densitycheck.h5part
      -T "${PLUGIN_TEST_DIR}"
      -neighbours 32
      -ghost_region 0.01
      -scalar SmoothedDensity
      -value_range "0.0178490821272135 624893248"
      -peak_position "0.557741165161133 0.581159114837646 0.552553355693817"
    )
    
ENDIF (${Trilinos_FOUND})
