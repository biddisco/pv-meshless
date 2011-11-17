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
ADD_EXECUTABLE(TestSPHProbeFilter TestSPHProbeFilter.cxx)
TARGET_LINK_LIBRARIES(TestSPHProbeFilter
  pv_meshless 
  ${Trilinos_LIBRARIES}
)
          
#--------------------------------------------------
# Test Probe Filter : Kernel Mode
# Serial / Parallel
#--------------------------------------------------
ADD_TEST(
  NAME TestSPHProbeKernel-small-serial
  COMMAND 
    $<TARGET_FILE:TestSPHProbeFilter> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F jet.h5part
    -T "${PLUGIN_TEST_DIR}"
    -particlesize 0.001
    -massScalars "mass"
    -scalar ShepardCoeff
    -value_range "0.228203 1.1817"
    -peak_position "-0.163415 0.129112 0.00150253"
)
  
IF (${Trilinos_FOUND})
  foreach(N 2;3;4)
    ADD_TEST(
      NAME TestSPHProbeKernel-small-parallel-${N}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
        $<TARGET_FILE:TestSPHProbeFilter> 
        -D "${PROJECT_SOURCE_DIR}/Testing/data"
        -F jet.h5part
        -T "${PLUGIN_TEST_DIR}"
        -ghost_region 0.002
        -particlesize 0.001
        -massScalars "mass"
        -scalar ShepardCoeff
        -value_range "0.228203 1.1817"
        -peak_position "-0.163415 0.129112 0.00150253"
    )
  endforeach(N)
ENDIF (${Trilinos_FOUND})
    
#--------------------------------------------------
# Test Probe Filter : Shepard Mode
# Serial / Parallel
#--------------------------------------------------
ADD_TEST(
  NAME TestSPHProbeShepard-small-serial
  COMMAND 
    $<TARGET_FILE:TestSPHProbeFilter> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F AGN_subset.h5part
    -T "${PLUGIN_TEST_DIR}"
    -neighbours 32
    -massScalars "mass"
    -scalar SmoothedDensity
    -value_range "0.019823620095849 1565.79528808594"
    -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
)

IF (${Trilinos_FOUND})
  foreach(N 2;3;4)
    ADD_TEST(
      NAME TestSPHProbeShepard-small-parallel-${N}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
        $<TARGET_FILE:TestSPHProbeFilter> 
        -D "${PROJECT_SOURCE_DIR}/Testing/data"
        -F AGN_subset.h5part
        -T "${PLUGIN_TEST_DIR}"
        -neighbours 32
        -ghost_region 0.01
        -massScalars "mass"
        -scalar SmoothedDensity
        -value_range "0.019823620095849 1565.79528808594"
        -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
    )
  endforeach(N)
  
  ADD_TEST(
    NAME TestSPHProbeShepardLarge
    COMMAND 
      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
      $<TARGET_FILE:TestSPHProbeFilter> 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F AGN_densitycheck.h5part
      -T "${PLUGIN_TEST_DIR}"
      -neighbours 32
      -ghost_region 0.01
      -massScalars "mass"
      -scalar SmoothedDensity
      -value_range "0.0178490821272135 624893248"
      -peak_position "0.557741165161133 0.581159114837646 0.552553355693817"
  )
    
ENDIF (${Trilinos_FOUND})

#--------------------------------------------------
# Test Probe Filter : Image Resample + Kernel Mode
# Serial / Parallel
#--------------------------------------------------
ADD_TEST(
  NAME TestSPHImageResample-small-serial
  COMMAND 
    $<TARGET_FILE:TestSPHProbeFilter> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F jet.h5part
    -T "${PLUGIN_TEST_DIR}"
    -particlesize 0.001
    -ghost_region 0.005
    -massScalars "mass"
    -scalar ShepardCoeff
    -contour 0.5
    -densityScalars "FLUID_density"
    -imagetest 1
    -cameraPosition "-0.160102 0.106656 -0.212755"
    -cameraFocus    "-0.160102 0.106656  0.023384"
    -cameraViewUp   "0 1 0"
    -V ${PROJECT_SOURCE_DIR}/Testing/baseline/TestSPHImageResample-small-serial.png
)

IF (${Trilinos_FOUND})
  foreach(N 2;3;4)
    ADD_TEST(
      NAME TestSPHImageResample-small-parallel-${N}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
        $<TARGET_FILE:TestSPHProbeFilter> 
        -D "${PROJECT_SOURCE_DIR}/Testing/data"
        -F jet.h5part
        -T "${PLUGIN_TEST_DIR}"
        -particlesize 0.001
        -ghost_region 0.005
        -scalar ShepardCoeff
        -contour 0.5
        -densityScalars "FLUID_density"
        -imagetest 1
        -cameraPosition "-0.160102 0.106656 -0.212755"
        -cameraFocus    "-0.160102 0.106656  0.023384"
        -cameraViewUp   "0 1 0"
        -V ${PROJECT_SOURCE_DIR}/Testing/baseline/TestSPHImageResample-small-parallel-${N}.png
    )
  endforeach(N)
ENDIF (${Trilinos_FOUND})
    
IF (${Trilinos_FOUND})
  foreach(resolution 0.05 0.01 0.005 0.0025)
    foreach(N 4 8 16 32 64)
      ADD_TEST(
        NAME TestSPHImageResample-large-parallel-${resolution}-${N}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F dam-17.h5part
          -T "${PLUGIN_TEST_DIR}"
          -particlesize 0.005
          -ghost_region 0.01
          -gridSpacing "${resolution} ${resolution} ${resolution}"
          -scalar ShepardCoeff
          -contour 0.1
          -densityScalars "FLUID_rho"
          -imagetest 1
          -windowSize "1024 768"
          -cameraPosition "3.92508375929083 -0.470409761751235 0.596415817759249"
          -cameraFocus    "-2.02775364406982 2.726745716382 -0.509656384551184"
          -cameraViewUp   "-0.11305097024463 0.130278141647241 0.98501121005585"
          -V ${PROJECT_SOURCE_DIR}/Testing/baseline/TestSPHImageResample-large-parallel-${resolution}-${N}.png
      )
    endforeach(N)
  endforeach(resolution)
ENDIF (${Trilinos_FOUND})
    
