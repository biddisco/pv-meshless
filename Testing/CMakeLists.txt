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
  NAME TestSPHProbeKernel-P1
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
      NAME TestSPHProbeKernel-P${N}
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
        -value_range   "0.228203 1.1817"
        -peak_position "-0.163415 0.129112 0.00150253"
    )
  endforeach(N)
ENDIF (${Trilinos_FOUND})
    
#--------------------------------------------------
# Test Probe Filter : Shepard Mode
# Serial / Parallel
#--------------------------------------------------
ADD_TEST(
  NAME TestSPHProbeShepard-P1
  COMMAND 
    $<TARGET_FILE:TestSPHProbeFilter> 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F AGN_subset.h5part
    -T "${PLUGIN_TEST_DIR}"
    -neighbours 32
    -massScalars "mass"
    -scalar SmoothedDensity
    -value_range   "0.019823620095849 1565.79528808594"
    -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
)

IF (${Trilinos_FOUND})
  foreach(N 2;3;4)
    ADD_TEST(
      NAME TestSPHProbeShepard-P${N}
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
        -value_range   "0.019823620095849 1565.79528808594"
        -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
    )
  endforeach(N)
  
#  ADD_TEST(
#    NAME TestSPHProbeShepardLarge
#    COMMAND 
#      ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS}
#      $<TARGET_FILE:TestSPHProbeFilter> 
#      -D "${PROJECT_SOURCE_DIR}/Testing/data"
#      -F AGN_densitycheck.h5part
#      -T "${PLUGIN_TEST_DIR}"
#      -neighbours 32
#      -ghost_region 0.01
#      -massScalars "mass"
#      -scalar SmoothedDensity
#      -value_range   "0.0178490821272135 624893248"
#      -peak_position "0.557741165161133 0.581159114837646 0.552553355693817"
#  )
    
ENDIF (${Trilinos_FOUND})

#--------------------------------------------------
# Test Probe Filter : Image Resample + Kernel Mode
# jet Dataset : Serial/Parallel 
#--------------------------------------------------
SET(counter 1)
foreach(resolution 0.005 0.001 0.0005)
  foreach(N 1;2;3;4)
    SET(test_name "TestSPHImageResample-jet-N${counter}-R${resolution}-P${N}")
    #
    # Serial N==1
    #
    IF(N EQUAL 1)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F jet.h5part
          -T "${PLUGIN_TEST_DIR}"
          -particlesize 0.001
          -ghost_region 0.003
          -gridSpacing "${resolution} ${resolution} ${resolution}"
          -massScalars "mass"
          -scalar ShepardCoeff
          -contour 0.5
          -densityScalars "FLUID_density"
          -imagetest 1
          -cameraPosition "-0.160102 0.106656 -0.212755"
          -cameraFocus    "-0.160102 0.106656  0.023384"
          -cameraViewUp   "0 1 0"
          -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
      )
    ENDIF(N EQUAL 1)
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F jet.h5part
          -T "${PLUGIN_TEST_DIR}"
          -particlesize 0.001
          -ghost_region 0.005
          -gridSpacing "${resolution} ${resolution} ${resolution}"
          -scalar ShepardCoeff
          -contour 0.5
          -densityScalars "FLUID_density"
          -imagetest 1
          -cameraPosition "-0.160102 0.106656 -0.212755"
          -cameraFocus    "-0.160102 0.106656  0.023384"
          -cameraViewUp   "0 1 0"
          -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
      )
    ENDIF(${Trilinos_FOUND} AND N GREATER 1)
  endforeach(N)
  MATH(EXPR counter "${counter} + 1")
endforeach(resolution)
    
#--------------------------------------------------
# Test Probe Filter : Image Resample + Kernel Mode
# dam Dataset : Parallel only for high res
#--------------------------------------------------
SET(counter 1)
foreach(resolution 0.05 0.01 0.005)
  foreach(N 1 4 8 16 32 64)
    SET(test_name "TestSPHImageResample-dam-N${counter}-R${resolution}-P${N}")
    #
    # Serial N==1
    #
    IF((N EQUAL 1) AND (resolution EQUAL 0.05))
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F dam-17.h5part
          -T "${PLUGIN_TEST_DIR}"
          -particlesize 0.01
          -ghost_region 0.03
          -gridSpacing "${resolution} ${resolution} ${resolution}"
          -scalar ShepardCoeff
          -contour 0.1
          -densityScalars "FLUID_rho"
          -imagetest 1
          -windowSize "1024 768"
          -cameraPosition " 3.9250837 -0.4704097  0.5964158"
          -cameraFocus    "-2.0277536  2.7267457 -0.5096563"
          -cameraViewUp   "-0.1130509  0.1302781  0.9850112"
          -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
      )
    ENDIF((N EQUAL 1) AND (resolution EQUAL 0.05))
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F dam-17.h5part
          -T "${PLUGIN_TEST_DIR}"
          -particlesize 0.01
          -ghost_region 0.03
          -gridSpacing "${resolution} ${resolution} ${resolution}"
          -scalar ShepardCoeff
          -contour 0.1
          -densityScalars "FLUID_rho"
          -imagetest 1
          -windowSize "1024 768"
          -cameraPosition " 3.9250837 -0.4704097  0.5964158"
          -cameraFocus    "-2.0277536  2.7267457 -0.5096563"
          -cameraViewUp   "-0.1130509  0.1302781  0.9850112"
          -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
      )
    ENDIF(${Trilinos_FOUND} AND N GREATER 1)
  endforeach(N)
  MATH(EXPR counter "${counter} + 1")
endforeach(resolution)
    
