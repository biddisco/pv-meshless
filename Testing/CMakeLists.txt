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
# convenience macro for a>=b for testing processor count
#--------------------------------------------------
MACRO(greater_equal a b boolresult) 
  IF((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 1)
  ELSE((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 0)
  ENDIF((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
ENDMACRO(greater_equal)

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
IF(VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_EXECUTABLE(TestH5PartParallelWriter TestH5PartParallelWriter.cxx)
  TARGET_LINK_LIBRARIES(TestH5PartParallelWriter
    pv_meshless 
    )
  ADD_TEST(TestH5PartParallelWriter ${PROJECT_BINARY_DIR}/../bin/TestH5PartParallelWriter)
ENDIF (VTK_USE_PARALLEL AND VTK_USE_MPI)

#--------------------------------------------------
# Slurm job manager require "-rmk" "slurm" arguments
# so split MPIEXEC_NUMPROC_FLAG into a list of args
# to stop cmake from putting quotes around all of it
#--------------------------------------------------
string(REPLACE " " ";" MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG})

#--------------------------------------------------
# Trilinos
#--------------------------------------------------
IF(${Trilinos_FOUND})

  ADD_EXECUTABLE(TestParticlePartitionFilter 
    TestParticlePartitionFilter.cxx
  )
  TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
    pv_meshless 
    ${Trilinos_LIBRARIES}
    )
  
  greater_equal(${MPIEXEC_MAX_NUMPROCS} 4 processors)
  IF(processors)
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
  ENDIF(processors)
 
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
FOREACH(N 1 2 3 4)
  SET(test_name "TestSPHProbeKernel-P${N}")
  SET(TestParams 
    -D "${PROJECT_SOURCE_DIR}/Testing/data"
    -F jet.h5part
    -T "${PLUGIN_TEST_DIR}"
    -testName ${test_name}
    -particlesize 0.001
    -ghost_region 0.003
    -massScalars "mass"
    -scalar ShepardCoeff
    -value_range "0.228203 1.1817"
    -peak_position "-0.163415 0.129112 0.00150253"
  )
  #
  # Serial N==1
  #
  IF(N EQUAL 1)
    ADD_TEST(
      NAME ${test_name}
      COMMAND 
        $<TARGET_FILE:TestSPHProbeFilter> 
        ${TestParams}
    )
  ENDIF(N EQUAL 1)
  #
  # Parallel N>1
  #
  IF(${Trilinos_FOUND} AND N GREATER 1)
    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF(processors)
      ADD_TEST(
      NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF(processors)
  ENDIF (${Trilinos_FOUND} AND N GREATER 1)
ENDFOREACH(N)
    
#--------------------------------------------------
# Test Probe Filter : Shepard Mode
# Serial / Parallel
#--------------------------------------------------
FOREACH(N 1 2 3 4)
  SET(test_name "TestSPHProbeShepard-subset-P${N}")
    #
    # Serial N==1
    #
    IF(N EQUAL 1)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F AGN_subset.h5part
          -T "${PLUGIN_TEST_DIR}"
          -testName ${test_name}
          -neighbours 32
          -ghost_region 0.01
          -massScalars "mass"
          -scalar SmoothedDensity
          -value_range   "0.019823620095849 1565.79528808594"
          -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
      )
    ENDIF(N EQUAL 1)
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
      IF(processors)
        ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          -D "${PROJECT_SOURCE_DIR}/Testing/data"
          -F AGN_subset.h5part
          -T "${PLUGIN_TEST_DIR}"
          -testName ${test_name}
          -neighbours 32
          -ghost_region 0.01
          -massScalars "mass"
          -scalar SmoothedDensity
          -value_range   "0.019823620095849 1565.79528808594"
          -peak_position "0.53753125667572 0.529183149337769 0.526147127151489"
      )
      ENDIF(processors)
  ENDIF (${Trilinos_FOUND} AND N GREATER 1)
ENDFOREACH(N)
  
FOREACH(N 16 32 64 128 192 256 276)
  SET(test_name "TestSPHProbeShepard-densitycheck-P${N}")
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
      IF(processors)
        ADD_TEST(
          NAME ${test_name}
          COMMAND 
            ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
            $<TARGET_FILE:TestSPHProbeFilter> 
            -D "${PROJECT_SOURCE_DIR}/Testing/data"
            -F AGN_densitycheck.h5part
            -T "${PLUGIN_TEST_DIR}"
            -testName ${test_name}
            -neighbours 32
            -ghost_region 0.01
            -massScalars "mass"
            -scalar SmoothedDensity
            -value_range   "0.0178490821272135 624893248"
            -peak_position "0.557741165161133 0.581159114837646 0.552553355693817"
        )
      ENDIF(processors)
  ENDIF (${Trilinos_FOUND} AND N GREATER 1)
ENDFOREACH(N)

#--------------------------------------------------
# Test Probe Filter : Image Resample + Kernel Mode
# jet Dataset : Serial/Parallel 
#--------------------------------------------------
SET(counter 1)
FOREACH(resolution 0.005 0.001 0.0005)
  FOREACH(N 1 2 3 4)
    SET(test_name "TestSPHImageResample-jet-N${counter}-R${resolution}-P${N}")
    SET(TestParams 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F jet.h5part
      -T "${PLUGIN_TEST_DIR}"
      -testName ${test_name}
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
      -imageScalars "ProcessId"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )
    #
    # Serial N==1
    #
    IF(N EQUAL 1)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF(N EQUAL 1)
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
      IF(processors)
        ADD_TEST(
          NAME ${test_name}
          COMMAND 
            ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
            $<TARGET_FILE:TestSPHProbeFilter> 
            ${TestParams}
        )
      ENDIF(processors)
    ENDIF(${Trilinos_FOUND} AND N GREATER 1)
  ENDFOREACH(N)
  MATH(EXPR counter "${counter} + 1")
ENDFOREACH(resolution)
    
#--------------------------------------------------
# Test Probe Filter : Image Resample + Kernel Mode
# dam Dataset : Parallel only for high res
#--------------------------------------------------
SET(counter 1)
FOREACH(resolution 0.05 0.01 0.005)
  FOREACH(N 1 4 8 16 32 64 128 192 256 276)
    SET(test_name "TestSPHImageResample-dam-N${counter}-R${resolution}-P${N}")
    SET(TestParams 
      -D "${PROJECT_SOURCE_DIR}/Testing/data"
      -F dam-17.h5part
      -T "${PLUGIN_TEST_DIR}"
      -testName ${test_name}
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
      -imageScalars "ProcessId"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )
    #
    # Serial N==1
    #
    IF((N EQUAL 1) AND (resolution EQUAL 0.05))
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF((N EQUAL 1) AND (resolution EQUAL 0.05))
    #
    # Parallel N>1
    #
    IF(${Trilinos_FOUND} AND N GREATER 1)
      greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
      IF(processors)
        ADD_TEST(
          NAME ${test_name}
          COMMAND 
            ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
            $<TARGET_FILE:TestSPHProbeFilter> 
            ${TestParams}
        )
      ENDIF(processors)
    ENDIF(${Trilinos_FOUND} AND N GREATER 1)
  ENDFOREACH(N)
  MATH(EXPR counter "${counter} + 1")
ENDFOREACH(resolution)
    
