#--------------------------------------------------
# Tests here are added using the 
# ADD_TEST(NAME ... COMMAND ...)
# form and will only run on windows when ctest is run using
# ctest -C Release / Debug / etc
#--------------------------------------------------

#--------------------------------------------------
# Setup some convenience flags depending on machine
#--------------------------------------------------
SET(COMPUTERNAME $ENV{COMPUTERNAME})
IF (COMPUTERNAME)
  SET(HOSTNAME $ENV{COMPUTERNAME})
ELSE(COMPUTERNAME)
  SET(HOSTNAME $ENV{HOSTNAME})
ENDIF (COMPUTERNAME)

#--------------------------------------------------
# Machines running OSMesa instead of OpenGL might need a larger threshold for image test comparisons
#--------------------------------------------------
set(IMAGE_THRESHOLD "0")

#--------------------------------------------------
# Some data paths which are machine dependent
#--------------------------------------------------
set(CORES_PER_NODE 1)
IF (HOSTNAME MATCHES "DINO")
  set(CORES_PER_NODE 2)
  set(LARGE_DATA_PATH "C:\\data\\nextmuse")
  set(SC12_DATA_PATH "C:/Data/sc12")
  set(SC12_DATA_FILE "wave_tank7_")
  set(SC12_DATA_NMIN 2)
  set(SC12_DATA_BOXN 0004)
  set(SC12_DATA_DX_N 0.5)
  set(SC12_DATA_DX_2 1.0)
  set(SC12_DATA_DX_4 2.0)
ELSEIF (HOSTNAME MATCHES "BRENO")
ELSEIF (HOSTNAME MATCHES "todi")
  set(CORES_PER_NODE 16)
  set(LARGE_DATA_PATH "/scratch/todi/biddisco/nextmuse")
  set(SC12_DATA_PATH "/scratch/todi/biddisco/nextmuse")
  set(SC12_DATA_FILE "wave_tank1_")
  set(SC12_DATA_NMIN 4)
  set(SC12_DATA_BOXN 00032    00064     00128    00256     00512     01024     02048     04096)
  set(SC12_DATA_DX_N 0.132832 0.105429 0.0836791 0.0664162 0.0527145 0.0418396 0.0332081 0.0263573)
  set(SC12_DATA_DX_2 0.265664 0.210858 0.1673582 0.1328324 0.1054290 0.0836792 0.0664162 0.0527146)
  set(SC12_DATA_DX_4 0.531328 0.421716 0.3347164 0.2656648 0.2108580 0.1673584 0.1328324 0.1054292)
ELSEIF (HOSTNAME MATCHES "rosa")
  set(CORES_PER_NODE 16)
  set(SC12_DATA_PATH "/scratch/weisshorn/biddisco")
  set(SC12_DATA_FILE "wave_tank1_")
  set(SC12_DATA_NMIN 4)
  set(SC12_DATA_BOXN 00032    00064     00128    00256     00512     01024     02048     04096)
  set(SC12_DATA_DX_N 0.132832 0.105429 0.0836791 0.0664162 0.0527145 0.0418396 0.0332081 0.0263573)
  set(SC12_DATA_DX_2 0.265664 0.210858 0.1673582 0.1328324 0.1054290 0.0836792 0.0664162 0.0527146)
  set(SC12_DATA_DX_4 0.531328 0.421716 0.3347164 0.2656648 0.2108580 0.1673584 0.1328324 0.1054292)
ELSEIF (HOSTNAME MATCHES "eiger")
  set(LARGE_DATA_PATH "/project/csvis/biddisco/nextmuse")
ENDIF ()
# MESSAGE("Host name is ${HOSTNAME} and large data path is ${LARGE_DATA_PATH}")

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
  IF ((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 1)
  ELSE((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
    SET(${boolresult} 0)
  ENDIF ((${a} GREATER ${b}) OR (${a} EQUAL ${b}))
ENDMACRO(greater_equal)

#--------------------------------------------------
# Simple non MPI Filter/Kernel Tests 
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
# Test H5PartWriter
#--------------------------------------------------
IF (USE_HDF5 AND VTK_USE_MPI)
  ADD_EXECUTABLE(TestH5PartParallelWriter TestH5PartParallelWriter.cxx)
  TARGET_LINK_LIBRARIES(TestH5PartParallelWriter
    pv_meshless 
    )
  FOREACH(N 1 2 3 4 8 16 32 64)
    SET(test_name "TestH5PartParallelWriter-P${N}")
      #
      # Serial N==1
      #
      IF (N EQUAL 1)
        ADD_TEST(
          NAME ${test_name}
          COMMAND 
            $<TARGET_FILE:TestH5PartParallelWriter> 
            -N 1000000
            -F temp.h5
            -D ${PLUGIN_TEST_DIR}
            -T "${PLUGIN_TEST_DIR}"
            -testName ${test_name}
        )
      ENDIF (N EQUAL 1)
      #
      # Parallel N>1
      #
      IF (N GREATER 1)
        greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
        IF (processors)
          ADD_TEST(
          NAME ${test_name}
          COMMAND 
            ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
            $<TARGET_FILE:TestH5PartParallelWriter> 
            -N 1000000
            -F temp.h5
            -D ${PLUGIN_TEST_DIR}
            -T "${PLUGIN_TEST_DIR}"
            -testName ${test_name}
        )
        ENDIF (processors)
    ENDIF (N GREATER 1)
  ENDFOREACH(N)
  
ENDIF (USE_HDF5 AND VTK_USE_MPI)

#--------------------------------------------------
# Slurm job manager require "-rmk" "slurm" arguments
# so split MPIEXEC_NUMPROC_FLAG into a list of args
# to stop cmake from putting quotes around all of it
#--------------------------------------------------
string(REPLACE " " ";" MPIEXEC_NUMPROC_FLAG ${MPIEXEC_NUMPROC_FLAG})

#--------------------------------------------------
# Trilinos
#--------------------------------------------------
IF (${TRILINOS_ENABLED} AND VTK_USE_MPI)

  ADD_LIBRARY(pvmeshlessTestLib STATIC
    TestUtils.cxx
  )

  ADD_EXECUTABLE(TestParticlePartitionFilter 
    TestParticlePartitionFilter.cxx
  )
  TARGET_LINK_LIBRARIES(TestParticlePartitionFilter
    pv_meshless 
    pvmeshlessTestLib
    ${Trilinos_LIBRARIES}
    )
  
  greater_equal(${MPIEXEC_MAX_NUMPROCS} 4 processors)
  SET(test_name "TestParticlePartitionFilter-P4")
  IF (processors)
    ADD_TEST(
      NAME ${test_name}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} 4
        $<TARGET_FILE:TestParticlePartitionFilter> 
        -testName ${test_name}
        -generateParticles 1000
        -V ${PROJECT_SOURCE_DIR}/Testing/baseline/ParticlePartition.png
        -doRender 1
    )
  ENDIF (processors)
 
ENDIF (${TRILINOS_ENABLED} AND VTK_USE_MPI)

#--------------------------------------------------
# All tests from here on make use of HDF5 for IO
#--------------------------------------------------
IF (NOT USE_HDF5) 
  return()
ENDIF (NOT USE_HDF5) 

#--------------------------------------------------
# All tests from here on make use of MPI
#--------------------------------------------------
IF (NOT VTK_USE_MPI) 
  return()
ENDIF (NOT VTK_USE_MPI) 

#--------------------------------------------------
# Test Probe Filter
#--------------------------------------------------
ADD_EXECUTABLE(TestSPHProbeFilter TestSPHProbeFilter.cxx)
TARGET_LINK_LIBRARIES(TestSPHProbeFilter
  pv_meshless 
  pvmeshlessTestLib
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
  greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
  IF (processors)
    ADD_TEST(
    NAME ${test_name}
      COMMAND 
        ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
        $<TARGET_FILE:TestSPHProbeFilter> 
        ${TestParams}
    )
  ENDIF (processors)
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
    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF (processors)
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
    ENDIF (processors)
ENDFOREACH(N)
  
FOREACH(N 16 32 64 128 192 256 276)
  SET(test_name "TestSPHProbeShepard-densitycheck-P${N}")
  #
  # Parallel N>1
  #
  greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
  IF (processors)
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
  ENDIF (processors)
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
      -scalar ShepardCoeff
      -contour 0.5
      -densityScalars "FLUID_density"
      -imageResample 1
      -cameraPosition "-0.160102 0.106656 -0.212755"
      -cameraFocus    "-0.160102 0.106656  0.023384"
      -cameraViewUp   "0 1 0"
      -imageScalars "ProcessId"
      -imageThreshold "${IMAGE_THRESHOLD}"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )
    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF (processors)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF (processors)
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
      -imageResample 1
      -windowSize "1024 768"
      -cameraPosition " 3.9250837 -0.4704097  0.5964158"
      -cameraFocus    "-2.0277536  2.7267457 -0.5096563"
      -cameraViewUp   "-0.1130509  0.1302781  0.9850112"
      -imageScalars "ProcessId"
      -imageThreshold "${IMAGE_THRESHOLD}"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )
    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF (processors)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF (processors)
  ENDFOREACH(N)
  MATH(EXPR counter "${counter} + 1")
ENDFOREACH(resolution)
    
#--------------------------------------------------
# Test Probe Filter for very large parallel data
# Image Resample + Kernel Mode
# dam Dataset(s) : Parallel only 
#--------------------------------------------------
SET(NUMPARTICLES    
  3234168 6449040 12873672 25672080 51371840 102832380 
  205376640 410490080 820490990 1642018560 3280404204
)
SET(SIM_CPUS 
  00032 00064 00128 00256 00512 01024
  02048 04096 08192 16384 32768
)
SET(RESOLUTIONS
    5.95392956532109E-003 4.72563702766763E-003 3.75074059446973E-003
    2.97696478266054E-003 2.36281851383381E-003 1.87537029723487E-003
    1.48848239133027E-003 1.18140925691691E-003 9.37685148617434E-004
    7.44241195665136E-004 5.90704628458454E-004
)
SET(GHOSTWIDTH
    0.018 0.015 0.012 0.009 0.007 0.006 0.0045 0.0036 0.0030 0.0023 0.002
)
SET(SAMPLERES
    0.02 0.015 0.012 0.009 0.007 0.006 0.0045 0.0035 0.0030 0.0025 0.002
)

#--------------------------------------------------
# Test 1 : Take one dataset and contour it to the same resolution
# on several different CPU counts to test strong scaling
# 1024^3 
#--------------------------------------------------
SET(counter 1)
FOREACH(index 5)
  LIST(GET NUMPARTICLES ${index} NUM)
  LIST(GET SIM_CPUS     ${index} CPU)
  LIST(GET RESOLUTIONS  ${index} RES)
  LIST(GET GHOSTWIDTH   ${index} GST)
  SET(GRIDSIZE 1024)
  SET(file_name "${CPU}_${RES}_dambreak0.h5part")
  FOREACH(N 32 64 128 192 256 384 512 768 1024 1536 2048)
    SET(test_name "benchmark-strong-N${counter}-R${index}-P${N}")
#    MESSAGE("${test_name} ${file_name} ${GST}")
    SET(TestParams 
      -D "${LARGE_DATA_PATH}"
      -F ${file_name}
      -T "${PLUGIN_TEST_DIR}"
      -testName ${test_name}
      -Xarray "X"
      -Yarray "Y"
      -Zarray "Z"
      -particlesize ${RES}
      -ghost_region ${GST}
      -gridResolution "${GRIDSIZE} ${GRIDSIZE} ${GRIDSIZE}"
      -scalar ShepardCoeff
      -contour 0.8
      -imageResample 1
      -skipImageTest 1
      -windowSize "1024 768"
      -cameraPosition " 1.48858 2.90908 1.08421"
      -cameraFocus    " 2.606   0       0.275"
      -cameraViewUp   "0.0790531 -0.238854 0.967832"
      -imageScalars "ProcessId"
      -imageThreshold "${IMAGE_THRESHOLD}"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )

    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF (processors)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF (processors)
  ENDFOREACH(N)
  MATH(EXPR counter "${counter} + 1")
ENDFOREACH(index)

#--------------------------------------------------
# Test 2 : Take one dataset and contour it to the same resolution
# on several different CPU counts to test strong scaling
#--------------------------------------------------
SET(counter 1)
FOREACH(index 0 1 2 3 4 5 6 7 8 9 10)
  LIST(GET NUMPARTICLES ${index} NUM)
  LIST(GET SIM_CPUS     ${index} CPU)
  LIST(GET RESOLUTIONS  ${index} RES)
  LIST(GET GHOSTWIDTH   ${index} GST)
  SET(GRIDSIZE 1024)
  SET(file_name "${CPU}_${RES}_dambreak0.h5part")
  FOREACH(N 2048)
    SET(test_name "benchmark-weak-N${counter}-R${index}-P${N}")
#    MESSAGE("${test_name} ${file_name} ${GST}")
    SET(TestParams 
      -D "${LARGE_DATA_PATH}"
      -F ${file_name}
      -T "${PLUGIN_TEST_DIR}"
      -testName ${test_name}
      -Xarray "X"
      -Yarray "Y"
      -Zarray "Z"
      -expectedparticles ${NUM}
      -particlesize ${RES}
      -ghost_region ${GST}
      -gridResolution "${GRIDSIZE} ${GRIDSIZE} ${GRIDSIZE}"
      -scalar ShepardCoeff
      -contour 0.8
      -imageResample 1
      -skipImageTest 1
      -windowSize "1024 768"
      -cameraPosition " 1.48858 2.90908 1.08421"
      -cameraFocus    " 2.606   0       0.275"
      -cameraViewUp   "0.0790531 -0.238854 0.967832"
      -imageScalars "ProcessId"
      -imageThreshold "${IMAGE_THRESHOLD}"
      -V ${PROJECT_SOURCE_DIR}/Testing/baseline/${test_name}.png
    )

    greater_equal(${MPIEXEC_MAX_NUMPROCS} ${N} processors)
    IF (processors)
      ADD_TEST(
        NAME ${test_name}
        COMMAND 
          ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${N}
          $<TARGET_FILE:TestSPHProbeFilter> 
          ${TestParams}
      )
    ENDIF (processors)
  ENDFOREACH(N)
  MATH(EXPR counter "${counter} + 1")
ENDFOREACH(index)

#/project/csvis/biddisco/nextmuse/00032_5.95392956532109E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/00064_4.72563702766763E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/00128_3.75074059446973E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/00256_2.97696478266054E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/00512_2.36281851383381E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/01024_1.87537029723487E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/02048_1.48848239133027E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/04096_1.18140925691691E-003_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/08192_9.37685148617434E-004_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/16384_7.44241195665136E-004_dambreak0.h5part
#/project/csvis/biddisco/nextmuse/32768_5.90704628458454E-004_dambreak0.h5part

#--------------------------------------------------
#--------------------------------------------------
# Benchmark SC12 : Take one dataset with bounding boxes
# saved and run partitioning filter on different CPU counts
# with/without boxes. Time the partitioning.
#--------------------------------------------------
#--------------------------------------------------
LIST(GET SC12_DATA_BOXN 0 SC12_TEST)
IF (EXISTS ${SC12_DATA_PATH}/${SC12_DATA_FILE}${SC12_TEST}.h5part)
  FOREACH(index 0 1)
    SET(ncount 0)
    FOREACH(N ${SC12_DATA_BOXN})
      LIST(GET SC12_DATA_DX_N ${ncount} DELTAX)
      LIST(GET SC12_DATA_DX_4 ${ncount} GST)
      SET(TestParams 
        -D "${SC12_DATA_PATH}"
        -F "${SC12_DATA_FILE}${N}.h5part"
        -T "${PLUGIN_TEST_DIR}"
        -Xarray "X"
        -Yarray "Y"
        -Zarray "Z"
        -particlesize ${DELTAX}
        -ghost_region ${GST}
        -ignorePartitions "${index}"
      )
      # remove 000 from 00032 for Nprocs etc
      STRING(REGEX REPLACE "(^0+)" "" N "${N}")
      #
      SET(M ${N})
      while (NOT ${M} LESS ${SC12_DATA_NMIN})
        IF (NOT ${M} GREATER ${MPIEXEC_MAX_NUMPROCS})
          SET(test_name "benchmark-PartitionBoxes-B${N}-P${M}-I${index}")
          MATH(EXPR MPI_DEPTH "(${N}/${CORES_PER_NODE})/${M}")
          IF (${MPI_DEPTH} GREATER 1)
            SET (MPI_DEPTH_FLAGS -d ${MPI_DEPTH})
          ELSE (${MPI_DEPTH} GREATER 1)
            SET (MPI_DEPTH_FLAGS )
          ENDIF (${MPI_DEPTH} GREATER 1)
          IF (${MPI_DEPTH} GREATER ${CORES_PER_NODE})
            SET (MPI_DEPTH_FLAGS -d ${CORES_PER_NODE})
          ENDIF (${MPI_DEPTH} GREATER ${CORES_PER_NODE})
          ADD_TEST(
            NAME ${test_name}
            COMMAND 
              ${MPIEXEC} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${M} ${MPI_DEPTH_FLAGS}
              $<TARGET_FILE:TestParticlePartitionFilter> 
              -testName ${test_name}
              ${TestParams}
          )
        ENDIF (NOT ${M} GREATER ${MPIEXEC_MAX_NUMPROCS})
        #
        MATH(EXPR M "${M} / 2")
      endwhile (NOT ${M} LESS ${SC12_DATA_NMIN})
      MATH(EXPR ncount "${ncount} + 1")
    ENDFOREACH(N)
  ENDFOREACH(index)
ENDIF (EXISTS ${SC12_DATA_PATH}/${SC12_DATA_FILE}${SC12_TEST}.h5part)
