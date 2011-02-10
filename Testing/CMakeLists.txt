INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/vtkCSCSMeshless/sph)

# a list of all tests
SET(TEST_LIST 
  TestKernel3D
  TestKernel2D 
  TestCenerOfMassFilter
  TestMomentsOfInertiaFilter
)

# a macro to add one test
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

LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS}/${CMAKE_CFG_INTDIR})

MESSAGE("${Trilinos_LIBRARY_DIRS}/${CMAKE_CFG_INTDIR}")

ADD_EXECUTABLE(TestH5PartParallelWriter TestH5PartParallelWriter.cxx)
TARGET_LINK_LIBRARIES(TestH5PartParallelWriter
  pvmeshless_plugin 
  )
ADD_TEST(TestH5PartParallelWriter.cxx ${PROJECT_BINARY_DIR}/../bin/TestH5PartParallelWriter.cxx)


ADD_EXECUTABLE(TestParticlePartition TestParticlePartition.cxx)
TARGET_LINK_LIBRARIES(TestParticlePartition
  pvmeshless_plugin 
  ${Trilinos_LIBRARIES}
  )
ADD_TEST(TestParticlePartition.cxx ${PROJECT_BINARY_DIR}/../bin/TestParticlePartition.cxx)
