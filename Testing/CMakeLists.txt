INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/vtkCSCSMeshless/sph)
ADD_EXECUTABLE(TestKernel TestKernel.cxx)
TARGET_LINK_LIBRARIES(TestKernel
  pvmeshless_plugin 
)
ADD_TEST(TestKernel ${PROJECT_BINARY_DIR}/../bin/TestKernel)
