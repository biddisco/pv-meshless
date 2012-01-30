PROJECT("pv-meshless-Tools")

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#--------------------------------------------------------
# Setup include paths
#--------------------------------------------------------
INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CSCS_COMMON_INCLUDE_DIRS}
)
 
SET(TOOL_VTK_LIBS
  pv_common
  pv_meshless
  pv_astro
  vtkWidgets
  vtkHybrid
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkFiltering
  vtkCommon
  vtksys
)

IF(VTK_USE_MPI)
  SET(TOOL_VTK_LIBS ${TOOL_VTK_LIBS} ${MPI_LIBRARY} )
  IF(MPI_EXTRA_LIBRARY)
    SET(TOOL_VTK_LIBS ${TOOL_VTK_LIBS} ${MPI_EXTRA_LIBRARY} )
  ENDIF(MPI_EXTRA_LIBRARY)
ENDIF(VTK_USE_MPI)

#-----------------------------------------------
# Create conversion utilities
#-----------------------------------------------
ADD_EXECUTABLE( 
  ConvertTo2H5Part 
  ConvertToH5Part.cxx
)

TARGET_LINK_LIBRARIES(
  ConvertTo2H5Part
  ${TOOL_VTK_LIBS}
)
