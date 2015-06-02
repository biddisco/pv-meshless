PROJECT(H5PartTools)

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
#  pv_common
  pv_meshless
  vtkIOInfovis
  vtkTestingRendering
#  vtkRendering
#  vtkGraphics
#  vtkImaging
#  vtkIO
#  vtkFiltering
#  vtkCommon
  vtksys
)

IF(PARAVIEW_USE_MPI)
  SET(TOOL_VTK_LIBS ${TOOL_VTK_LIBS} ${MPI_LIBRARY} )
  IF(MPI_EXTRA_LIBRARY)
    SET(TOOL_VTK_LIBS ${TOOL_VTK_LIBS} ${MPI_EXTRA_LIBRARY} )
  ENDIF(MPI_EXTRA_LIBRARY)
ENDIF(PARAVIEW_USE_MPI)

#-----------------------------------------------
# Create vtk conversion utilities
#-----------------------------------------------
ADD_EXECUTABLE( 
  VTM2H5Part 
  VTM2H5Part.cpp
)

TARGET_LINK_LIBRARIES(
  VTM2H5Part
  ${TOOL_VTK_LIBS}
)

#-----------------------------------------------
# Create Ensight conversion utility
#-----------------------------------------------
ADD_EXECUTABLE( 
  Ensight2H5Part 
  Ensight2H5Part.cpp 
)

TARGET_LINK_LIBRARIES(
  Ensight2H5Part 
  ${TOOL_VTK_LIBS}
)

#-----------------------------------------------
# Create CSV conversion utility
#-----------------------------------------------
ADD_EXECUTABLE( 
  CSV2H5Part 
  CSV2H5Part.cpp 
)

TARGET_LINK_LIBRARIES(
  CSV2H5Part
  ${TOOL_VTK_LIBS}
)

#-----------------------------------------------
# Create ASCII conversion utility
#-----------------------------------------------
if(DO_THIS)
ADD_EXECUTABLE( 
  Ascii2H5Part 
  vtkASCIIParticleReader.cpp
  Ascii2H5Part.cpp 
)

TARGET_LINK_LIBRARIES(
  Ascii2H5Part
  ${TOOL_VTK_LIBS}
)

SET(TOOL_CFG_FILES
  ASCII2H5Part.ecn.cfg
  ASCII2H5Part.edf.cfg
  ASCII2H5Part.manchester.cfg
  ASCII2H5Part.plymouth.cfg
)
endif()

INSTALL(
  FILES ${TOOL_CFG_FILES}
  DESTINATION bin/
)
