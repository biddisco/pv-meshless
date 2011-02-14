
#--------------------------------------------------------
# Parallel compilation flag
#--------------------------------------------------------
IF (VTK_USE_MPI)
  ADD_DEFINITIONS("-DH5PART_HAS_MPI")
ENDIF (VTK_USE_MPI)

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

#--------------------------------------------------------
# Source files
#--------------------------------------------------------
SET(H5Part_SRCS
  ${CMAKE_CURRENT_SOURCE_DIR}/H5Part.c
  ${CMAKE_CURRENT_SOURCE_DIR}/h5tools_type.c
)

#
# Make sure msvc builds as cxx not C
# This resolves some struct typedef issues
#
IF(WIN32)
  SET_SOURCE_FILES_PROPERTIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/H5Part.c
    ${CMAKE_CURRENT_SOURCE_DIR}/H5Block.c
    PROPERTIES LANGUAGE CXX 
  )
ENDIF (WIN32)
  
#--------------------------------------------------------
# List the libraries used
#--------------------------------------------------------
SET(H5Part_LIBS  
  ${PARAVIEW_HDF5_LIBRARIES}
)

IF (VTK_USE_MPI)
  SET(H5Part_LIBS
    ${H5Part_LIBS} 
    ${MPI_LIBRARY}
  )
ENDIF (VTK_USE_MPI)

#--------------------------------------------------------
# Include class, configuration headers
#--------------------------------------------------------
SET(H5PART_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_BINARY_DIR}
  ${HDF5_INCLUDE_DIR}
)

INCLUDE_DIRECTORIES(${H5PART_INCLUDE_DIRS})

#--------------------------------------------------------
# Create the Library
#--------------------------------------------------------
ADD_LIBRARY(H5Part STATIC
  ${H5Part_SRCS}
)

TARGET_LINK_LIBRARIES(H5Part
  ${H5Part_LIBS}
)

# Make sure static lib can be linked into dynamic one
IF(BUILD_SHARED_LIBS AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_COMPILER_IS_GNUCC)
  SET_TARGET_PROPERTIES( H5Part PROPERTIES COMPILE_FLAGS -fPIC)  
ENDIF(BUILD_SHARED_LIBS AND CMAKE_SYSTEM_NAME STREQUAL "Linux" AND CMAKE_COMPILER_IS_GNUCC)

#--------------------------------------------------------
# Add H5Part conversion tools
#--------------------------------------------------------
ADD_SUBDIRECTORY(tools)