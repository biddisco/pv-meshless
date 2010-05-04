project(TipsyLib)

INCLUDE_DIRECTORIES (
  ${PROJECT_SOURCE_DIR} 
  ${PROJECT_BINARY_DIR}
)

#--------------------------------------------------
# Tipsy reader uses code from Doug Potter's Tipsy lib.
#--------------------------------------------------
ADD_LIBRARY(TipsyLib STATIC 
  adapter.cpp 
  binner.cpp 
	ftipsy.cpp 
	native.cpp 
	standard.cpp
	vtipsy.cpp
)

# Make sure static lib can be linked into dynamic one
IF(BUILD_SHARED_LIBS AND NOT WIN32)
  SET_TARGET_PROPERTIES( TipsyLib PROPERTIES COMPILE_FLAGS -fPIC)
ENDIF(BUILD_SHARED_LIBS AND NOT WIN32)
