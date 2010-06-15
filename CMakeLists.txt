PROJECT("pv-meshless")

CMAKE_MINIMUM_REQUIRED (VERSION 2.8)

IF(COMMAND CMAKE_POLICY)
  CMAKE_POLICY(SET CMP0003 NEW)
ENDIF(COMMAND CMAKE_POLICY)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF(NOT PARAVIEW_LOADED)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
  SET(PARAVIEW_LOADED 1)
ELSE(NOT PARAVIEW_LOADED)
  IF (PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF (PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${PARAVIEW_INCLUDE_DIRS}
    ${ParaView_SOURCE_DIR}/VTK/Geovis
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${PARAVIEW_KWSYS_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    ${MPI_INCLUDE_PATH}
  )
ENDIF(NOT PARAVIEW_LOADED)

#--------------------------------------------------
# Find and Use vtkCSCSCommon
#--------------------------------------------------
FIND_PATH(CSCS_COMMON_DIR 
  UsevtkCSCSCommon.cmake 
  PATHS 
    "${CSCS_COMMON_DIR}"
    "${PROJECT_BINARY_DIR}/../Common"  
)
IF (EXISTS ${CSCS_COMMON_DIR}/UsevtkCSCSCommon.cmake)
  INCLUDE(${CSCS_COMMON_DIR}/UsevtkCSCSCommon.cmake)
ELSE (EXISTS ${CSCS_COMMON_DIR}/UsevtkCSCSCommon.cmake)
  MESSAGE(FATAL_ERROR "We need vtkCSCSCommon for Project ${PROJECT_NAME}")
ENDIF(EXISTS ${CSCS_COMMON_DIR}/UsevtkCSCSCommon.cmake)

#--------------------------------------------------
# Set Includes
#--------------------------------------------------
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sph)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/h5part)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/cosmo)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/tipsylib)

#------------------------------------------------------
# Parallel IO (H5Part Readers and Writers)
#------------------------------------------------------
IF (VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_DEFINITIONS(-DPARALLEL_IO)
ELSE (VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_DEFINITIONS(-DUSE_SERIAL_COSMO)
ENDIF (VTK_USE_PARALLEL AND VTK_USE_MPI)

#--------------------------------------------------
# Tipsy reader uses code from Doug Potter's Tipsy lib.
#--------------------------------------------------
ADD_SUBDIRECTORY(tipsylib)

#--------------------------------------------------
# Build minimal H5Part lib for support
#--------------------------------------------------
ADD_SUBDIRECTORY(h5part)

#--------------------------------------------------
# Build minimal Cosmo lib for support
#--------------------------------------------------
ADD_SUBDIRECTORY(cosmo)

#--------------------------------------------------
# Set Definitions 
#--------------------------------------------------
ADD_DEFINITIONS(-DMPICH_SKIP_MPICXX)
IF (WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

#--------------------------------------------------
# Setup GUI custom Qt panel sources and wrapping
#--------------------------------------------------
IF(PARAVIEW_BUILD_QT_GUI)

  #--------------------------------------------------
  # make sure all *.ui dialogs/etc are wrapped
  #--------------------------------------------------
  QT4_WRAP_UI(
    MESHLESS_UI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.ui
  )

  #--------------------------------------------------
  # make sure all *.h files for ui are wrapped by moc
  #--------------------------------------------------
  QT4_WRAP_CPP(
    MESHLESS_MOC_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.h
  )
  
  #--------------------------------------------------
  # invoke macro that adds our grid source 3D widget
  #--------------------------------------------------
  ADD_3DWIDGET(
    MESHLESS_IFACE_W
    MESHLESS_IFACE_SRCS_W
    CLASS_NAME 
      pqRegularGridSourceWidget
    WIDGET_TYPE 
      # string found in <hints> section of xml
      "RegularGridSource"
  )

  #--------------------------------------------------
  # invoke macro that adds our SPH manager dock window
  #--------------------------------------------------
  ADD_PARAVIEW_DOCK_WINDOW(
    MESHLESS_IFACE_D
    MESHLESS_IFACE_SRCS_D
    CLASS_NAME 
      pqSPHManagerPanel 
    DOCK_AREA 
      Left Right Top Bottom
  )

ENDIF(PARAVIEW_BUILD_QT_GUI)

#--------------------------------------------------
# Define Wrapping hints
#--------------------------------------------------
SET(VTK_WRAP_HINTS "${CMAKE_CURRENT_SOURCE_DIR}/hints.txt" )
 
#--------------------------------------------------
# Tell Visual Studio to group these together 
#--------------------------------------------------
SET(KernelFiles_SRCS 
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/Vector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/Kernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelGaussian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelSpline3rdOrder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelSpline5thOrder.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelQuadratic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelBox.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelCusp.cpp
)

SOURCE_GROUP(KernelFiles FILES 
  ${KernelFiles_SRCS}
)

#--------------------------------------------------
# Define Plugin 
#--------------------------------------------------
SET(LIBRARY_TARGET_NAME pvmeshless_plugin)
ADD_PARAVIEW_PLUGIN(
  ${LIBRARY_TARGET_NAME}
  "1.0" 
  REQUIRED_PLUGINS CSCSCommon_plugin

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTipsyReader.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPTipsyReader.xml
  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSMCustomBoxRepresentationProxy.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTree.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTipsyReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartWriter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5SPHReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkPTipsyReader.cxx
  SERVER_SOURCES
    ${KernelFiles_SRCS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AstroVizHelpers.cxx
  GUI_INTERFACES 
    ${MESHLESS_IFACE_W}
    ${MESHLESS_IFACE_D}
  GUI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.cxx
    ${MESHLESS_IFACE_SRCS_W}
    ${MESHLESS_IFACE_SRCS_D}
    ${MESHLESS_MOC_SRCS}
    ${MESHLESS_UI_SOURCES}
  GUI_RESOURCES
    # none for now  
  GUI_RESOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Sources.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Filters.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Readers.xml
)

TARGET_LINK_LIBRARIES(pvmeshless_plugin 
  vtkWidgets
  H5Part 
  TipsyLib 
  Cosmo2
  CSCSCommon_plugin
)

#--------------------------------------------------------
# Create the configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/vtkLocalConfigure.h.in
               ${PROJECT_BINARY_DIR}/vtkCSCSMeshlessConfigure.h
               @ONLY IMMEDIATE)

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
#  ADD_SUBDIRECTORY(Testing)
ENDIF (BUILD_TESTING)  

#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH 
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
)

#INSTALL(
#  FILES ${TOOL_CFG_FILES}
#  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
#)

INSTALL(
  TARGETS
    pvmeshless_plugin
  DESTINATION 
    ${INSTALL_PATH}
)
