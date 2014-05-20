PROJECT("pv_meshless")

#--------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------
cmake_minimum_required(VERSION 2.8)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
  cmake_policy(SET CMP0023 OLD)
endif(COMMAND cmake_policy)

#-----------------------------------------------------------------------------
# Basic naming/config stuff here
#-----------------------------------------------------------------------------
set (pv_meshless_VERS_MAJOR   0)
set (pv_meshless_VERS_MINOR   9)
set (pv_meshless_VERS_RELEASE 0)

set (pv_meshless_PACKAGE                "pv_meshless")
set (pv_meshless_PACKAGE_NAME           "pv_meshless")
set (pv_meshless_PACKAGE_VERSION        "${pv_meshless_VERS_MAJOR}.${pv_meshless_VERS_MINOR}.${pv_meshless_VERS_RELEASE}")
set (pv_meshless_PACKAGE_VERSION_MAJOR  "${pv_meshless_VERS_MAJOR}.${pv_meshless_VERS_MINOR}")
set (pv_meshless_PACKAGE_VERSION_MINOR  "${pv_meshless_VERS_RELEASE}")
set (pv_meshless_PACKAGE_VERSION_STRING "${pv_meshless_PACKAGE_VERSION}")
set (pv_meshless_PACKAGE_SOVERSION      "${pv_meshless_SOVERS_MAJOR}.${pv_meshless_SOVERS_RELEASE}.${pv_meshless_SOVERS_MINOR}")
set (pv_meshless_PACKAGE_STRING         "${pv_meshless_PACKAGE_NAME} ${pv_meshless_PACKAGE_VERSION_STRING}")
set (pv_meshless_PACKAGE_TARNAME        "${pv_meshless_PACKAGE}")
set (pv_meshless_PACKAGE_URL            "http://www.cscs.ch")
set (pv_meshless_PACKAGE_BUGREPORT      "help@cscs.ch")

#--------------------------------------------------
# Location of cmake files
#--------------------------------------------------
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  IF(PARAVIEW_BUILD_QT_GUI)
    INCLUDE(${QT_USE_FILE})
  ENDIF(PARAVIEW_BUILD_QT_GUI)
  INCLUDE_DIRECTORIES(
    ${HDF5_INCLUDE_DIR}
    ${PARAVIEW_INCLUDE_DIRS}
    ${PARAVIEW_GUI_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

#--------------------------------------------------
# pv-zoltan : Required for parallel sampling
#--------------------------------------------------
find_package(pv_zoltan NO_MODULE)
if(NOT pv_zoltan_FOUND)
  message(FATAL_ERROR "pv-meshless requires pv_zoltan")
endif()
include_directories(${pv_zoltan_INCLUDE_DIR})

#--------------------------------------------------
# Set project include directories 
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
)

#--------------------------------------------------
# Set Includes
#--------------------------------------------------
  INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/sph
  ${CMAKE_CURRENT_SOURCE_DIR}/h5part
  ${CMAKE_CURRENT_SOURCE_DIR}/cosmo
  ${CMAKE_CURRENT_SOURCE_DIR}/tipsylib
)

#------------------------------------------------------
# MPI support if ParaView compiled with MPI
#------------------------------------------------------
IF (PARAVIEW_USE_MPI)
  FIND_PACKAGE(MPI)   
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ENDIF (PARAVIEW_USE_MPI)

#-----------------------------------------------------------------------------
# OpenMP
#-----------------------------------------------------------------------------
OPTION(PV_MESHLESS_USE_OPENMP "Compile pv-meshless with OpenMP support" ON)
IF (PV_MESHLESS_USE_OPENMP)
  ADD_DEFINITIONS(-DHAVE_OPENMP)
  SET(OpenMP_CXX_FLAGS "")
  IF (WIN32 AND MSVC) 
    SET(OpenMP_CXX_FLAGS "/openmp")
  ENDIF (WIN32 AND MSVC) 
ENDIF (PV_MESHLESS_USE_OPENMP)

#--------------------------------------------------
# pv-zoltan
#--------------------------------------------------
find_package(pv_zoltan)
if (pv_zoltan_FOUND)
  set(ZOLTAN_SUPPORT 1)
  include_directories(
    ${pv_zoltan_INCLUDE_DIR}
  )
  if(WIN32)
    include_directories(
#      ${ZOLTAN_DIR}/trilinos/cmake/tribits/common_tools/win_interface/include
   )
  endif(WIN32)
  ADD_DEFINITIONS(-DPV_MESHLESS_ZOLTAN_SUPPORT)
  set(ZOLTAN_LIBS pv_zoltan)
  message("Zoltan support enabled, libs are ${ZOLTAN_LIBS}")
endif()

#--------------------------------------------------
# Set Definitions 
#--------------------------------------------------
ADD_DEFINITIONS(-DMPICH_SKIP_MPICXX)
IF (WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

# --------------------------------------------------
# Option to disable hdf5 readers etc
# --------------------------------------------------
OPTION(PV_MESHLESS_USE_HDF5 "Use HDF5 readers/writers" ON)

# --------------------------------------------------
# Configure all HDF5 options
# --------------------------------------------------
IF(PV_MESHLESS_USE_HDF5) 
  #-----------------------------------------------
  # If (paraview) using system hdf5, find it
  #-----------------------------------------------
#  IF (NOT TARGET "vtkhdf5")
    FIND_PACKAGE(HDF5 REQUIRED NO_MODULE PATHS "$ENV{ProgramFiles}/HDF5")
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
#  ENDIF (NOT TARGET "vtkhdf5")
  
  #-----------------------------------------------
  # Check if HDF5 was compiled with Parallel support
  #-----------------------------------------------
  INCLUDE (${CMAKE_ROOT}/Modules/CheckSymbolExists.cmake)
  UNSET(HAVE_H5_HAVE_PARALLEL CACHE)
  SET(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIR} ${VTK_INCLUDE_DIRS})
  CHECK_SYMBOL_EXISTS(H5_HAVE_PARALLEL "H5pubconf.h" HAVE_H5_HAVE_PARALLEL)

  #-----------------------------------------------
  # Enable parallel IO if available
  #-----------------------------------------------
  IF(HAVE_H5_HAVE_PARALLEL)
    ADD_DEFINITIONS(-DPARALLEL_IO)
  ELSE(HAVE_H5_HAVE_PARALLEL)
    MESSAGE("H5_HAVE_PARALLEL not found, Parallel HDF5 IO will not be enabled")
  ENDIF(HAVE_H5_HAVE_PARALLEL)
  
  #-----------------------------------------------
  # Build minimal H5Part lib for support
  #-----------------------------------------------
  ADD_SUBDIRECTORY(h5part)
  SET(H5PART_LIBS H5Part)
  
  SET(HDF5_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.xml
  )
  SET(HDF5_CXX
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartWriter.cxx
  )
ENDIF(PV_MESHLESS_USE_HDF5) 
  
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
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridPanel.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.ui
#    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.ui
  )

  #--------------------------------------------------
  # make sure all *.h files for ui are wrapped by moc
  #--------------------------------------------------
  QT4_WRAP_CPP(
    MESHLESS_MOC_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridPanel.h      
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridFilterPanel.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.h
  )
  
  #--------------------------------------------------
  # invoke macro that adds our grid source 3D widget
  #--------------------------------------------------
#  ADD_3DWIDGET(
#    MESHLESS_IFACE_W
#    MESHLESS_IFACE_SRCS_W
#    CLASS_NAME 
#      pqRegularGridSourceWidget
#    WIDGET_TYPE 
#      # string found in <hints> section of xml
#      "RegularGridSource"
#  )

  #--------------------------------------------------
  # invoke macro that adds our SPH manager panel
  #--------------------------------------------------
  ADD_PARAVIEW_OBJECT_PANEL(
    MESHLESS_IFACE_O
    MESHLESS_IFACE_SRCS_O
    CLASS_NAME 
      pqSPHManagerPanel
    XML_NAME 
     "SPHManager"
    XML_GROUP 
     meshless_helpers
  )

  #--------------------------------------------------
  # invoke macro that adds our Sampling grid panel
  #--------------------------------------------------
  ADD_PARAVIEW_OBJECT_PANEL(
    MESHLESS_IFACE_O
    MESHLESS_IFACE_SRCS_O
    CLASS_NAME 
      pqSamplingGridPanel
    XML_NAME 
     "SamplingGridSource"
    XML_GROUP 
     sources
  )

  #--------------------------------------------------
  # invoke macro that adds our Sampling grid panel
  #--------------------------------------------------
  ADD_PARAVIEW_OBJECT_PANEL(
    MESHLESS_IFACE_O2
    MESHLESS_IFACE_SRCS_O2
    CLASS_NAME 
      pqSamplingGridFilterPanel
    XML_NAME 
     "SamplingGridFilter"
    XML_GROUP 
     filters
  )

  #--------------------------------------------------
  # invoke macro that adds our SPH manager dock window
  #--------------------------------------------------
  ADD_PARAVIEW_DOCK_WINDOW(
    MESHLESS_IFACE_D
    MESHLESS_IFACE_SRCS_D
    CLASS_NAME 
      pqSPHManagerDockWindow 
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
    ${CMAKE_CURRENT_SOURCE_DIR}/sph/KernelWendland.cpp
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
SET(PLUGIN_NAME pv_meshless)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 
  REQUIRED_PLUGINS 
# pv_common

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSamplingGridGenerator.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.xml
#    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter3.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleIdFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractValueFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractDataOverTime.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkImageSamplerSource.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHImageResampler.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkClipPolyData.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalPlotValueFilter.xml
    ${HDF5_XML}
  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSamplingGridGenerator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSMCustomBoxRepresentationProxy.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.cxx
#    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter3.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTree.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleIdFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractValueFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkImageSamplerSource.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHImageResampler.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalPlotValueFilter.cxx
    ${HDF5_CXX}
  SERVER_SOURCES
    ${KernelFiles_SRCS}
  GUI_INTERFACES 
    ${MESHLESS_IFACE_W}
    ${MESHLESS_IFACE_O}
    ${MESHLESS_IFACE_O2}
    ${MESHLESS_IFACE_D}
  GUI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridPanel.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridFilterPanel.cxx
    ${MESHLESS_IFACE_SRCS_W}
    ${MESHLESS_IFACE_SRCS_O}
    ${MESHLESS_IFACE_SRCS_O2}
    ${MESHLESS_IFACE_SRCS_D}
    ${MESHLESS_MOC_SRCS}
    ${MESHLESS_UI_SOURCES}
  GUI_RESOURCES
    # none for now  
  GUI_RESOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Sources.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Filters.xml
#    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Readers.xml
)

SET_TARGET_PROPERTIES(${PLUGIN_NAME} PROPERTIES COMPILE_FLAGS "${OpenMP_CXX_FLAGS}")

TARGET_LINK_LIBRARIES(${PLUGIN_NAME}
  PUBLIC
    vtksys
    vtkFiltersFlowPaths
    ${H5PART_LIBS}
#    pv_common
    ${ZOLTAN_LIBS}
    vtkInteractionWidgets
)

# --------------------------------------------------
# Tools 
# --------------------------------------------------
IF(EXISTS "${PROJECT_SOURCE_DIR}/../pv-astro/")
  OPTION(USE_ASTRO "Use pv-astro library for ramses reader/converter" OFF)
  IF(USE_ASTRO)
    ADD_DEFINITIONS("-DUSE_PV_ASTRO")
    INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/../pv-astro/")
    ADD_SUBDIRECTORY(tools)
  ENDIF(USE_ASTRO)
ENDIF(EXISTS "${PROJECT_SOURCE_DIR}/../pv-astro/")

#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(Testing)
ENDIF (BUILD_TESTING)  

#--------------------------------------------------------
# Create the configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/module_configure.h.in
               ${PROJECT_BINARY_DIR}/pv_meshless_configure.h
               @ONLY IMMEDIATE)

#-----------------------------------------------------------------------------
# Configure the project-config.cmake file for the build directory
#-----------------------------------------------------------------------------
set (${PROJECT_NAME}_INCLUDES_BUILD_TIME
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
)
include(CMakeInstallation)

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
    ${PLUGIN_NAME}
  DESTINATION 
    ${INSTALL_PATH}
)
