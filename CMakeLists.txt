PROJECT("pv-meshless")

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

#--------------------------------------------------
# Set project include directories 
# pv-common is assumed to be checkout out too
#--------------------------------------------------
INCLUDE_DIRECTORIES(
  "${PROJECT_SOURCE_DIR}"
  "${PROJECT_BINARY_DIR}"
  "${PROJECT_SOURCE_DIR}/../pv-common"
  "${PROJECT_BINARY_DIR}/../pv-common"
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

# --------------------------------------------------
# Option to disable hdf5 readers etc
# --------------------------------------------------
OPTION(USE_HDF5 "Use HDF5 readers/writers" ON)

# --------------------------------------------------
# Check if HDF5 was compiled with Parallel support
# --------------------------------------------------
IF(USE_HDF5) 
  INCLUDE (${CMAKE_ROOT}/Modules/CheckSymbolExists.cmake)
  CHECK_SYMBOL_EXISTS(H5_HAVE_PARALLEL H5pubconf.h HAVE_H5_HAVE_PARALLEL)
  IF(NOT HAVE_H5_HAVE_PARALLEL)
    MESSAGE("H5_HAVE_PARALLEL was not found, Parallel IO will not be enabled")
  ENDIF(NOT HAVE_H5_HAVE_PARALLEL)
ENDIF(USE_HDF5) 
  
#------------------------------------------------------
# Parallel IO (H5Part Readers and Writers)
#------------------------------------------------------
IF (VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_DEFINITIONS(-DPARALLEL_IO)
  FIND_PACKAGE(MPI REQUIRED)
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
ELSE (VTK_USE_PARALLEL AND VTK_USE_MPI)
  ADD_DEFINITIONS(-DUSE_SERIAL_COSMO)
ENDIF (VTK_USE_PARALLEL AND VTK_USE_MPI)

# --------------------------------------------------
# Optional Trilinos support 
# --------------------------------------------------
OPTION(USE_TRILINOS "Use Trilinos library" OFF)
IF(USE_TRILINOS)
  SET(Trilinos_ENABLE_Zoltan 1)
  SET(Trilinos_ENABLE_Tests 0)
  SET(TPL_ENABLE_MPI ON)
  ADD_SUBDIRECTORY(trilinos)
  SET(Trilinos_FOUND 1)
#  FIND_PACKAGE(Trilinos)
  IF(${Trilinos_FOUND})
    ADD_DEFINITIONS(-DPVMESHLESS_USE_TRILINOS)
    INCLUDE_DIRECTORIES("${Trilinos_INCLUDE_DIRS}")
    LINK_DIRECTORIES(${Trilinos_LIBRARY_DIRS})
    SET(TRILINOS_CXX ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.cxx)
    SET(TRILINOS_XML ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionFilter.xml)
#    MESSAGE("Trilinos_LIBRARY_DIRS ${Trilinos_LIBRARY_DIRS}/${CMAKE_CFG_INTDIR}") 
  ENDIF(${Trilinos_FOUND})
ELSE(USE_TRILINOS)
  SET(Trilinos_FOUND 0)
ENDIF(USE_TRILINOS)

#--------------------------------------------------
# Build minimal H5Part lib for support
#-----------------------------------------------
IF(USE_HDF5)
  ADD_SUBDIRECTORY(h5part)
  SET(H5PART_LIBS H5Part)
ENDIF(USE_HDF5)

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
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridPanel.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.ui
#    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.ui
  )

  #--------------------------------------------------
  # make sure all *.h files for ui are wrapped by moc
  #--------------------------------------------------
  QT4_WRAP_CPP(
    MESHLESS_MOC_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridPanel.h      
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSamplingGridFilterPanel.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.h
      ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.h
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

IF(USE_HDF5)
  SET(HDF5_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.xml
  )
  SET(HDF5_CXX
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5PartWriter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkH5SPHReader.cxx
  )
ENDIF(USE_HDF5)
  
#--------------------------------------------------
# Define Plugin 
#--------------------------------------------------
SET(PLUGIN_NAME pv_meshless)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 
  REQUIRED_PLUGINS pv_common

  SERVER_MANAGER_XML 
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSamplingGridGenerator.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter3.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleIdFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractValueFilter.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractDataOverTime.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMultiBlockExtractDataOverTime.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkProbeFileReader.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkImageSamplerSource.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHImageResampler.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionRepresentation.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkClipPolyData.xml
    ${TRILINOS_XML}
    ${HDF5_XML}
    ${TIPSY_XML}
  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkBoundsExtentTranslator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkRegularGridSource.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxWidget.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCustomBoxRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSMCustomBoxRepresentationProxy.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHManager.cxx
#    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSMSPHManagerProxy.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHProbeFilter3.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalParticleDataInterpolator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTree.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleBoxTreeRepresentation.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticleIdFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkCenterOfMassFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMomentsOfInertiaFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSamplingGridGenerator.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExtractValueFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkProbeFileReader.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkMultiBlockExtractDataOverTime.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkImageSamplerSource.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkSPHImageResampler.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkParticlePartitionRepresentation.cxx
    ${TRILINOS_CXX}
    ${HDF5_CXX}
    ${TIPSY_CXX}
  SERVER_SOURCES
    ${KernelFiles_SRCS}
#    ${CMAKE_CURRENT_SOURCE_DIR}/AstroVizHelpers.cxx  
  GUI_INTERFACES 
    ${MESHLESS_IFACE_W}
    ${MESHLESS_IFACE_O}
    ${MESHLESS_IFACE_O2}
    ${MESHLESS_IFACE_D}
  GUI_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerPanel.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqSPHManagerDockWindow.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.cxx
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
    ${CMAKE_CURRENT_SOURCE_DIR}/pqRegularGridSourceWidget.ui
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Sources.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Filters.xml
    ${CMAKE_CURRENT_SOURCE_DIR}/Meshless_Readers.xml
)

TARGET_LINK_LIBRARIES(${PLUGIN_NAME}
  ${H5PART_LIBS}
  pv_common
  vtkWidgets
  ${TIPSY_LIB}
  ${Trilinos_LIBRARIES}
)

#--------------------------------------------------------
# Create the configuration header.
#--------------------------------------------------------
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/pv_meshless_configure.h.in
               ${PROJECT_BINARY_DIR}/vtkCSCSMeshlessConfigure.h
               @ONLY IMMEDIATE)

# --------------------------------------------------
# Tools 
# --------------------------------------------------
OPTION(USE_ASTRO "Use pv-astro library for ramses reader/converter" OFF)
IF(USE_ASTRO)
  ADD_DEFINITIONS("-DUSE_PV_ASTRO")
  INCLUDE_DIRECTORIES(${pv-astro_SOURCE_DIR})
  ADD_SUBDIRECTORY(tools)
ENDIF(USE_ASTRO)
  
#--------------------------------------------------
# Testing
#--------------------------------------------------
IF (BUILD_TESTING)
  SET(PLUGIN_TEST_DIR ${PROJECT_BINARY_DIR}/Testing/Temporary)
  MAKE_DIRECTORY(${PLUGIN_TEST_DIR})
  ADD_SUBDIRECTORY(Testing)
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
    ${PLUGIN_NAME}
  DESTINATION 
    ${INSTALL_PATH}
)
