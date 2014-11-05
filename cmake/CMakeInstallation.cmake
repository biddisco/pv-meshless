#-----------------------------------------------------------------------------
# Add Target(s) to CMake Install for import into other projects
#-----------------------------------------------------------------------------
if (NOT (${PROJECT_NAME}_EXTERNALLY_CONFIGURED OR ${PROJECT_NAME}_IS_SUBPROJECT))
  install (
      EXPORT ${${PROJECT_NAME}_EXPORTED_TARGETS}
      DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}/${${PROJECT_NAME}_PACKAGE}
      FILE ${${PROJECT_NAME}_PACKAGE}-targets.cmake
      COMPONENT configinstall
  )
endif()

#-----------------------------------------------------------------------------
# Export all exported targets to the build tree for use by parent project
#-----------------------------------------------------------------------------
if (NOT (${PROJECT_NAME}_EXTERNALLY_CONFIGURED OR ${PROJECT_NAME}_IS_SUBPROJECT))
  EXPORT (
      TARGETS ${${PROJECT_NAME}_LIBRARIES_TO_EXPORT} ${${PROJECT_NAME}_LIB_DEPENDENCIES}
      FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake
  )
endif()

#-----------------------------------------------------------------------------
# Configure the project-config.cmake file for the build directory
#-----------------------------------------------------------------------------
set (${PROJECT_NAME}_VERSION_STRING @${PROJECT_NAME}_PACKAGE_VERSION@)
set (${PROJECT_NAME}_VERSION_MAJOR  @${PROJECT_NAME}_PACKAGE_VERSION_MAJOR@)
set (${PROJECT_NAME}_VERSION_MINOR  @${PROJECT_NAME}_PACKAGE_VERSION_MINOR@)

configure_file (
    ${PROJECT_SOURCE_DIR}/cmake/config.cmake.build.in 
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake @ONLY
)

#-----------------------------------------------------------------------------
# Configure the FindXXX.cmake file for the install directory
#-----------------------------------------------------------------------------
#if (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)
#  configure_file (
#     ${PROJECT_SOURCE_DIR}/cmake/Find${PROJECT_NAME}.cmake.in 
#     ${PROJECT_BINARY_DIR}/CMakeFiles/Find${PROJECT_NAME}.cmake @ONLY
#  )
#  install (
#      FILES  ${PROJECT_BINARY_DIR}/CMakeFiles/Find${PROJECT_NAME}.cmake
#      DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}/${${PROJECT_NAME}_PACKAGE}
#      COMPONENT configinstall
#  )
#endif (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)

#-----------------------------------------------------------------------------
# Configure the project-config.cmake file for the install directory
#-----------------------------------------------------------------------------
if (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)
  configure_file (
     ${PROJECT_SOURCE_DIR}/cmake/config.cmake.install.in
       ${PROJECT_BINARY_DIR}/CMakeFiles/${${PROJECT_NAME}_PACKAGE}-config.cmake @ONLY
  )
  install (
      FILES  ${PROJECT_BINARY_DIR}/CMakeFiles/${${PROJECT_NAME}_PACKAGE}-config.cmake
      DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}/${${PROJECT_NAME}_PACKAGE}
      COMPONENT configinstall
  )
endif (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)

#-----------------------------------------------------------------------------
# Configure the project-config-version .cmake file for the install directory
#-----------------------------------------------------------------------------
if (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)
  configure_file (
     ${PROJECT_SOURCE_DIR}/cmake/config-version.cmake.in
       ${PROJECT_BINARY_DIR}/CMakeFiles/${${PROJECT_NAME}_PACKAGE}-config-version.cmake @ONLY
  )
  install (
      FILES  ${PROJECT_BINARY_DIR}/CMakeFiles/${${PROJECT_NAME}_PACKAGE}-config-version.cmake
      DESTINATION ${${PROJECT_NAME}_INSTALL_CMAKE_DIR}/${${PROJECT_NAME}_PACKAGE}
      COMPONENT configinstall
  )
endif (NOT ${PROJECT_NAME}_EXTERNALLY_CONFIGURED)
