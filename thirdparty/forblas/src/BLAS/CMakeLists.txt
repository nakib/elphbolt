enable_language(Fortran)

# Check for any necessary platform specific compiler flags
include(CheckLAPACKCompilerFlags)
CheckLAPACKCompilerFlags()

add_subdirectory(SRC)
if(BUILD_TESTING)
  add_subdirectory(TESTING)
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/blas.pc.in ${CMAKE_CURRENT_BINARY_DIR}/${BLASLIB}.pc @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${BLASLIB}.pc
  DESTINATION ${PKG_CONFIG_DIR}
  COMPONENT Development
  )
