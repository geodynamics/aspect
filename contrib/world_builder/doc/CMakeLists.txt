
set(WB_BUILD_DOCUMENTATION OFF CACHE BOOL "Whether or not to build the documentation. Requires doxygen, sphynx and breathe to be installed.")

if(WB_BUILD_DOCUMENTATION)
  MESSAGE(STATUS "Setting up documentation")
  find_package(Doxygen REQUIRED)

  # Find all the public headers
  get_target_property(GWB_PUBLIC_HEADER_DIR WorldBuilder INTERFACE_INCLUDE_DIRECTORIES)
  file(GLOB_RECURSE CAT_CUTIFIER_PUBLIC_HEADERS ${GWB_PUBLIC_HEADER_DIR}/*.h)

  set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR})
  set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
  set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
  set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxfile.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  
  #Replace variables inside @@ with the current values
  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
  
  file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
  add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                     DEPENDS ${CAT_CUTIFIER_PUBLIC_HEADERS}
                     COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                     MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     COMMENT "Generating docs")
  
 

  add_custom_target(doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

  find_package(Sphinx REQUIRED)


  set(SPHINXFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/conf.in.py)
  set(SPHINXFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/sphinx/conf.py)

  #Replace variables inside @@ with the current values
  configure_file(${SPHINXFILE_IN} ${SPHINXFILE_OUT_RTD} @ONLY)
  configure_file(${SPHINXFILE_IN} ${SPHINXFILE_OUT} @ONLY)

  set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
  set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

  add_custom_target(manual_html
                    COMMAND
                    ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/_static/ ${CMAKE_CURRENT_BINARY_DIR}/sphinx/_static/ &&
                    ${SPHINX_EXECUTABLE} -b html  -c ${CMAKE_CURRENT_BINARY_DIR}/sphinx -v
                    ${SPHINX_SOURCE} ${SPHINX_BUILD} 
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx
                    COMMENT "Generating documentation with Sphinx")


  add_custom_target(manual_pdf
                    COMMAND
                    ${CMAKE_COMMAND} -E copy_directory
                    ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/_static/ ${CMAKE_CURRENT_BINARY_DIR}/sphinx/_static/ &&
                    ${SPHINX_EXECUTABLE} -b latex  -c ${CMAKE_CURRENT_BINARY_DIR}/sphinx -v
                    ${SPHINX_SOURCE} ${SPHINX_BUILD} && cd ${SPHINX_BUILD} &&  make
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/sphinx
                    COMMENT "Generating documentation with Sphinx")



  add_custom_target(manual
                    COMMAND
                    ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target manual_html)

endif()
