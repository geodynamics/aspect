# find diff program
FIND_PROGRAM(DIFF_EXECUTABLE
	     NAMES numdiff diff FC
	     HINTS ${DIFF_DIR}
	     PATH_SUFFIXES bin
	     )

 IF(NOT DIFF_EXECUTABLE MATCHES "-NOTFOUND")
	 SET(TEST_DIFF ${DIFF_EXECUTABLE})
	 MESSAGE(STATUS "Found ${DIFF_EXECUTABLE} to test the difference between files.")
 ELSE()
	     MESSAGE(FATAL_ERROR
		     "Could not find numdiff, diff or fc. This is required for running the testsuite.\n"
		     "Please specify TEST_DIFF by hand."
		     )
ENDIF()

# Unit test testing
#find all the unit test files
file(GLOB_RECURSE UNIT_TEST_SOURCES_CXX "unit_tests/*.cc")

if(CMAKE_Fortran_COMPILER)
  file(GLOB_RECURSE UNIT_TEST_SOURCES_FORTAN "unit_tests/*.f90")
endif()

set(UNIT_TEST_SOURCES ${UNIT_TEST_SOURCES_CXX} ${UNIT_TEST_SOURCES_FORTAN})


#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/)

# Run through each source
foreach(test_source ${UNIT_TEST_SOURCES})
  # Extract the filename without an extension (NAME_WE)
  get_filename_component(test_name ${test_source} NAME_WE)

  # Add compile target
  add_executable(${test_name} ${test_source})

  if(NOT ${CMAKE_VERSION} VERSION_LESS "3.12.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
    target_compile_options(${test_name} INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE})
  endif()

  if(${CMAKE_VERSION} VERSION_LESS "3.13.0") 
    SET(CMAKE_EXE_LINKER_FLAGS "${WB_LINKER_OPTIONS}")
  else()
    target_link_options(${test_name} INTERFACE ${WB_LINKER_OPTIONS})
  endif()

  # Make sure that the whole library is loaded, so the registration is done correctly.
  if(NOT APPLE)
    SET(GWB_LIBRARY_WHOLE -Wl,--whole-archive WorldBuilder -Wl,--no-whole-archive)  
  elseif(MSVC)
    SET(GWB_LIBRARY_WHOLE WorldBuilder)
    SET(CMAKE_EXE_LINKER_FLAGS "/WHOLEARCHIVE:libworldbuilder")
  else()
    SET(GWB_LIBRARY_WHOLE -Wl,-force_load WorldBuilder)
  endif()
	
  target_link_libraries(${test_name} ${GWB_LIBRARY_WHOLE})

  #Move testing binaries into a testBin directory
  set_target_properties(${test_name} PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/bin/)

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${test_name}
           WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
           COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name})
endforeach(test_source)

# App testing
IF(CMAKE_RUN_APP_TESTS)
  # Create directory for the test results
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/app)

  # Test help
  add_test(testing_help
           ${CMAKE_COMMAND} 
  	 -D TEST_NAME=testing_help 
     -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/WorldBuilderApp${CMAKE_EXECUTABLE_SUFFIX} 
  	 -D TEST_ARGS=--help 
  	 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/app/testing_help/screen-output.log 
     -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/app/testing_help/screen-output.log
     -D TEST_DIFF=${TEST_DIFF}
  	 -P ${CMAKE_SOURCE_DIR}/tests/app/run_app_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/app/)
  
  # Test no files provided
  add_test(testing_no_file
           ${CMAKE_COMMAND} 
  	 -D TEST_NAME=testing_no_file 
  	 -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/WorldBuilderApp${CMAKE_EXECUTABLE_SUFFIX} 
  	 -D TEST_ARGS= 
  	 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/app/testing_no_file/screen-output.log 
  	 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/app/testing_no_file/screen-output.log
     -D TEST_DIFF=${TEST_DIFF}
  	 -P ${CMAKE_SOURCE_DIR}/tests/app/run_app_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/app/)
  
  
  # Test one file provided 
  add_test(testing_one_file
           ${CMAKE_COMMAND} 
  	 -D TEST_NAME=testing_one_file 
  	 -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/WorldBuilderApp${CMAKE_EXECUTABLE_SUFFIX} 
  	 -D TEST_ARGS=non_existend_file
  	 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/app/testing_one_file/screen-output.log 
  	 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/app/testing_one_file/screen-output.log
     -D TEST_DIFF=${TEST_DIFF}
  	 -P ${CMAKE_SOURCE_DIR}/tests/app/run_app_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/app/)
  
  #find all the integration test files
  file(GLOB_RECURSE APP_TEST_SOURCES "app/*.wb")
  
  # Run through each sourceUforeach(test_source ${APP_TEST_SOURCES})
  foreach(test_source ${APP_TEST_SOURCES})
          get_filename_component(test_name ${test_source} NAME_WE)
  	set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/app/${test_name}.wb\;${CMAKE_SOURCE_DIR}/tests/app/${test_name}.dat")
          add_test(${test_name}
                   ${CMAKE_COMMAND} 
  	         -D TEST_NAME=${test_name}
  	         -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/WorldBuilderApp${CMAKE_EXECUTABLE_SUFFIX} 
  	         -D TEST_ARGS=${TEST_ARGUMENTS}
  	         -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/app/${test_name}/screen-output.log 
  	         -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/app/${test_name}/screen-output.log
             -D TEST_DIFF=${TEST_DIFF}
  	         -P ${CMAKE_SOURCE_DIR}/tests/app/run_app_tests.cmake
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/app/)
  endforeach(test_source)
ENDIF()

# visualization tests
 # Create directory for the test results
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/visualization)
#find all the integration visu files
file(GLOB_RECURSE VISU_TEST_SOURCES "visualization/*.wb")
# Run through each sourceUforeach(test_source ${VISU_TEST_SOURCES})
foreach(test_source ${VISU_TEST_SOURCES})
        get_filename_component(test_name ${test_source} NAME_WE)
	set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/visualization/${test_name}.wb\;${CMAKE_SOURCE_DIR}/tests/visualization/${test_name}.grid")
        add_test(${test_name}
                 ${CMAKE_COMMAND} 
	         -D TEST_NAME=${test_name}
	         -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/WorldBuilderVisualization${CMAKE_EXECUTABLE_SUFFIX} 
	         -D TEST_ARGS=${TEST_ARGUMENTS}
	         -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/visualization/${test_name}.vtu 
	         -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/visualization/${test_name}.vtu
	         -P ${CMAKE_SOURCE_DIR}/tests/visualization/run_visualization_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/visualization/)
endforeach(test_source)

#test fortran compilation and wrapper if compiler found
if(CMAKE_Fortran_COMPILER)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran)
  if(NOT APPLE)
    add_test(NAME compile_simple_fortran_test
             COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test.f90 -L../../lib/ -Wl,--whole-archive -lWorldBuilder -Wl,--no-whole-archive -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lstdc++ 
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  else()
    add_test(NAME compile_simple_fortran_test
	    COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test.f90 -L../../lib/ -Wl,-force_load,../../lib/libWorldBuilder.a -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lc++ 
	    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  endif()
  add_test(run_simple_fortran_test
           ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_fortran_test
                 -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/fortran/test${CMAKE_EXECUTABLE_SUFFIX}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/fortran/run_simple_fortran_test.log
		 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/fortran/run_simple_fortran_test.log
                 -P ${CMAKE_SOURCE_DIR}/tests/fortran/run_fortran_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/) 
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran)

  if(NOT APPLE)
    add_test(NAME compile_simple_fortran_example 
	     COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/example.f90 -L../../lib/ -Wl,--whole-archive -lWorldBuilder -Wl,--no-whole-archive -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lstdc++ 
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  else()
    add_test(NAME compile_simple_fortran_example 
	     COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/example.f90 -L../../lib/ -Wl,-force_load,../../lib/libWorldBuilder.a -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lc++ 
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  endif()
  add_test(run_simple_fortran_example
	   ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_fortran_example
                 -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/fortran/example${CMAKE_EXECUTABLE_SUFFIX}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/fortran/run_simple_fortran_example.log
		 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/fortran/run_simple_fortran_example.log
                 -P ${CMAKE_SOURCE_DIR}/tests/fortran/run_fortran_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/) 
 endif()

#test python compilation and wrapper if compiler found
if(MAKE_PYTHON_WRAPPER)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/python)

   add_test(run_simple_python_test
           ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_python_test
                 -D TEST_PROGRAM=${Python_EXECUTABLE}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/python/test.py
		 -D TEST_INPUT=
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/python/run_simple_python_test.log
		 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/python/run_simple_python_test.log
                 -P ${CMAKE_SOURCE_DIR}/tests/python/run_python_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/python/) 

    add_test(run_simple_python_example
	    ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_python_example
                 -D TEST_PROGRAM=${Python_EXECUTABLE}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/python/example.py 
		             -D TEST_INPUT=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/python/run_simple_python_example.log
		             -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/python/run_simple_python_example.log
                 -P ${CMAKE_SOURCE_DIR}/tests/python/run_python_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/python/) 
endif()
