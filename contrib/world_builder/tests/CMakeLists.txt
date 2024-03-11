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
          "Please specify TEST_DIFF by hand.")
ENDIF()

# Unit test testing
#find all the unit test files
SET(WB_INCLUDE_UNIT_TEST ON
    CACHE BOOL "Whether to include unit tests or not in the build.")

SET(WB_RUN_TESTS_WITH_GDB OFF
    CACHE BOOL "Whether to also run the test with GDB")

SET(WB_RUN_TESTS_WITH_VALGRIND OFF
    CACHE BOOL "Whether to also run the test with Valgrind")

if(WB_INCLUDE_UNIT_TEST)
  file(GLOB_RECURSE UNIT_TEST_SOURCES_CXX "unit_tests/*.cc")

  if(CMAKE_Fortran_COMPILER)
    file(GLOB_RECURSE UNIT_TEST_SOURCES_FORTRAN "unit_tests/*.f90")
  endif()

  set(UNIT_TEST_SOURCES ${UNIT_TEST_SOURCES_CXX} ${UNIT_TEST_SOURCES_FORTRAN})

  set(test_name "unit_tests")
  
  # Add compile target
  add_executable(${test_name} ${UNIT_TEST_SOURCES})
  
  if(NOT ${CMAKE_VERSION} VERSION_LESS "3.9.0") # Preventing issues with older cmake compilers which do not support VERSION_GREATER_EQUAL
    target_compile_options(${test_name} INTERFACE ${WB_COMPILER_OPTIONS_INTERFACE} PRIVATE ${WB_COMPILER_OPTIONS_PRIVATE})
  endif()

  if(${CMAKE_VERSION} VERSION_LESS "3.13.0")
    SET(CMAKE_EXE_LINKER_FLAGS "${WB_LINKER_OPTIONS}")
  else()
    target_link_options(${test_name} INTERFACE ${WB_LINKER_OPTIONS})
  endif()

  # Make sure that the whole library is loaded, so the registration is done correctly.
  if(NOT APPLE)
    SET(GWB_LIBRARY_WHOLE -Wl,--whole-archive ${WB_TARGET} -Wl,--no-whole-archive)
  elseif(MSVC)
    SET(GWB_LIBRARY_WHOLE ${WB_TARGET})
    SET(CMAKE_EXE_LINKER_FLAGS "/WHOLEARCHIVE:lib${WB_TARGET}")
  else()
    SET(GWB_LIBRARY_WHOLE -Wl,-force_load ${WB_TARGET})
  endif()

  if(${USE_MPI})
    target_link_libraries(${test_name} PUBLIC MPI::MPI_CXX ${GWB_LIBRARY_WHOLE})
  else()
    target_link_libraries(${test_name} ${GWB_LIBRARY_WHOLE})
  endif()

  #Move testing binaries into a testBin directory
  set_target_properties(${test_name} PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/bin/)

  #Finally add it to test execution -
  #Notice the WORKING_DIRECTORY and COMMAND
  add_test(NAME ${test_name}
          WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
          COMMAND ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name})

  if(WB_RUN_TESTS_WITH_GDB)
    add_test(NAME ${test_name}_gdb
            WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
            COMMAND gdb --ex r --ex bt -ex "set confirm off" --ex exit($_exitcode) --args ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name})
  endif()

  if(WB_RUN_TESTS_WITH_VALGRIND)
    add_test(NAME ${test_name}_valgrind
             WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
             COMMAND valgrind -v --leak-check=full --error-exitcode=1  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${test_name})
  endif()
endif()

# App testing
IF(WB_ENABLE_APPS AND WB_RUN_APP_TESTS)
  # Create directory for the test results
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/gwb-dat)

  # Test help
  add_test(testing_help
           ${CMAKE_COMMAND}
           -D TEST_NAME=testing_help
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=--help
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/testing_help/screen-output.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/testing_help/screen-output.log
           -D TEST_DIFF=${TEST_DIFF}
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)
  
  # Test no files provided
  add_test(testing_no_file
          ${CMAKE_COMMAND}
          -D TEST_NAME=testing_no_file
          -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
          -D TEST_ARGS=
          -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/testing_no_file/screen-output.log
          -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/testing_no_file/screen-output.log
          -D TEST_DIFF=${TEST_DIFF}
          -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)
  
  
  # Test one file provided 
  add_test(testing_one_file
           ${CMAKE_COMMAND}
           -D TEST_NAME=testing_one_file
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=non_existent_file
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/testing_one_file/screen-output.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/testing_one_file/screen-output.log
           -D TEST_DIFF=${TEST_DIFF}
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)


  # Test two many arguments provided 
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.dat\;lalala\;nonono")
  add_test(testing_too_many_arguments
           ${CMAKE_COMMAND}
           -D TEST_NAME=testing_too_many_arguments
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/testing_too_many_arguments/screen-output.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/testing_too_many_arguments/screen-output.log
           -D TEST_DIFF=${TEST_DIFF}
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)

  # Add tests which generates the declarations and schema files in the docs folder. Do not care about the output.
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/doc/generate_decl_schema.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-dat/app_wb2.dat\;--output-json-files")
  add_test(generate_declarations_and_schema
           ${CMAKE_COMMAND}
           -D TEST_NAME=generate_declarations_and_schema
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/generate_declarations_and_schema/screen-output.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/generate_declarations_and_schema/screen-output.log
           -D TEST_DIFF=${TEST_DIFF}
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)


  # Add tests which tests the help output
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_filtered.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid\;--help")
  add_test(grid_help
           ${CMAKE_COMMAND}
	         -D TEST_NAME=${test_name}
	         -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
	         -D TEST_ARGS=${TEST_ARGUMENTS}
	         -D TEST_DIFF=${TEST_DIFF}
	         -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/grid_help.log
	         -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/grid_help.log
	         -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
	         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)


  # Add tests which test what happens when too many arguments are provided.
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_filtered.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid\;too;many;arguments")
  add_test(grid_too_many_arguments
           ${CMAKE_COMMAND}
           -D TEST_NAME=${test_name}
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_DIFF=${TEST_DIFF}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/grid_too_many_arguments.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/grid_too_many_arguments.log
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)


  # Add tests which test the grid version output
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_filtered.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid\;--version")
  add_test(grid_version
           ${CMAKE_COMMAND}
           -D TEST_NAME=${test_name}
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_DIFF=${TEST_DIFF}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/grid_version.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/grid_version.log
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)



  # Add tests which input a non-existing input file
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_non-existing.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid")
  add_test(grid_non_existing_wb
           ${CMAKE_COMMAND}
           -D TEST_NAME=${test_name}
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_DIFF=${TEST_DIFF}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/grid_non_existing_wb.log
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/grid_non_existing_wb.log
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)

  # Add tests which filters the output.
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_filtered.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid\;--filtered")
  add_test(grid_filtered  
           ${CMAKE_COMMAND}
	         -D TEST_NAME=${test_name}
	         -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
	         -D TEST_ARGS=${TEST_ARGUMENTS}
	         -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/subducting_plate_composition_smooth_filtered.filtered.vtu
	         -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/subducting_plate_composition_smooth_filtered.filtered.vtu
	         -P ${CMAKE_SOURCE_DIR}/tests/gwb-grid/run_gwb-grid_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)


  # Add tests which filters the by tag.
  set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/data/subducting_plate_composition_smooth_by_tag.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/subducting_plate_composition_smooth.grid\;--by-tag")
  add_test(grid_by_tag  
           ${CMAKE_COMMAND}
           -D TEST_NAME=${test_name}
           -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
           -D TEST_ARGS=${TEST_ARGUMENTS}
           -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/subducting_plate_composition_smooth_by_tag.1.vtu
           -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/subducting_plate_composition_smooth_by_tag.1.vtu
           -P ${CMAKE_SOURCE_DIR}/tests/gwb-grid/run_gwb-grid_tests.cmake
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)
  
  #find all the integration test files
  file(GLOB_RECURSE APP_TEST_SOURCES "gwb-dat/*.wb")
  
  # Run through each sourceUforeach(test_source ${APP_TEST_SOURCES})
  foreach(test_source ${APP_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.dat\;--limit-debug-consistency-checks")
    add_test(dat_${test_name}
             ${CMAKE_COMMAND}
             -D TEST_NAME=${test_name}
             -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}
             -D TEST_ARGS=${TEST_ARGUMENTS}
             -D TEST_DAT=${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.dat
             -D TEST_SOURCE_DIR=${CMAKE_SOURCE_DIR}
             -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-dat/${test_name}/screen-output.log
             -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/${test_name}/screen-output.log
             -D TEST_DIFF=${TEST_DIFF}
             -P ${CMAKE_SOURCE_DIR}/tests/gwb-dat/run_gwb-dat_tests.cmake
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-dat/)

    if(WB_RUN_TESTS_WITH_GDB)
      add_test(NAME dat_${test_name}_gdb
              WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
              COMMAND gdb --return-child-result --ex "set confirm off" --ex r  --ex bt --ex exit  --args ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX} ${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.wb ${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.dat --limit-debug-consistency-checks)
    endif()

    if(WB_RUN_TESTS_WITH_VALGRIND)
      add_test(NAME dat_${test_name}_valgrind
              WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
              COMMAND valgrind -v --leak-check=full --error-exitcode=1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-dat${CMAKE_EXECUTABLE_SUFFIX}  ${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.wb ${CMAKE_SOURCE_DIR}/tests/gwb-dat/${test_name}.dat --limit-debug-consistency-checks)
    endif()
  endforeach(test_source)
ENDIF()

# gwb-grid tests
if (WB_ENABLE_APPS)
  # Create directory for the test results
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/gwb-grid)
  # find all the integration visu files
  file(GLOB_RECURSE VISU_TEST_SOURCES "gwb-grid/*.wb")
  # Run through each sourceUforeach(test_source ${VISU_TEST_SOURCES})
  foreach(test_source ${VISU_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.wb\;${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.grid\;-j\;1")
    add_test(grid_${test_name}
            ${CMAKE_COMMAND}
            -D TEST_NAME=${test_name}
            -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
            -D TEST_ARGS=${TEST_ARGUMENTS}
            -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/gwb-grid/${test_name}.vtu
            -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/${test_name}.vtu
            -P ${CMAKE_SOURCE_DIR}/tests/gwb-grid/run_gwb-grid_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/gwb-grid/)

    if(WB_RUN_TESTS_WITH_GDB)
      add_test(NAME grid_${test_name}_gdb
              WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
              COMMAND gdb --return-child-result --ex "set confirm off" --ex r  --ex bt --ex exit  --args ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX} ${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.wb ${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.grid)
    endif()

    if(WB_RUN_TESTS_WITH_VALGRIND)
      add_test(NAME grid_${test_name}_valgrind
              WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/
              COMMAND valgrind -v --leak-check=full --error-exitcode=1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}  ${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.wb ${CMAKE_SOURCE_DIR}/tests/gwb-grid/${test_name}.grid)
    endif()
  endforeach(test_source)
endif()

add_custom_target(update_test_references 
                  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                  COMMAND cp -r ${CMAKE_BINARY_DIR}/tests/ ${CMAKE_SOURCE_DIR}/
                  COMMENT "Updating the World Builder test results...")

# Documentation tests
if (WB_ENABLE_APPS)
  # Create directory for the test results
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/documentation)
  #find all the integration visu files
  file(GLOB_RECURSE DOCU_TEST_SOURCES "${CMAKE_SOURCE_DIR}/doc/sphinx/_static/gwb_input_files/*.wb")
  # Run through each source and create a test for each one
  foreach(test_source ${DOCU_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    set(TEST_ARGUMENTS "${CMAKE_SOURCE_DIR}/doc/sphinx/_static/gwb_input_files/${test_name}.wb\;${CMAKE_SOURCE_DIR}/doc/sphinx/_static/gwb_input_files/${test_name}.grid")
    add_test(doc_${test_name}
            ${CMAKE_COMMAND}
            -D TEST_NAME=${test_name}
            -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
            -D TEST_ARGS=${TEST_ARGUMENTS}
            -P ${CMAKE_SOURCE_DIR}/tests/documentation/run_documentation_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation/)
  endforeach(test_source)
 endif()

# Cookbook tests
if (WB_ENABLE_APPS)
  # Create directory for the test results
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/cookbooks)
  #find all the integration visu files
  file(GLOB_RECURSE COOKBOOK_TEST_SOURCES "${CMAKE_SOURCE_DIR}/cookbooks/*/*.wb")
  # Run through each source and add a test for each one
  foreach(test_source ${COOKBOOK_TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    get_filename_component(test_dir ${test_source} DIRECTORY)
    set(TEST_ARGUMENTS "${test_dir}/${test_name}.wb\;${test_dir}/${test_name}.grid\;--resolution-limit\;50")
    add_test(cookbooks_${test_name}
            ${CMAKE_COMMAND}
            -D TEST_NAME=${test_name}
            -D TEST_PROGRAM=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/gwb-grid${CMAKE_EXECUTABLE_SUFFIX}
            -D TEST_ARGS=${TEST_ARGUMENTS}
            -P ${CMAKE_SOURCE_DIR}/tests/documentation/run_documentation_tests.cmake
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/cookbooks/)
  endforeach(test_source)
endif()

#test C compilation and wrapper if compiler found
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C)
  if(NOT MSVC AND NOT APPLE)
    add_test(NAME compile_simple_C_test
             COMMAND ${CMAKE_C_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/C/test.c -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -lm  ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
  elseif(APPLE)
    add_test(NAME compile_simple_C_test
	    COMMAND ${CMAKE_C_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/C/test.c -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -lm -I/usr/include/ -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/ -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)  
  else()
    #MSVS
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      add_test(NAME compile_simple_C_test 
        COMMAND ${CMAKE_C_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/C/test.c /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
    else()
      add_test(NAME compile_simple_C_test 
      COMMAND ${CMAKE_C_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/C/test.c /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX}  ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
    endif()
  endif()
  set_tests_properties(compile_simple_C_test PROPERTIES TIMEOUT 60) 
  add_test(run_simple_C_test
           ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_C_test
                 -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/C/test${CMAKE_EXECUTABLE_SUFFIX}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/C/run_simple_C_test.log
		 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/C/run_simple_C_test.log
                 -P ${CMAKE_SOURCE_DIR}/tests/C/run_C_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/) 
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C)
  set_tests_properties(run_simple_C_test PROPERTIES TIMEOUT 60) 
  set_tests_properties(run_simple_C_test PROPERTIES DEPENDS compile_simple_C_test)

  if(NOT MSVC AND NOT APPLE)
    add_test(NAME compile_simple_C_example 
	     COMMAND ${CMAKE_C_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/C/example.c -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -lm  -L${MPI_C_LIBRARIES} -L${MPI_CXX_LIBRARIES}
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
  elseif(APPLE)
    add_test(NAME compile_simple_C_example 
	     COMMAND ${CMAKE_C_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/C/example.c -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -lm -I/usr/include/ -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/ -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
  else()
    #MSVS
    if(CMAKE_BUILD_TYPE STREQUAL Debug)
      add_test(NAME compile_simple_C_example 
      COMMAND ${CMAKE_C_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/C/example.c /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} -L${MPI_C_LIBRARIES} -L${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
    else()
      add_test(NAME compile_simple_C_example 
      COMMAND ${CMAKE_C_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/C/example.c /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/)
    endif()
  endif()
  set_tests_properties(compile_simple_C_example PROPERTIES TIMEOUT 60) 
  add_test(run_simple_C_example
	   ${CMAKE_COMMAND}
                 -D TEST_NAME=run_simple_C_example
                 -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/C/example${CMAKE_EXECUTABLE_SUFFIX}
                 -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
                 -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/C/run_simple_C_example.log
		 -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/C/run_simple_C_example.log
                 -P ${CMAKE_SOURCE_DIR}/tests/C/run_C_tests.cmake
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/C/) 
  set_tests_properties(run_simple_C_example PROPERTIES TIMEOUT 60) 
  set_tests_properties(run_simple_C_example PROPERTIES DEPENDS compile_simple_C_example)


#test CPP compilation and wrapper
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP)
list(GET 1 ${MPI_INCLUDE_PATH}, MPI_single_include_path)
if(NOT MSVC AND NOT APPLE)
  add_test(NAME compile_simple_CPP_test
  COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP/test.cpp -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -std=c++14 -lm -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
elseif(APPLE)
  add_test(NAME compile_simple_CPP_test
    COMMAND clang ${CMAKE_CURRENT_SOURCE_DIR}/CPP/test.cpp -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -std=c++14 -lm -I/usr/include/ -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/  ${LDFLAGS} -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)  
else()
  #MSVS
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_test(NAME compile_simple_CPP_test 
      COMMAND ${CMAKE_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP/test.cpp /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
  else()
    add_test(NAME compile_simple_CPP_test 
    COMMAND ${CMAKE_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP/test.cpp /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
  endif()
endif()
set_tests_properties(compile_simple_CPP_test PROPERTIES TIMEOUT 60) 
add_test(run_simple_CPP_test
         ${CMAKE_COMMAND}
               -D TEST_NAME=run_simple_CPP_test
               -D MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
               -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/CPP/test${CMAKE_EXECUTABLE_SUFFIX} 
               -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
               -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/CPP/run_simple_CPP_test.log
   -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/CPP/run_simple_CPP_test.log
               -P ${CMAKE_SOURCE_DIR}/tests/CPP/run_CPP_tests.cmake
               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP)
set_tests_properties(run_simple_CPP_test PROPERTIES TIMEOUT 60)
set_tests_properties(run_simple_CPP_test PROPERTIES DEPENDS compile_simple_CPP_test)

if(NOT MSVC AND NOT APPLE)
  add_test(NAME compile_simple_CPP_example 
     COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP/example.cpp -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -std=c++14 -lm  -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
elseif(APPLE)
  add_test(NAME compile_simple_CPP_example 
     COMMAND clang ${CMAKE_CURRENT_SOURCE_DIR}/CPP/example.cpp -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -std=c++14 -lm -I/usr/include/ -I/Library/Developer/CPPommandLineTools/SDKs/MacOSX.sdk/usr/include/ -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
else()
  #MSVS
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_test(NAME compile_simple_CPP_example 
    COMMAND ${CMAKE_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP/example.cpp /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
  else()
    add_test(NAME compile_simple_CPP_example 
    COMMAND ${CMAKE_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP/example.cpp /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/)
  endif()
endif()
set_tests_properties(compile_simple_CPP_example PROPERTIES TIMEOUT 60) 
add_test(run_simple_CPP_example
   ${CMAKE_COMMAND}
               -D TEST_NAME=run_simple_CPP_example
               -D MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
               -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/CPP/example${CMAKE_EXECUTABLE_SUFFIX}
               -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
               -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/CPP/run_simple_CPP_example.log
   -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/CPP/run_simple_CPP_example.log
               -P ${CMAKE_SOURCE_DIR}/tests/CPP/run_CPP_tests.cmake
               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP/) 
set_tests_properties(run_simple_CPP_example PROPERTIES TIMEOUT 60) 
set_tests_properties(run_simple_CPP_example PROPERTIES DEPENDS compile_simple_CPP_example)

#test CPP MPI compilation and wrapper if MPI found
if(USE_MPI)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI)
list(GET 1 ${MPI_INCLUDE_PATH}, MPI_single_include_path)
if(NOT MSVC AND NOT APPLE)
  add_test(NAME compile_simple_CPP_MPI_test
  COMMAND ${MPI_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/test.cpp -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -std=c++14 -lm -L${MPI_C_LIBRARIES} -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
elseif(APPLE)
  add_test(NAME compile_simple_CPP_MPI_test
    COMMAND ${MPI_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/test.cpp -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I../../inlcude/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -std=c++14 -lm -I/usr/include/ -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include/  ${LDFLAGS} -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib  -L${MPI_C_LIBRARIES} -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)  
else()
  #MSVS
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_test(NAME compile_simple_CPP_MPI_test 
      COMMAND ${MPI_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/test.cpp /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
  else()
    add_test(NAME compile_simple_CPP_MPI_test 
    COMMAND ${MPI_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/test.cpp /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:test${CMAKE_EXECUTABLE_SUFFIX}  ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
  endif()
endif()
set_tests_properties(compile_simple_CPP_MPI_test PROPERTIES TIMEOUT 60) 
add_test(run_simple_CPP_MPI_test
         ${CMAKE_COMMAND}
               -D TEST_NAME=run_simple_CPP_MPI_test
               -D MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
               -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/CPP_MPI/test${CMAKE_EXECUTABLE_SUFFIX} 
               -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
               -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/CPP_MPI/run_simple_CPP_MPI_test.log
   -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/run_simple_CPP_MPI_test.log
               -P ${CMAKE_SOURCE_DIR}/tests/CPP_MPI/run_CPP_MPI_tests.cmake
               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/) 
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI)
set_tests_properties(run_simple_CPP_MPI_test PROPERTIES TIMEOUT 60)
set_tests_properties(run_simple_CPP_MPI_test PROPERTIES DEPENDS compile_simple_CPP_MPI_test)

if(NOT MSVC AND NOT APPLE)
  add_test(NAME compile_simple_CPP_MPI_example 
     COMMAND ${MPI_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/example.cpp -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lstdc++ -std=c++14 -lm  -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
elseif(APPLE)
  add_test(NAME compile_simple_CPP_MPI_example 
     COMMAND ${MPI_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/example.cpp -L${CMAKE_BINARY_DIR}/lib/ -Wl,-force_load,${CMAKE_BINARY_DIR}/lib/lib${WB_TARGET}.a -I${CMAKE_SOURCE_DIR}/include/ ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_COMPILER_OPTIONS_PRIVATE_COVERAGE_NEW} -lc++ -std=c++14 -lm -I/usr/include/ -I/Library/Developer/CPP_MPIommandLineTools/SDKs/MacOSX.sdk/usr/include/  -L${MPI_CXX_LIBRARIES} -I${MPI_single_include_path}
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
else()
  #MSVS
  if(CMAKE_BUILD_TYPE STREQUAL Debug)
    add_test(NAME compile_simple_CPP_MPI_example 
    COMMAND ${MPI_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/example.cpp /W3 /EHsc /RTC1 /MDd /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} ${MPI_CXX_LIBRARIES} 
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
  else()
    add_test(NAME compile_simple_CPP_MPI_example 
    COMMAND ${MPI_CXX_COMPILER} /I${CMAKE_SOURCE_DIR}\\include ${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/example.cpp /W3 /EHsc /RTC1 /MD /link /INCREMENTAL /NOLOGO ${CMAKE_BINARY_DIR}\\lib\\${WB_TARGET}.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib /IMPLIB:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib /WHOLEARCHIVE:${CMAKE_BINARY_DIR}/lib/${WB_TARGET}.lib ${WB_C_COMPILER_FLAGS_COVERAGE} /out:example${CMAKE_EXECUTABLE_SUFFIX} ${MPI_CXX_LIBRARIES}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/)
  endif()
endif()
set_tests_properties(compile_simple_CPP_MPI_example PROPERTIES TIMEOUT 60) 
add_test(run_simple_CPP_MPI_example
   ${CMAKE_COMMAND}
               -D TEST_NAME=run_simple_CPP_MPI_example
               -D MPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
               -D TEST_PROGRAM=${CMAKE_BINARY_DIR}/tests/CPP_MPI/example${CMAKE_EXECUTABLE_SUFFIX}
               -D TEST_ARGS=${CMAKE_CURRENT_SOURCE_DIR}/data/continental_plate.wb
               -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/CPP_MPI/run_simple_CPP_MPI_example.log
   -D TEST_REFERENCE=${CMAKE_CURRENT_SOURCE_DIR}/CPP_MPI/run_simple_CPP_MPI_example.log
               -P ${CMAKE_SOURCE_DIR}/tests/CPP_MPI/run_CPP_MPI_tests.cmake
               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/CPP_MPI/) 
set_tests_properties(run_simple_CPP_MPI_example PROPERTIES TIMEOUT 60) 
set_tests_properties(run_simple_CPP_MPI_example PROPERTIES DEPENDS compile_simple_CPP_MPI_example)
endif()

#test fortran compilation and wrapper if compiler found
if(CMAKE_Fortran_COMPILER)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran)
  if(NOT APPLE)
    add_test(NAME compile_simple_fortran_test
             COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test.f90 -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lstdc++ ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} 
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  else()
    add_test(NAME compile_simple_fortran_test
	    COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test.f90 -L../../lib/ -Wl,-force_load,../../lib/lib${WB_TARGET}.a -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o test${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lc++ ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} 
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
  set_tests_properties(run_simple_fortran_test PROPERTIES DEPENDS compile_simple_fortran_test)

  if(NOT APPLE)
    add_test(NAME compile_simple_fortran_example 
	     COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/example.f90 -L../../lib/ -Wl,--whole-archive -l${WB_TARGET} -Wl,--no-whole-archive -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lstdc++ ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} 
	     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/fortran/)
  else()
    add_test(NAME compile_simple_fortran_example 
	     COMMAND ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/fortran/example.f90 -L../../lib/ -Wl,-force_load,../../lib/lib${WB_TARGET}.a -I../../mod/ ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -o example${CMAKE_EXECUTABLE_SUFFIX} ${WB_FORTRAN_COMPILER_FLAGS_COVERAGE} -lc++ ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES} 
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
  set_tests_properties(run_simple_fortran_example PROPERTIES DEPENDS compile_simple_fortran_example)
 endif()

#test python compilation and wrapper if compiler found
if(MAKE_PYTHON_WRAPPER)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests/python)


  add_test(run_python_version_test
  ${CMAKE_COMMAND}
      -D TEST_NAME=run_python_version_test
      -D TEST_PROGRAM=${Python_EXECUTABLE}
      -D TEST_ARGS=--version
      -D TEST_INPUT=
      -D TEST_OUTPUT=${CMAKE_BINARY_DIR}/tests/python/run_python_version_test.log
      -D TEST_REFERENCE=${CMAKE_BINARY_DIR}/tests/python/run_python_version_test.log # the output is not important
      -D TEST_NO_DIFF=TRUE
      -P ${CMAKE_SOURCE_DIR}/tests/python/run_python_tests.cmake
     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/python/)   


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
