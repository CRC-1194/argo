# wrapper around cmake's add_test
function(add_argo_test)
  include(CMakeParseArguments)
  set(OPTIONS)
  set(SINGLEARGS NAME TARGET TIMEOUT)
  set(MULTIARGS SOURCES COMPILE_DEFINITIONS COMPILE_FLAGS CMD_ARGS COMMAND CONDITIONS LINK_LIBRARIES)
  cmake_parse_arguments(ADDARGOTEST "${OPTIONS}" "${SINGLEARGS}" "${MULTIARGS}" ${ARGN})

  if (ADDARGOTEST_UNPARSED_ARGUMENTS)
      message(WARNING "Some arguments have not been parsed in call to add_argo_test: '${ADDARGOTEST_UNPARSED_ARGUMENTS}'!")
  endif ()

  # check input arguments validity
  if (NOT ADDARGOTEST_SOURCES AND NOT ADDARGOTEST_TARGET)
      message(FATAL_ERROR "Either the SOURCES or the TARGET option must be specied when adding a test")
  endif ()

  if (ADDARGOTEST_SOURCES AND ADDARGOTEST_TARGET)
      message(FATAL_ERROR "Cannot specify both SOURCES and TARGET when adding a test")
  endif()

  if (NOT ADDARGOTEST_NAME)
      message(FATAL_ERROR "No name given for the new test")
  endif ()

  # set default values
  if (NOT ADDARGOTEST_COMMAND)
      set(ADDARGOTEST_COMMAND ${ADDARGOTEST_NAME})
  endif ()

  if (NOT ADDARGOTEST_TIMEOUT)
      set(ADDARGOTEST_TIMEOUT 300)
  endif ()

  # Check conditions, skip test if not fulfilled
  set(CONDS_FULFILLED true)
  foreach (cond IN LISTS ADDARGOTEST_CONDITIONS)
    if (NOT cond)
        set(CONDS_FULFILLED false)
    endif ()
  endforeach ()

  if (ADDARGOTEST_SOURCES AND CONDS_FULFILLED)
    add_executable(${ADDARGOTEST_NAME} ${ADDARGOTEST_SOURCES})
    target_compile_definitions(
        ${ADDARGOTEST_NAME}
        PUBLIC ${ADDARGOTEST_COMPILE_DEFINITIONS}
    )

    if (ADDARGOTEST_COMPILE_FLAGS)
        target_compile_options(
            ${ADDARGOTEST_NAME}
            PUBLIC ${ADDARGOTEST_COMPILE_FLAGS}
        )
    endif ()

    if (ADDARGOTEST_LINK_LIBRARIES)
        target_link_libraries(${ADDARGOTEST_NAME} ${ADDARGOTEST_LINK_LIBRARIES})
    endif ()

    set(ADDARGOTEST_TARGET ${ADDARGOTEST_NAME})

  elseif (NOT CONDS_FULFILLED)
    message(STATUS "Test ${ADDARGOTEST_NAME} is disabled due to unmet requirements")

    # write dummy main file with message that test is skipped,
    # returning the error code that we set as property later
    set(DUMMY_MAIN "${CMAKE_CURRENT_BINARY_DIR}/test_dummy_${ADDARGOTEST_NAME}.cc")
    file(WRITE ${DUMMY_MAIN}
               "#include <iostream>\n\n"
               "int main()\n"
               "{\n"
               "    std::cout << \"Test is skipped due to unmet condition\" << std::endl;\n"
               "    return 255;\n"
               "}")
    add_executable(${ADDARGOTEST_NAME} ${DUMMY_MAIN})
    set(ADDARGOTEST_TARGET ${ADDARGOTEST_NAME})
    set(ADDARGOTEST_COMMAND ${ADDARGOTEST_NAME})
    set(ADDARGOTEST_CMD_ARGS "")
  endif ()

  # Now add the actual test
  add_test(NAME ${ADDARGOTEST_NAME} COMMAND "${ADDARGOTEST_COMMAND}" ${ADDARGOTEST_CMD_ARGS})
  set_tests_properties(${ADDARGOTEST_NAME} PROPERTIES SKIP_RETURN_CODE 255)
  set_tests_properties(${ADDARGOTEST_NAME} PROPERTIES TIMEOUT ${ADDARGOTEST_TIMEOUT})
endfunction(add_argo_test)
