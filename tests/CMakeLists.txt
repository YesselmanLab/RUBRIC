# RUBRIC test suite. Build and run from the build directory:
#
#   cmake -S . -B build && cmake --build build -j && ctest --test-dir build --output-on-failure
#
# Each test program runs offline against the bundled resources/ and the small
# structures in tests/data/. Labels: unit, cli, regression (ctest -L unit).

set(RUBRIC_TEST_ENV "RUBRIC_RESOURCES=${PROJECT_SOURCE_DIR}/resources")

function(rubric_add_test name label)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE pairfinder_core)
    target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(${name} PRIVATE
        RUBRIC_TESTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
        RUBRIC_EXE="$<TARGET_FILE:rubric>"
        RUBRIC_VERSION_EXPECTED="${PROJECT_VERSION}")
    if(MSVC)
        target_compile_options(${name} PRIVATE /W4 $<$<BOOL:${PAIRFINDER_WARNINGS_AS_ERRORS}>:/WX>)
    else()
        target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic
            $<$<BOOL:${PAIRFINDER_WARNINGS_AS_ERRORS}>:-Werror>)
    endif()
    add_dependencies(${name} rubric)
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES LABELS ${label} ENVIRONMENT "${RUBRIC_TEST_ENV}" TIMEOUT 120)
endfunction()

rubric_add_test(test_geometry   unit unit/test_geometry.cpp)
rubric_add_test(test_suites     unit unit/test_suites.cpp)
rubric_add_test(test_backbone   unit unit/test_backbone.cpp)
rubric_add_test(test_basepairs  unit unit/test_basepairs.cpp)
rubric_add_test(test_chemistry  unit unit/test_chemistry.cpp)
rubric_add_test(test_cli        cli  cli/test_cli.cpp)
rubric_add_test(test_inputs     cli  cli/test_inputs.cpp)
rubric_add_test(test_panel      regression regression/test_panel.cpp)
