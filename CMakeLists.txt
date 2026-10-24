cmake_minimum_required(VERSION 3.28)  # FetchContent EXCLUDE_FROM_ALL/SYSTEM
project(pairfinder
    VERSION 0.1.0
    DESCRIPTION "Fast RNA base-pair finder + scorer"
    LANGUAGES CXX)

# ---- Standard / options -----------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PAIRFINDER_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)

# ---- Core library -----------------------------------------------------------
add_library(pairfinder_core
    src/core/atom_type.cpp
    src/core/residue.cpp
    src/core/nucleotide_registry.cpp
    src/core/base_typing.cpp
    src/core/resource_locator.cpp
    src/geometry/superposition.cpp
    src/algorithms/base_frame_calculator.cpp
    src/algorithms/candidate_finder.cpp
    src/algorithms/edge_classifier.cpp
    src/algorithms/saenger.cpp
    src/algorithms/template_aligner.cpp
    src/algorithms/hbond_patterns.cpp
    src/algorithms/confidence.cpp
    src/algorithms/quality_scorer.cpp
    src/algorithms/classify_pair.cpp
    src/algorithms/carbonyl_margin.cpp
    src/algorithms/score_candidates.cpp
    src/algorithms/rigid_body.cpp
    src/algorithms/rna_chains.cpp
    src/algorithms/selection.cpp
    src/algorithms/hbond/chemistry.cpp
    src/algorithms/hbond/finder.cpp
    src/scoring/pair_parameters.cpp
    src/scoring/torsions.cpp
    src/scoring/richardson.cpp
    src/scoring/scorer.cpp
    src/scoring/clashes.cpp
    src/io/pdb_parser.cpp
    src/pipeline.cpp
)
add_library(pairfinder::core ALIAS pairfinder_core)
target_include_directories(pairfinder_core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/third_party>)
target_compile_features(pairfinder_core PUBLIC cxx_std_20)
# Dev fallback for the resource locator: the in-source resources/ dir, used only
# when the locator can't find resources exe-relative, cwd-relative, or via
# $PARSE_RESOURCES (so running from a dev tree always works).
target_compile_definitions(pairfinder_core PRIVATE
    PAIRFINDER_SOURCE_RESOURCES="${CMAKE_CURRENT_SOURCE_DIR}/resources")

# ---- GEMMI (PDB/CIF/mmCIF parsing) ------------------------------------------
# Self-contained: fetch + build GEMMI v0.6.5 IN-TREE (own libgemmi_cpp under our
# build/) by cloning from upstream git.
include(FetchContent)
FetchContent_Declare(gemmi
    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
    GIT_TAG v0.6.5
    EXCLUDE_FROM_ALL    # only the gemmi targets we link get built
    SYSTEM)             # treat gemmi headers as system (no warnings from them)
# GEMMI's CMake unconditionally calls add_test(cpptest), but its executable is
# EXCLUDE_FROM_ALL (never built) -> it would pollute our suite as a missing-exe
# failure. Intercept add_test() to drop GEMMI's cpptest while forwarding every
# other registration (our own tests) to the builtin _add_test.
function(add_test)
    if("${ARGV1}" STREQUAL "cpptest")
        return()
    endif()
    _add_test(${ARGV})
endfunction()
# Build GEMMI as a static lib.
set(_pf_shared_saved ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "" FORCE)
FetchContent_MakeAvailable(gemmi)
set(BUILD_SHARED_LIBS ${_pf_shared_saved} CACHE INTERNAL "" FORCE)
# gemmi_cpp propagates its own include dirs (and zlib) transitively.
target_link_libraries(pairfinder_core PUBLIC gemmi_cpp)
target_compile_options(pairfinder_core PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<BOOL:${PAIRFINDER_WARNINGS_AS_ERRORS}>:-Werror>)

# ---- CLI app ----------------------------------------------------------------
add_executable(parse
    apps/main.cpp
    apps/support/cli_support.cpp
    apps/differential/differential_commands.cpp)
target_link_libraries(parse PRIVATE pairfinder_core)
# apps/ on the include path so the CLI TUs include their siblings as
# "support/cli_support.hpp" / "differential/differential_commands.hpp" (no ../ paths).
target_include_directories(parse PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/apps)

# After each build: (1) regenerate the precomputed resource caches (H-bond
# patterns from the idealized templates, plus base-typing + H-bond-chemistry
# tables from the ligand DB) so they can never go stale, then (2) copy resources/
# next to the binary so the exe-relative locator finds it (build/<cfg>/resources).
# Step 1 writes into the source tree (output is deterministic, so it only changes
# when the underlying data does) and step 2 stages the fresh copies for the build.
add_custom_command(TARGET parse POST_BUILD
    COMMAND parse dump-caches "${CMAKE_CURRENT_SOURCE_DIR}/resources"
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/resources" "$<TARGET_FILE_DIR:parse>/resources"
    COMMENT "Regenerating resource caches + staging bundled resources next to parse")

# ---- Install ----------------------------------------------------------------
# Installed layout: <prefix>/bin/parse + <prefix>/share/parse/resources
# (resolved by the locator's <exe>/../share/parse/resources rule).
install(TARGETS parse RUNTIME DESTINATION bin)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/resources/"
    DESTINATION share/parse/resources)

# Optional PyMOL viewer: the launcher goes on PATH (bin/), its plugin alongside
# the shared data. parse-view resolves the plugin from ../share/parse and
# finds `parse` on PATH when installed.
install(PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/integrations/pymol/parse-view"
    DESTINATION bin)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/integrations/pymol/parse_pymol.py"
    DESTINATION share/parse)

# Unit tests are developer-only and are not part of this package.
