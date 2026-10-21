cmake_minimum_required(VERSION 3.16)
project(al_maleh_vector_efg LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ZLIB REQUIRED)
execute_process(COMMAND git rev-parse HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  OUTPUT_VARIABLE AL_MALEH_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT AL_MALEH_GIT_COMMIT)
  set(AL_MALEH_GIT_COMMIT "unknown")
endif()
add_library(al_maleh_mls STATIC src/mls.cpp)
target_include_directories(al_maleh_mls PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(al_maleh_mls PUBLIC Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_mls PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(mls_probe apps/mls_probe.cpp)
target_link_libraries(mls_probe PRIVATE al_maleh_mls)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_mls tests/test_mls.cpp)
  target_link_libraries(test_mls PRIVATE al_maleh_mls)
  add_test(NAME mls_consistency COMMAND test_mls)
endif()

# M2: vector Cartesian EFG with a background integration grid.
add_library(al_maleh_cartesian STATIC src/cartesian.cpp)
target_include_directories(al_maleh_cartesian PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(al_maleh_cartesian PUBLIC al_maleh_mls Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_cartesian PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(cartesian_rect apps/cartesian_rect.cpp)
target_link_libraries(cartesian_rect PRIVATE al_maleh_cartesian)
if(BUILD_TESTING)
  add_executable(test_cartesian tests/test_cartesian.cpp)
  target_link_libraries(test_cartesian PRIVATE al_maleh_cartesian)
  add_test(NAME cartesian_assembly COMMAND test_cartesian)
endif()

# M3: generalized symmetric spectral solver (small dense reference only).
add_library(al_maleh_spectral STATIC src/spectral.cpp src/spectral_export.cpp)
target_include_directories(al_maleh_spectral PUBLIC include)
target_link_libraries(al_maleh_spectral PUBLIC al_maleh_cartesian Eigen3::Eigen ZLIB::ZLIB)
target_compile_definitions(al_maleh_spectral PRIVATE
  AL_MALEH_GIT_COMMIT="${AL_MALEH_GIT_COMMIT}"
  AL_MALEH_BUILD_TYPE="${CMAKE_BUILD_TYPE}"
  AL_MALEH_COMPILER="${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_spectral PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(cartesian_spectrum apps/cartesian_spectrum.cpp)
target_link_libraries(cartesian_spectrum PRIVATE al_maleh_spectral)
if(BUILD_TESTING)
  add_executable(test_spectral tests/test_spectral.cpp)
  target_link_libraries(test_spectral PRIVATE al_maleh_spectral)
  add_test(NAME generalized_spectrum COMMAND test_spectral)
  add_executable(test_spectral_export tests/test_spectral_export.cpp)
  target_link_libraries(test_spectral_export PRIVATE al_maleh_spectral)
  add_test(NAME spectral_export_roundtrip COMMAND test_spectral_export)
endif()

# Run the complete spectrum CLI and verify that it exported diagnostics.
if(BUILD_TESTING)
  add_test(NAME m3_spectrum_cli COMMAND cartesian_spectrum 5 5 1 3.3 "${CMAKE_CURRENT_BINARY_DIR}/m3_ci")
  add_test(NAME m3_spectrum_artifacts COMMAND ${CMAKE_COMMAND}
           -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/m3_ci
           -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_m3_csv.cmake)
  set_tests_properties(m3_spectrum_artifacts PROPERTIES DEPENDS m3_spectrum_cli)
endif()

# M4: analytic TM curl-shape matching and controlled validation campaign.
add_library(al_maleh_tm_validation STATIC src/tm_validation.cpp)
target_include_directories(al_maleh_tm_validation PUBLIC include)
target_link_libraries(al_maleh_tm_validation PUBLIC al_maleh_spectral Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_tm_validation PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(tm_validation apps/tm_validation.cpp)
target_link_libraries(tm_validation PRIVATE al_maleh_tm_validation)
if(BUILD_TESTING)
  add_executable(test_tm_validation tests/test_tm_validation.cpp)
  target_link_libraries(test_tm_validation PRIVATE al_maleh_tm_validation)
  add_test(NAME tm_overlap_validation COMMAND test_tm_validation)
  add_test(NAME m4_smoke COMMAND tm_validation 5 5 1 3.3 3
    "${CMAKE_CURRENT_BINARY_DIR}/m4_smoke")
endif()

# Scientific post-processing is optional for core C++ builds. CI installs it.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(BUILD_TESTING AND Python3_Interpreter_FOUND)
  # Independent stdlib-only audit: verify that low positive EFG modes
  # and the failing doublet remain visible in the archived campaign.
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/docs/results/independent_validation/LATEST")
    add_test(NAME archived_spectral_pollution
      COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -p test_spectral_pollution.py -v)
    set_tests_properties(archived_spectral_pollution PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_test(NAME validation_manifest_integrity
    COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -p test_validation_manifest.py -v)
  set_tests_properties(validation_manifest_integrity PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  execute_process(COMMAND ${CMAKE_COMMAND} -E env PYTHONNOUSERSITE=1
    ${Python3_EXECUTABLE} -c "import numpy, scipy, matplotlib"
    RESULT_VARIABLE AL_MALEH_SCIENTIFIC_PYTHON OUTPUT_QUIET ERROR_QUIET)
  if(AL_MALEH_SCIENTIFIC_PYTHON EQUAL 0)
    add_test(NAME independent_validation_algorithms
      COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests -p test_independent_validation.py -v)
    set_tests_properties(independent_validation_algorithms PROPERTIES
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      ENVIRONMENT "PYTHONNOUSERSITE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/docs/results/independent_validation/LATEST")
      file(READ "${CMAKE_CURRENT_SOURCE_DIR}/docs/results/independent_validation/LATEST" AL_MALEH_CAMPAIGN)
      string(STRIP "${AL_MALEH_CAMPAIGN}" AL_MALEH_CAMPAIGN)
      add_test(NAME independent_published_results
        COMMAND ${Python3_EXECUTABLE} scripts/independent_validation.py --check-only
          --campaign "docs/results/independent_validation/${AL_MALEH_CAMPAIGN}")
      add_test(NAME independent_published_hashes
        COMMAND ${Python3_EXECUTABLE} scripts/validation_manifest.py --verify
          --campaign "docs/results/independent_validation/${AL_MALEH_CAMPAIGN}")
      set_tests_properties(independent_published_results independent_published_hashes PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        ENVIRONMENT "PYTHONNOUSERSITE=1;OPENBLAS_NUM_THREADS=1;OMP_NUM_THREADS=1")
    endif()
  endif()
endif()

# M6b: independent divergence/normal-boundary forms (no penalty to K/M).
add_library(al_maleh_constraints STATIC src/constraints.cpp)
target_include_directories(al_maleh_constraints PUBLIC include)
target_link_libraries(al_maleh_constraints PUBLIC al_maleh_spectral Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_constraints PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(constraint_subspace apps/constraint_subspace.cpp)
target_link_libraries(constraint_subspace PRIVATE al_maleh_constraints)
if(BUILD_TESTING)
  add_executable(test_constraints tests/test_constraints.cpp)
  target_link_libraries(test_constraints PRIVATE al_maleh_constraints)
  add_test(NAME m6b_constraint_forms COMMAND test_constraints)
  add_test(NAME m6b_subspace_smoke COMMAND constraint_subspace 9 7 70 71
      "${CMAKE_CURRENT_BINARY_DIR}/m6b_smoke")
endif()
