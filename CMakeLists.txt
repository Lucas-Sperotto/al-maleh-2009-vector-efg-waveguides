cmake_minimum_required(VERSION 3.16)
project(al_maleh_vector_efg LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
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
add_library(al_maleh_spectral STATIC src/spectral.cpp)
target_include_directories(al_maleh_spectral PUBLIC include)
target_link_libraries(al_maleh_spectral PUBLIC al_maleh_cartesian Eigen3::Eigen)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(al_maleh_spectral PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_executable(cartesian_spectrum apps/cartesian_spectrum.cpp)
target_link_libraries(cartesian_spectrum PRIVATE al_maleh_spectral)
if(BUILD_TESTING)
  add_executable(test_spectral tests/test_spectral.cpp)
  target_link_libraries(test_spectral PRIVATE al_maleh_spectral)
  add_test(NAME generalized_spectrum COMMAND test_spectral)
endif()

# Run the complete spectrum CLI and verify that it exported diagnostics.
if(BUILD_TESTING)
  add_test(NAME m3_spectrum_cli COMMAND cartesian_spectrum 5 5 1 3.3 "${CMAKE_CURRENT_BINARY_DIR}/m3_ci")
  add_test(NAME m3_spectrum_artifacts COMMAND ${CMAKE_COMMAND}
           -DOUTPUT_DIR=${CMAKE_CURRENT_BINARY_DIR}/m3_ci
           -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_m3_csv.cmake)
  set_tests_properties(m3_spectrum_artifacts PROPERTIES DEPENDS m3_spectrum_cli)
endif()
