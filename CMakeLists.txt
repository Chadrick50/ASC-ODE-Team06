cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)

include(FetchContent)
FetchContent_Declare(Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.4)
FetchContent_MakeAvailable(Catch2)


include_directories(src nanoblas/src)

add_subdirectory (nanoblas)


add_executable (test_mass_spring demos/test_mass_spring.cpp)
target_link_libraries (test_mass_spring PUBLIC nanoblas)

add_executable (test_RC_circuit demos/test_RC_circuit.cpp)
target_link_libraries (test_RC_circuit PUBLIC nanoblas)

add_executable (demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries (demo_autodiff PUBLIC nanoblas)

enable_testing()
add_executable(tests tests/testConvergence.cpp)
target_link_libraries(tests PRIVATE nanoblas Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(tests)
