cmake_minimum_required(VERSION 3.16)
include(FetchContent)
FetchContent_Declare(
        doctest
        GIT_REPOSITORY https://github.com/onqtam/doctest.git
        GIT_TAG v2.4.12
)
FetchContent_MakeAvailable(doctest)
enable_testing()

###########
# test_io #
###########

add_executable(test_io test_io.cc)
target_compile_features(test_io PRIVATE cxx_std_20)
target_include_directories(test_io PRIVATE ${doctest_SOURCE_DIR})
target_link_libraries(test_io PRIVATE ${PROJECT_NAME})
set_target_properties(test_io PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_test(NAME test_io COMMAND test_io)
set_property(TEST test_io PROPERTY ENVIRONMENT CUES_TESTS_ROOT=${CMAKE_SOURCE_DIR})

################
# test_signals #
################

add_executable(test_signals test_signals.cc)
target_compile_features(test_signals PRIVATE cxx_std_20)
target_include_directories(test_signals PRIVATE ${doctest_SOURCE_DIR})
target_link_libraries(test_signals PRIVATE ${PROJECT_NAME})
set_target_properties(test_signals PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
add_test(NAME test_signals COMMAND test_signals)
set_property(TEST test_signals PROPERTY ENVIRONMENT CUES_TESTS_ROOT=${CMAKE_SOURCE_DIR})
