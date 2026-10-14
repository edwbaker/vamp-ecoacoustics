cmake_minimum_required(VERSION 3.10)
project(vamp-ecoacoustics)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to Release build for smaller binary size
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Strip symbols in Release mode for smaller binary
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(MINGW OR CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")
    endif()
endif()

# Source files
file(GLOB PLUGIN_SOURCES src/*.cpp)
file(GLOB PLUGIN_HEADERS src/*.h)

# Test runner: the plugin classes compiled straight into an executable (no Vamp host
# needed), so leave out the files that only exist to export the shared library
set(TEST_SOURCES ${PLUGIN_SOURCES})
list(FILTER TEST_SOURCES EXCLUDE REGEX "(plugins|PluginAdapter)\\.cpp$")

# Add module definition file for Windows exports
if(WIN32)
    set(PLUGIN_SOURCES ${PLUGIN_SOURCES} src/vamp-ecoacoustics.def)
endif()

add_library(vamp-ecoacoustics SHARED ${PLUGIN_SOURCES} ${PLUGIN_HEADERS})

target_include_directories(vamp-ecoacoustics PRIVATE src ext ext/vamp ext/vamp-sdk)

# Link math library on Unix (not needed on Windows)
if(NOT WIN32)
    target_link_libraries(vamp-ecoacoustics PRIVATE m)
endif()

# Set output name for the plugin
set_target_properties(vamp-ecoacoustics PROPERTIES 
    OUTPUT_NAME "vamp-ecoacoustics"
    PREFIX ""
)

# Tests: cmake --build build --target test_indices && ctest --test-dir build -C Release
enable_testing()
add_executable(test_indices EXCLUDE_FROM_ALL tests/test_indices.cpp ${TEST_SOURCES})
target_include_directories(test_indices PRIVATE src ext ext/vamp ext/vamp-sdk tests)
if(NOT WIN32)
    target_link_libraries(test_indices PRIVATE m)
endif()
add_test(NAME indices COMMAND test_indices)

# Optionally, install target
install(TARGETS vamp-ecoacoustics DESTINATION lib)
