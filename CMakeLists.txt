project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(WITH_PYTHON no)
set(DEBUG no)
set(RELEASE yes)
set(WITH_TESTS yes)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)


if(DEBUG)
        set(HERMES_BIN hermes1d-debug)
endif(DEBUG)

if(RELEASE)
        set(HERMES_BIN hermes1d)
endif(RELEASE)


add_subdirectory(src)
if(WITH_PYTHON)
    add_subdirectory(python)
endif(WITH_PYTHON)
add_subdirectory(examples)

if(WITH_TESTS)
    enable_testing()
	add_subdirectory(tests)
endif(WITH_TESTS)

message("\nConfiguration results:")
message("----------------------")
message("Build with Python: ${WITH_PYTHON}")
message("Build with debug: ${DEBUG}")
message("Build with release: ${RELEASE}")
message("Build with tests: ${WITH_TESTS}")
message("\n")
