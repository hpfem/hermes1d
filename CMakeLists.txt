project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(WITH_PYTHON no)

# allow to override the default values in CMake.vars
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)


set(HERMES_BIN hermes1d)

add_subdirectory(src)
if(WITH_PYTHON)
    add_subdirectory(python)
endif(WITH_PYTHON)
add_subdirectory(examples)
