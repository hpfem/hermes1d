project(hermes1d)

cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(HERMES_BIN hermes1d)


add_subdirectory(src)
add_subdirectory(python)
add_subdirectory(examples)
