include_directories(${hermes_common_SOURCE_DIR})

project(solvers)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES} ${HERMES_COMMON})
set_source_files_properties(main.cpp PROPERTIES
    COMPILE_FLAGS "-std=c++0x"
    )



# tests:
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(solvers ${BIN})
