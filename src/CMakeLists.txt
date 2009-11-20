set(SRC
    common.cpp iterator.cpp qsort.cpp
    lobatto.cpp legendre.cpp  matrix.cpp discrete.cpp mesh.cpp linearizer.cpp
    quad_std.cpp transforms.cpp adapt.cpp graph.cpp h1_polys.cpp
    )

add_library(${HERMES_BIN} SHARED ${SRC})
if(DEBUG)
    set_target_properties(hermes1d-debug PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_DEBUG})
endif(DEBUG)
if(RELEASE)
    set_target_properties(hermes1d PROPERTIES COMPILE_FLAGS
        ${CMAKE_CXX_FLAGS_RELEASE})
endif(RELEASE)
