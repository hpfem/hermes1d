set(SRC
    common.cpp
    lobatto.cpp  matrix.cpp discrete.cpp mesh.cpp quad_std.cpp
    )

add_library(${HERMES_BIN} SHARED ${SRC})
