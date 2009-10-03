set(SRC
    common.cpp
    lobatto.cpp  matrix.cpp discrete.cpp
    )

add_library(${HERMES_BIN} SHARED ${SRC})
