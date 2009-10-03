set(SRC
    common.cpp
    lobatto.cpp  matrix.cpp
    )

add_library(${HERMES_BIN} SHARED ${SRC})
