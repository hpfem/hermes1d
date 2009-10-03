set(SRC
    common.cpp
    lobatto.cpp  matrix.cpp weakform.cpp
    )

add_library(${HERMES_BIN} SHARED ${SRC})
