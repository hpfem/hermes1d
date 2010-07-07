SET(CYTHON_INCLUDE_DIRECTORIES . -I../../hermes_common/)
cython_add_module(_hermes1d)
target_link_libraries(_hermes1d hermes_common ${HERMES_BIN})
set_property(TARGET _hermes1d PROPERTY COMPILE_FLAGS "-g")

install(TARGETS _hermes1d
            RUNTIME DESTINATION lib/python/hermes1d
            ARCHIVE DESTINATION lib/python/hermes1d
            LIBRARY DESTINATION lib/python/hermes1d
        )
set_target_properties(_hermes1d PROPERTIES
    INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib/python/hermes1d"
    )

install(FILES __init__.py
    DESTINATION lib/python/hermes1d
    )

install(FILES config.h stdcython.h _hermes2d_api.h
    DESTINATION lib/python/hermes1d/include
    )

install(FILES _hermes1d.pxd __init__.py
    DESTINATION lib/python/hermes1d/include/hermes2d
    )
