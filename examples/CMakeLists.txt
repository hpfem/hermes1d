enable_language(Fortran)
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

add_subdirectory(laplace)
add_subdirectory(laplace_bc_neumann)
add_subdirectory(laplace_bc_newton)
add_subdirectory(laplace_bc_newton2)
add_subdirectory(system)

if(WITH_PYTHON)
    add_subdirectory(schroedinger)
endif(WITH_PYTHON)
