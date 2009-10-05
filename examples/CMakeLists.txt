add_subdirectory(laplace)

enable_language(Fortran)
find_package(UMFPACK REQUIRED)
find_package(BLAS REQUIRED)

if(WITH_PYTHON)
    add_subdirectory(schroedinger)
endif(WITH_PYTHON)
