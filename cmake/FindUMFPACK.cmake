#
# UMFPACK
#

FIND_PATH(UMFPACK_INCLUDE_DIR umfpack.h ${UMFPACK_ROOT}/Include /usr/include /usr/include/umfpack /usr/local/include/UMFPACK /usr/include/suitesparse /opt/local/include/ufsparse)
FIND_PATH(AMD_INCLUDE_DIR amd.h ${AMD_ROOT}/Include /usr/include /usr/local/include/AMD /usr/include/suitesparse /opt/local/include/ufsparse)

FIND_LIBRARY(UMFPACK_LIBRARY umfpack ${UMFPACK_ROOT}/Lib /usr/lib /usr/local/lib/UMFPACK /opt/local/lib)
FIND_LIBRARY(AMD_LIBRARY amd ${AMD_ROOT}/Lib /usr/lib /usr/local/lib/AMD /opt/local/lib)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(UMFPACK DEFAULT_MSG UMFPACK_INCLUDE_DIR AMD_INCLUDE_DIR UMFPACK_LIBRARY AMD_LIBRARY)
