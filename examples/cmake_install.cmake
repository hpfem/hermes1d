# Install script for directory: /home/petr/Repositories/hermes-1d/examples

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/petr/Repositories/hermes-1d/examples/first_order_general/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_dirichlet/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_neumann/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_newton/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_newton2/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_exp/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_sin/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_sin_eigen/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_sin_adapt/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_pendulum/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_chaotic/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_car_model/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_dirichlet_adapt/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/first_order_general_adapt/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_dirichlet_adapt_new/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/laplace_bc_dirichlet_adapt_goal/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/neutronics/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_neutronics_eigenvalue/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_neutronics_fixedsrc/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/system_neutronics_fixedsrc2/cmake_install.cmake")
  include("/home/petr/Repositories/hermes-1d/examples/schroedinger/cmake_install.cmake")

endif()

