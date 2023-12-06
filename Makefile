# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/petr/Repositories/hermes-1d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/petr/Repositories/hermes-1d

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "No interactive CMake dialog available..."
	/usr/bin/cmake -E echo No\ interactive\ CMake\ dialog\ available.
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local/fast

# Special rule for the target install/strip
install/strip: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip

# Special rule for the target install/strip
install/strip/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing the project stripped..."
	/usr/bin/cmake -DCMAKE_INSTALL_DO_STRIP=1 -P cmake_install.cmake
.PHONY : install/strip/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/petr/Repositories/hermes-1d/CMakeFiles /home/petr/Repositories/hermes-1d//CMakeFiles/progress.marks
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/petr/Repositories/hermes-1d/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named doc

# Build rule for target.
doc: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 doc
.PHONY : doc

# fast build rule for target.
doc/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/doc.dir/build.make CMakeFiles/doc.dir/build
.PHONY : doc/fast

#=============================================================================
# Target rules for targets named doc-tex

# Build rule for target.
doc-tex: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 doc-tex
.PHONY : doc-tex

# fast build rule for target.
doc-tex/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/doc-tex.dir/build.make CMakeFiles/doc-tex.dir/build
.PHONY : doc-tex/fast

#=============================================================================
# Target rules for targets named test-quick

# Build rule for target.
test-quick: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 test-quick
.PHONY : test-quick

# fast build rule for target.
test-quick/fast:
	$(MAKE) $(MAKESILENT) -f CMakeFiles/test-quick.dir/build.make CMakeFiles/test-quick.dir/build
.PHONY : test-quick/fast

#=============================================================================
# Target rules for targets named hermes_common

# Build rule for target.
hermes_common: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 hermes_common
.PHONY : hermes_common

# fast build rule for target.
hermes_common/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/CMakeFiles/hermes_common.dir/build.make hermes_common/CMakeFiles/hermes_common.dir/build
.PHONY : hermes_common/fast

#=============================================================================
# Target rules for targets named _hermes_common

# Build rule for target.
_hermes_common: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 _hermes_common
.PHONY : _hermes_common

# fast build rule for target.
_hermes_common/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/CMakeFiles/_hermes_common.dir/build.make hermes_common/CMakeFiles/_hermes_common.dir/build
.PHONY : _hermes_common/fast

#=============================================================================
# Target rules for targets named sparse

# Build rule for target.
sparse: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 sparse
.PHONY : sparse

# fast build rule for target.
sparse/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/sparselib/CMakeFiles/sparse.dir/build.make hermes_common/sparselib/CMakeFiles/sparse.dir/build
.PHONY : sparse/fast

#=============================================================================
# Target rules for targets named spblas

# Build rule for target.
spblas: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 spblas
.PHONY : spblas

# fast build rule for target.
spblas/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/sparselib/spblas/CMakeFiles/spblas.dir/build.make hermes_common/sparselib/spblas/CMakeFiles/spblas.dir/build
.PHONY : spblas/fast

#=============================================================================
# Target rules for targets named mv

# Build rule for target.
mv: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 mv
.PHONY : mv

# fast build rule for target.
mv/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/sparselib/mv/CMakeFiles/mv.dir/build.make hermes_common/sparselib/mv/CMakeFiles/mv.dir/build
.PHONY : mv/fast

#=============================================================================
# Target rules for targets named python-in-cpp

# Build rule for target.
python-in-cpp: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 python-in-cpp
.PHONY : python-in-cpp

# fast build rule for target.
python-in-cpp/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/python-in-cpp/CMakeFiles/python-in-cpp.dir/build.make hermes_common/tests/python-in-cpp/CMakeFiles/python-in-cpp.dir/build
.PHONY : python-in-cpp/fast

#=============================================================================
# Target rules for targets named numpy-in-cpp

# Build rule for target.
numpy-in-cpp: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 numpy-in-cpp
.PHONY : numpy-in-cpp

# fast build rule for target.
numpy-in-cpp/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/numpy-in-cpp/CMakeFiles/numpy-in-cpp.dir/build.make hermes_common/tests/numpy-in-cpp/CMakeFiles/numpy-in-cpp.dir/build
.PHONY : numpy-in-cpp/fast

#=============================================================================
# Target rules for targets named matrix

# Build rule for target.
matrix: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 matrix
.PHONY : matrix

# fast build rule for target.
matrix/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/matrix/CMakeFiles/matrix.dir/build.make hermes_common/tests/matrix/CMakeFiles/matrix.dir/build
.PHONY : matrix/fast

#=============================================================================
# Target rules for targets named vector

# Build rule for target.
vector: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 vector
.PHONY : vector

# fast build rule for target.
vector/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/vector/CMakeFiles/vector.dir/build.make hermes_common/tests/vector/CMakeFiles/vector.dir/build
.PHONY : vector/fast

#=============================================================================
# Target rules for targets named matrix-io

# Build rule for target.
matrix-io: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 matrix-io
.PHONY : matrix-io

# fast build rule for target.
matrix-io/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/matrix-io/CMakeFiles/matrix-io.dir/build.make hermes_common/tests/matrix-io/CMakeFiles/matrix-io.dir/build
.PHONY : matrix-io/fast

#=============================================================================
# Target rules for targets named solvers

# Build rule for target.
solvers: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 solvers
.PHONY : solvers

# fast build rule for target.
solvers/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/solvers/CMakeFiles/solvers.dir/build.make hermes_common/tests/solvers/CMakeFiles/solvers.dir/build
.PHONY : solvers/fast

#=============================================================================
# Target rules for targets named leaks

# Build rule for target.
leaks: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 leaks
.PHONY : leaks

# fast build rule for target.
leaks/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/leaks/CMakeFiles/leaks.dir/build.make hermes_common/tests/leaks/CMakeFiles/leaks.dir/build
.PHONY : leaks/fast

#=============================================================================
# Target rules for targets named my_api2

# Build rule for target.
my_api2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 my_api2
.PHONY : my_api2

# fast build rule for target.
my_api2/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/leaks/CMakeFiles/my_api2.dir/build.make hermes_common/tests/leaks/CMakeFiles/my_api2.dir/build
.PHONY : my_api2/fast

#=============================================================================
# Target rules for targets named cpp-callbacks

# Build rule for target.
cpp-callbacks: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 cpp-callbacks
.PHONY : cpp-callbacks

# fast build rule for target.
cpp-callbacks/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/cpp-callbacks/CMakeFiles/cpp-callbacks.dir/build.make hermes_common/tests/cpp-callbacks/CMakeFiles/cpp-callbacks.dir/build
.PHONY : cpp-callbacks/fast

#=============================================================================
# Target rules for targets named my_api

# Build rule for target.
my_api: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 my_api
.PHONY : my_api

# fast build rule for target.
my_api/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/cpp-callbacks/CMakeFiles/my_api.dir/build.make hermes_common/tests/cpp-callbacks/CMakeFiles/my_api.dir/build
.PHONY : my_api/fast

#=============================================================================
# Target rules for targets named timer-1

# Build rule for target.
timer-1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 timer-1
.PHONY : timer-1

# fast build rule for target.
timer-1/fast:
	$(MAKE) $(MAKESILENT) -f hermes_common/tests/timer/CMakeFiles/timer-1.dir/build.make hermes_common/tests/timer/CMakeFiles/timer-1.dir/build
.PHONY : timer-1/fast

#=============================================================================
# Target rules for targets named hermes1d-debug

# Build rule for target.
hermes1d-debug: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 hermes1d-debug
.PHONY : hermes1d-debug

# fast build rule for target.
hermes1d-debug/fast:
	$(MAKE) $(MAKESILENT) -f src/CMakeFiles/hermes1d-debug.dir/build.make src/CMakeFiles/hermes1d-debug.dir/build
.PHONY : hermes1d-debug/fast

#=============================================================================
# Target rules for targets named h1d_wrapper

# Build rule for target.
h1d_wrapper: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 h1d_wrapper
.PHONY : h1d_wrapper

# fast build rule for target.
h1d_wrapper/fast:
	$(MAKE) $(MAKESILENT) -f hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/build.make hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/build
.PHONY : h1d_wrapper/fast

#=============================================================================
# Target rules for targets named first_order_general

# Build rule for target.
first_order_general: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 first_order_general
.PHONY : first_order_general

# fast build rule for target.
first_order_general/fast:
	$(MAKE) $(MAKESILENT) -f examples/first_order_general/CMakeFiles/first_order_general.dir/build.make examples/first_order_general/CMakeFiles/first_order_general.dir/build
.PHONY : first_order_general/fast

#=============================================================================
# Target rules for targets named laplace_bc_dirichlet

# Build rule for target.
laplace_bc_dirichlet: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_dirichlet
.PHONY : laplace_bc_dirichlet

# fast build rule for target.
laplace_bc_dirichlet/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_dirichlet/CMakeFiles/laplace_bc_dirichlet.dir/build.make examples/laplace_bc_dirichlet/CMakeFiles/laplace_bc_dirichlet.dir/build
.PHONY : laplace_bc_dirichlet/fast

#=============================================================================
# Target rules for targets named laplace_bc_neumann

# Build rule for target.
laplace_bc_neumann: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_neumann
.PHONY : laplace_bc_neumann

# fast build rule for target.
laplace_bc_neumann/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_neumann/CMakeFiles/laplace_bc_neumann.dir/build.make examples/laplace_bc_neumann/CMakeFiles/laplace_bc_neumann.dir/build
.PHONY : laplace_bc_neumann/fast

#=============================================================================
# Target rules for targets named laplace_bc_newton

# Build rule for target.
laplace_bc_newton: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_newton
.PHONY : laplace_bc_newton

# fast build rule for target.
laplace_bc_newton/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_newton/CMakeFiles/laplace_bc_newton.dir/build.make examples/laplace_bc_newton/CMakeFiles/laplace_bc_newton.dir/build
.PHONY : laplace_bc_newton/fast

#=============================================================================
# Target rules for targets named laplace_bc_newton2

# Build rule for target.
laplace_bc_newton2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_newton2
.PHONY : laplace_bc_newton2

# fast build rule for target.
laplace_bc_newton2/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_newton2/CMakeFiles/laplace_bc_newton2.dir/build.make examples/laplace_bc_newton2/CMakeFiles/laplace_bc_newton2.dir/build
.PHONY : laplace_bc_newton2/fast

#=============================================================================
# Target rules for targets named system_exp

# Build rule for target.
system_exp: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_exp
.PHONY : system_exp

# fast build rule for target.
system_exp/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_exp/CMakeFiles/system_exp.dir/build.make examples/system_exp/CMakeFiles/system_exp.dir/build
.PHONY : system_exp/fast

#=============================================================================
# Target rules for targets named system_sin

# Build rule for target.
system_sin: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_sin
.PHONY : system_sin

# fast build rule for target.
system_sin/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_sin/CMakeFiles/system_sin.dir/build.make examples/system_sin/CMakeFiles/system_sin.dir/build
.PHONY : system_sin/fast

#=============================================================================
# Target rules for targets named system_sin_eigen

# Build rule for target.
system_sin_eigen: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_sin_eigen
.PHONY : system_sin_eigen

# fast build rule for target.
system_sin_eigen/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_sin_eigen/CMakeFiles/system_sin_eigen.dir/build.make examples/system_sin_eigen/CMakeFiles/system_sin_eigen.dir/build
.PHONY : system_sin_eigen/fast

#=============================================================================
# Target rules for targets named system_sin_adapt

# Build rule for target.
system_sin_adapt: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_sin_adapt
.PHONY : system_sin_adapt

# fast build rule for target.
system_sin_adapt/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_sin_adapt/CMakeFiles/system_sin_adapt.dir/build.make examples/system_sin_adapt/CMakeFiles/system_sin_adapt.dir/build
.PHONY : system_sin_adapt/fast

#=============================================================================
# Target rules for targets named system_pendulum

# Build rule for target.
system_pendulum: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_pendulum
.PHONY : system_pendulum

# fast build rule for target.
system_pendulum/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_pendulum/CMakeFiles/system_pendulum.dir/build.make examples/system_pendulum/CMakeFiles/system_pendulum.dir/build
.PHONY : system_pendulum/fast

#=============================================================================
# Target rules for targets named system_chaotic

# Build rule for target.
system_chaotic: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_chaotic
.PHONY : system_chaotic

# fast build rule for target.
system_chaotic/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_chaotic/CMakeFiles/system_chaotic.dir/build.make examples/system_chaotic/CMakeFiles/system_chaotic.dir/build
.PHONY : system_chaotic/fast

#=============================================================================
# Target rules for targets named system_car_model

# Build rule for target.
system_car_model: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_car_model
.PHONY : system_car_model

# fast build rule for target.
system_car_model/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_car_model/CMakeFiles/system_car_model.dir/build.make examples/system_car_model/CMakeFiles/system_car_model.dir/build
.PHONY : system_car_model/fast

#=============================================================================
# Target rules for targets named laplace_bc_dirichlet_adapt

# Build rule for target.
laplace_bc_dirichlet_adapt: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_dirichlet_adapt
.PHONY : laplace_bc_dirichlet_adapt

# fast build rule for target.
laplace_bc_dirichlet_adapt/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_dirichlet_adapt/CMakeFiles/laplace_bc_dirichlet_adapt.dir/build.make examples/laplace_bc_dirichlet_adapt/CMakeFiles/laplace_bc_dirichlet_adapt.dir/build
.PHONY : laplace_bc_dirichlet_adapt/fast

#=============================================================================
# Target rules for targets named first_order_general_adapt

# Build rule for target.
first_order_general_adapt: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 first_order_general_adapt
.PHONY : first_order_general_adapt

# fast build rule for target.
first_order_general_adapt/fast:
	$(MAKE) $(MAKESILENT) -f examples/first_order_general_adapt/CMakeFiles/first_order_general_adapt.dir/build.make examples/first_order_general_adapt/CMakeFiles/first_order_general_adapt.dir/build
.PHONY : first_order_general_adapt/fast

#=============================================================================
# Target rules for targets named laplace_bc_dirichlet_adapt_new

# Build rule for target.
laplace_bc_dirichlet_adapt_new: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_dirichlet_adapt_new
.PHONY : laplace_bc_dirichlet_adapt_new

# fast build rule for target.
laplace_bc_dirichlet_adapt_new/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_dirichlet_adapt_new/CMakeFiles/laplace_bc_dirichlet_adapt_new.dir/build.make examples/laplace_bc_dirichlet_adapt_new/CMakeFiles/laplace_bc_dirichlet_adapt_new.dir/build
.PHONY : laplace_bc_dirichlet_adapt_new/fast

#=============================================================================
# Target rules for targets named first_order_general_adapt_new

# Build rule for target.
first_order_general_adapt_new: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 first_order_general_adapt_new
.PHONY : first_order_general_adapt_new

# fast build rule for target.
first_order_general_adapt_new/fast:
	$(MAKE) $(MAKESILENT) -f examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/build.make examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/build
.PHONY : first_order_general_adapt_new/fast

#=============================================================================
# Target rules for targets named laplace_bc_dirichlet_adapt_goal

# Build rule for target.
laplace_bc_dirichlet_adapt_goal: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 laplace_bc_dirichlet_adapt_goal
.PHONY : laplace_bc_dirichlet_adapt_goal

# fast build rule for target.
laplace_bc_dirichlet_adapt_goal/fast:
	$(MAKE) $(MAKESILENT) -f examples/laplace_bc_dirichlet_adapt_goal/CMakeFiles/laplace_bc_dirichlet_adapt_goal.dir/build.make examples/laplace_bc_dirichlet_adapt_goal/CMakeFiles/laplace_bc_dirichlet_adapt_goal.dir/build
.PHONY : laplace_bc_dirichlet_adapt_goal/fast

#=============================================================================
# Target rules for targets named neutronics

# Build rule for target.
neutronics: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 neutronics
.PHONY : neutronics

# fast build rule for target.
neutronics/fast:
	$(MAKE) $(MAKESILENT) -f examples/neutronics/CMakeFiles/neutronics.dir/build.make examples/neutronics/CMakeFiles/neutronics.dir/build
.PHONY : neutronics/fast

#=============================================================================
# Target rules for targets named system_neutronics_eigenvalue

# Build rule for target.
system_neutronics_eigenvalue: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_neutronics_eigenvalue
.PHONY : system_neutronics_eigenvalue

# fast build rule for target.
system_neutronics_eigenvalue/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_neutronics_eigenvalue/CMakeFiles/system_neutronics_eigenvalue.dir/build.make examples/system_neutronics_eigenvalue/CMakeFiles/system_neutronics_eigenvalue.dir/build
.PHONY : system_neutronics_eigenvalue/fast

#=============================================================================
# Target rules for targets named system_neutronics_fixedsrc

# Build rule for target.
system_neutronics_fixedsrc: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_neutronics_fixedsrc
.PHONY : system_neutronics_fixedsrc

# fast build rule for target.
system_neutronics_fixedsrc/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_neutronics_fixedsrc/CMakeFiles/system_neutronics_fixedsrc.dir/build.make examples/system_neutronics_fixedsrc/CMakeFiles/system_neutronics_fixedsrc.dir/build
.PHONY : system_neutronics_fixedsrc/fast

#=============================================================================
# Target rules for targets named system_neutronics_fixedsrc2

# Build rule for target.
system_neutronics_fixedsrc2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 system_neutronics_fixedsrc2
.PHONY : system_neutronics_fixedsrc2

# fast build rule for target.
system_neutronics_fixedsrc2/fast:
	$(MAKE) $(MAKESILENT) -f examples/system_neutronics_fixedsrc2/CMakeFiles/system_neutronics_fixedsrc2.dir/build.make examples/system_neutronics_fixedsrc2/CMakeFiles/system_neutronics_fixedsrc2.dir/build
.PHONY : system_neutronics_fixedsrc2/fast

#=============================================================================
# Target rules for targets named schroedinger

# Build rule for target.
schroedinger: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 schroedinger
.PHONY : schroedinger

# fast build rule for target.
schroedinger/fast:
	$(MAKE) $(MAKESILENT) -f examples/schroedinger/CMakeFiles/schroedinger.dir/build.make examples/schroedinger/CMakeFiles/schroedinger.dir/build
.PHONY : schroedinger/fast

#=============================================================================
# Target rules for targets named legendre-1

# Build rule for target.
legendre-1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 legendre-1
.PHONY : legendre-1

# fast build rule for target.
legendre-1/fast:
	$(MAKE) $(MAKESILENT) -f tests/legendre-1/CMakeFiles/legendre-1.dir/build.make tests/legendre-1/CMakeFiles/legendre-1.dir/build
.PHONY : legendre-1/fast

#=============================================================================
# Target rules for targets named legendre-2

# Build rule for target.
legendre-2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 legendre-2
.PHONY : legendre-2

# fast build rule for target.
legendre-2/fast:
	$(MAKE) $(MAKESILENT) -f tests/legendre-2/CMakeFiles/legendre-2.dir/build.make tests/legendre-2/CMakeFiles/legendre-2.dir/build
.PHONY : legendre-2/fast

#=============================================================================
# Target rules for targets named legendre-3

# Build rule for target.
legendre-3: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 legendre-3
.PHONY : legendre-3

# fast build rule for target.
legendre-3/fast:
	$(MAKE) $(MAKESILENT) -f tests/legendre-3/CMakeFiles/legendre-3.dir/build.make tests/legendre-3/CMakeFiles/legendre-3.dir/build
.PHONY : legendre-3/fast

#=============================================================================
# Target rules for targets named lobatto-1

# Build rule for target.
lobatto-1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 lobatto-1
.PHONY : lobatto-1

# fast build rule for target.
lobatto-1/fast:
	$(MAKE) $(MAKESILENT) -f tests/lobatto-1/CMakeFiles/lobatto-1.dir/build.make tests/lobatto-1/CMakeFiles/lobatto-1.dir/build
.PHONY : lobatto-1/fast

#=============================================================================
# Target rules for targets named lobatto-2

# Build rule for target.
lobatto-2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 lobatto-2
.PHONY : lobatto-2

# fast build rule for target.
lobatto-2/fast:
	$(MAKE) $(MAKESILENT) -f tests/lobatto-2/CMakeFiles/lobatto-2.dir/build.make tests/lobatto-2/CMakeFiles/lobatto-2.dir/build
.PHONY : lobatto-2/fast

#=============================================================================
# Target rules for targets named adapt-exact-quadr-L2

# Build rule for target.
adapt-exact-quadr-L2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-quadr-L2
.PHONY : adapt-exact-quadr-L2

# fast build rule for target.
adapt-exact-quadr-L2/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-quadr-L2/CMakeFiles/adapt-exact-quadr-L2.dir/build.make tests/adapt-exact-quadr-L2/CMakeFiles/adapt-exact-quadr-L2.dir/build
.PHONY : adapt-exact-quadr-L2/fast

#=============================================================================
# Target rules for targets named adapt-exact-sin-L2

# Build rule for target.
adapt-exact-sin-L2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-sin-L2
.PHONY : adapt-exact-sin-L2

# fast build rule for target.
adapt-exact-sin-L2/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-sin-L2/CMakeFiles/adapt-exact-sin-L2.dir/build.make tests/adapt-exact-sin-L2/CMakeFiles/adapt-exact-sin-L2.dir/build
.PHONY : adapt-exact-sin-L2/fast

#=============================================================================
# Target rules for targets named adapt-exact-system-sin-L2

# Build rule for target.
adapt-exact-system-sin-L2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-system-sin-L2
.PHONY : adapt-exact-system-sin-L2

# fast build rule for target.
adapt-exact-system-sin-L2/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-system-sin-L2/CMakeFiles/adapt-exact-system-sin-L2.dir/build.make tests/adapt-exact-system-sin-L2/CMakeFiles/adapt-exact-system-sin-L2.dir/build
.PHONY : adapt-exact-system-sin-L2/fast

#=============================================================================
# Target rules for targets named adapt-exact-quadr-H1-solver0

# Build rule for target.
adapt-exact-quadr-H1-solver0: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-quadr-H1-solver0
.PHONY : adapt-exact-quadr-H1-solver0

# fast build rule for target.
adapt-exact-quadr-H1-solver0/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-quadr-H1-solver0/CMakeFiles/adapt-exact-quadr-H1-solver0.dir/build.make tests/adapt-exact-quadr-H1-solver0/CMakeFiles/adapt-exact-quadr-H1-solver0.dir/build
.PHONY : adapt-exact-quadr-H1-solver0/fast

#=============================================================================
# Target rules for targets named adapt-exact-quadr-H1-solver1

# Build rule for target.
adapt-exact-quadr-H1-solver1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-quadr-H1-solver1
.PHONY : adapt-exact-quadr-H1-solver1

# fast build rule for target.
adapt-exact-quadr-H1-solver1/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-quadr-H1-solver1/CMakeFiles/adapt-exact-quadr-H1-solver1.dir/build.make tests/adapt-exact-quadr-H1-solver1/CMakeFiles/adapt-exact-quadr-H1-solver1.dir/build
.PHONY : adapt-exact-quadr-H1-solver1/fast

#=============================================================================
# Target rules for targets named adapt-exact-quadr-H1-solver2

# Build rule for target.
adapt-exact-quadr-H1-solver2: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-quadr-H1-solver2
.PHONY : adapt-exact-quadr-H1-solver2

# fast build rule for target.
adapt-exact-quadr-H1-solver2/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-quadr-H1-solver2/CMakeFiles/adapt-exact-quadr-H1-solver2.dir/build.make tests/adapt-exact-quadr-H1-solver2/CMakeFiles/adapt-exact-quadr-H1-solver2.dir/build
.PHONY : adapt-exact-quadr-H1-solver2/fast

#=============================================================================
# Target rules for targets named adapt-exact-sin-H1

# Build rule for target.
adapt-exact-sin-H1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-sin-H1
.PHONY : adapt-exact-sin-H1

# fast build rule for target.
adapt-exact-sin-H1/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/build.make tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/build
.PHONY : adapt-exact-sin-H1/fast

#=============================================================================
# Target rules for targets named adapt-exact-system-sin-H1

# Build rule for target.
adapt-exact-system-sin-H1: cmake_check_build_system
	$(MAKE) $(MAKESILENT) -f CMakeFiles/Makefile2 adapt-exact-system-sin-H1
.PHONY : adapt-exact-system-sin-H1

# fast build rule for target.
adapt-exact-system-sin-H1/fast:
	$(MAKE) $(MAKESILENT) -f tests/adapt-exact-system-sin-H1/CMakeFiles/adapt-exact-system-sin-H1.dir/build.make tests/adapt-exact-system-sin-H1/CMakeFiles/adapt-exact-system-sin-H1.dir/build
.PHONY : adapt-exact-system-sin-H1/fast

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... install/strip"
	@echo "... list_install_components"
	@echo "... rebuild_cache"
	@echo "... test"
	@echo "... doc"
	@echo "... doc-tex"
	@echo "... test-quick"
	@echo "... _hermes_common"
	@echo "... adapt-exact-quadr-H1-solver0"
	@echo "... adapt-exact-quadr-H1-solver1"
	@echo "... adapt-exact-quadr-H1-solver2"
	@echo "... adapt-exact-quadr-L2"
	@echo "... adapt-exact-sin-H1"
	@echo "... adapt-exact-sin-L2"
	@echo "... adapt-exact-system-sin-H1"
	@echo "... adapt-exact-system-sin-L2"
	@echo "... cpp-callbacks"
	@echo "... first_order_general"
	@echo "... first_order_general_adapt"
	@echo "... first_order_general_adapt_new"
	@echo "... h1d_wrapper"
	@echo "... hermes1d-debug"
	@echo "... hermes_common"
	@echo "... laplace_bc_dirichlet"
	@echo "... laplace_bc_dirichlet_adapt"
	@echo "... laplace_bc_dirichlet_adapt_goal"
	@echo "... laplace_bc_dirichlet_adapt_new"
	@echo "... laplace_bc_neumann"
	@echo "... laplace_bc_newton"
	@echo "... laplace_bc_newton2"
	@echo "... leaks"
	@echo "... legendre-1"
	@echo "... legendre-2"
	@echo "... legendre-3"
	@echo "... lobatto-1"
	@echo "... lobatto-2"
	@echo "... matrix"
	@echo "... matrix-io"
	@echo "... mv"
	@echo "... my_api"
	@echo "... my_api2"
	@echo "... neutronics"
	@echo "... numpy-in-cpp"
	@echo "... python-in-cpp"
	@echo "... schroedinger"
	@echo "... solvers"
	@echo "... sparse"
	@echo "... spblas"
	@echo "... system_car_model"
	@echo "... system_chaotic"
	@echo "... system_exp"
	@echo "... system_neutronics_eigenvalue"
	@echo "... system_neutronics_fixedsrc"
	@echo "... system_neutronics_fixedsrc2"
	@echo "... system_pendulum"
	@echo "... system_sin"
	@echo "... system_sin_adapt"
	@echo "... system_sin_eigen"
	@echo "... timer-1"
	@echo "... vector"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

