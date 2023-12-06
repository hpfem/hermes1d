# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

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

# Include any dependencies generated for this target.
include tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/progress.make

# Include the compile flags for this target's objects.
include tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/flags.make

tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o: tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/flags.make
tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o: tests/adapt-exact-sin-H1/main.cpp
tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o: tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o"
	cd /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o -MF CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o.d -o CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o -c /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1/main.cpp

tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.i"
	cd /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1/main.cpp > CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.i

tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.s"
	cd /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1/main.cpp -o CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.s

# Object files for target adapt-exact-sin-H1
adapt__exact__sin__H1_OBJECTS = \
"CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o"

# External object files for target adapt-exact-sin-H1
adapt__exact__sin__H1_EXTERNAL_OBJECTS =

tests/adapt-exact-sin-H1/adapt-exact-sin-H1: tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/main.cpp.o
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/build.make
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: src/libhermes1d-debug.so
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: hermes_common/libhermes_common.so
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: /usr/bin/python2/library
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: hermes_common/sparselib/mv/libmv.a
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: hermes_common/sparselib/libsparse.a
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: hermes_common/sparselib/spblas/libspblas.a
tests/adapt-exact-sin-H1/adapt-exact-sin-H1: tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable adapt-exact-sin-H1"
	cd /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/adapt-exact-sin-H1.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/build: tests/adapt-exact-sin-H1/adapt-exact-sin-H1
.PHONY : tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/build

tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/clean:
	cd /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 && $(CMAKE_COMMAND) -P CMakeFiles/adapt-exact-sin-H1.dir/cmake_clean.cmake
.PHONY : tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/clean

tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/depend:
	cd /home/petr/Repositories/hermes-1d && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1 /home/petr/Repositories/hermes-1d/tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/adapt-exact-sin-H1/CMakeFiles/adapt-exact-sin-H1.dir/depend
