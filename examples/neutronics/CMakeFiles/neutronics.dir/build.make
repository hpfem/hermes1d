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
include examples/neutronics/CMakeFiles/neutronics.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/neutronics/CMakeFiles/neutronics.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/neutronics/CMakeFiles/neutronics.dir/progress.make

# Include the compile flags for this target's objects.
include examples/neutronics/CMakeFiles/neutronics.dir/flags.make

examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o: examples/neutronics/CMakeFiles/neutronics.dir/flags.make
examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o: examples/neutronics/main.cpp
examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o: examples/neutronics/CMakeFiles/neutronics.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o"
	cd /home/petr/Repositories/hermes-1d/examples/neutronics && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o -MF CMakeFiles/neutronics.dir/main.cpp.o.d -o CMakeFiles/neutronics.dir/main.cpp.o -c /home/petr/Repositories/hermes-1d/examples/neutronics/main.cpp

examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/neutronics.dir/main.cpp.i"
	cd /home/petr/Repositories/hermes-1d/examples/neutronics && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/petr/Repositories/hermes-1d/examples/neutronics/main.cpp > CMakeFiles/neutronics.dir/main.cpp.i

examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/neutronics.dir/main.cpp.s"
	cd /home/petr/Repositories/hermes-1d/examples/neutronics && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/petr/Repositories/hermes-1d/examples/neutronics/main.cpp -o CMakeFiles/neutronics.dir/main.cpp.s

# Object files for target neutronics
neutronics_OBJECTS = \
"CMakeFiles/neutronics.dir/main.cpp.o"

# External object files for target neutronics
neutronics_EXTERNAL_OBJECTS =

examples/neutronics/neutronics: examples/neutronics/CMakeFiles/neutronics.dir/main.cpp.o
examples/neutronics/neutronics: examples/neutronics/CMakeFiles/neutronics.dir/build.make
examples/neutronics/neutronics: src/libhermes1d-debug.so
examples/neutronics/neutronics: hermes_common/libhermes_common.so
examples/neutronics/neutronics: /usr/bin/python2/library
examples/neutronics/neutronics: hermes_common/sparselib/mv/libmv.a
examples/neutronics/neutronics: hermes_common/sparselib/libsparse.a
examples/neutronics/neutronics: hermes_common/sparselib/spblas/libspblas.a
examples/neutronics/neutronics: examples/neutronics/CMakeFiles/neutronics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable neutronics"
	cd /home/petr/Repositories/hermes-1d/examples/neutronics && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/neutronics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/neutronics/CMakeFiles/neutronics.dir/build: examples/neutronics/neutronics
.PHONY : examples/neutronics/CMakeFiles/neutronics.dir/build

examples/neutronics/CMakeFiles/neutronics.dir/clean:
	cd /home/petr/Repositories/hermes-1d/examples/neutronics && $(CMAKE_COMMAND) -P CMakeFiles/neutronics.dir/cmake_clean.cmake
.PHONY : examples/neutronics/CMakeFiles/neutronics.dir/clean

examples/neutronics/CMakeFiles/neutronics.dir/depend:
	cd /home/petr/Repositories/hermes-1d && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/examples/neutronics /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/examples/neutronics /home/petr/Repositories/hermes-1d/examples/neutronics/CMakeFiles/neutronics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/neutronics/CMakeFiles/neutronics.dir/depend

