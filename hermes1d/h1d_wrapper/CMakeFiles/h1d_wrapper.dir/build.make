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
include hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/compiler_depend.make

# Include the progress variables for this target.
include hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/progress.make

# Include the compile flags for this target's objects.
include hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/flags.make

hermes1d/h1d_wrapper/h1d_wrapper.cpp: hermes1d/h1d_wrapper/h1d_wrapper.pyx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Cythonizing h1d_wrapper.pyx"
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && cython --cplus -I . -I../../../hermes_common/ -o h1d_wrapper.cpp /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper/h1d_wrapper.pyx

hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o: hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/flags.make
hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o: hermes1d/h1d_wrapper/h1d_wrapper.cpp
hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o: hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o"
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o -MF CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o.d -o CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o -c /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper/h1d_wrapper.cpp

hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.i"
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper/h1d_wrapper.cpp > CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.i

hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.s"
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper/h1d_wrapper.cpp -o CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.s

# Object files for target h1d_wrapper
h1d_wrapper_OBJECTS = \
"CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o"

# External object files for target h1d_wrapper
h1d_wrapper_EXTERNAL_OBJECTS =

hermes1d/h1d_wrapper/h1d_wrapper.so: hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/h1d_wrapper.cpp.o
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/build.make
hermes1d/h1d_wrapper/h1d_wrapper.so: src/libhermes1d-debug.so
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes_common/libhermes_common.so
hermes1d/h1d_wrapper/h1d_wrapper.so: /usr/bin/python2/library
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes_common/sparselib/mv/libmv.a
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes_common/sparselib/libsparse.a
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes_common/sparselib/spblas/libspblas.a
hermes1d/h1d_wrapper/h1d_wrapper.so: hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library h1d_wrapper.so"
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/h1d_wrapper.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/build: hermes1d/h1d_wrapper/h1d_wrapper.so
.PHONY : hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/build

hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/clean:
	cd /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper && $(CMAKE_COMMAND) -P CMakeFiles/h1d_wrapper.dir/cmake_clean.cmake
.PHONY : hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/clean

hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/depend: hermes1d/h1d_wrapper/h1d_wrapper.cpp
	cd /home/petr/Repositories/hermes-1d && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper /home/petr/Repositories/hermes-1d/hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : hermes1d/h1d_wrapper/CMakeFiles/h1d_wrapper.dir/depend
