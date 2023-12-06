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
include examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/progress.make

# Include the compile flags for this target's objects.
include examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/flags.make

examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o: examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/flags.make
examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o: examples/first_order_general_adapt_new/main.cpp
examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o: examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o"
	cd /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o -MF CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o.d -o CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o -c /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new/main.cpp

examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/first_order_general_adapt_new.dir/main.cpp.i"
	cd /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new/main.cpp > CMakeFiles/first_order_general_adapt_new.dir/main.cpp.i

examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/first_order_general_adapt_new.dir/main.cpp.s"
	cd /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new/main.cpp -o CMakeFiles/first_order_general_adapt_new.dir/main.cpp.s

# Object files for target first_order_general_adapt_new
first_order_general_adapt_new_OBJECTS = \
"CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o"

# External object files for target first_order_general_adapt_new
first_order_general_adapt_new_EXTERNAL_OBJECTS =

examples/first_order_general_adapt_new/first_order_general_adapt_new: examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/main.cpp.o
examples/first_order_general_adapt_new/first_order_general_adapt_new: examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/build.make
examples/first_order_general_adapt_new/first_order_general_adapt_new: src/libhermes1d-debug.so
examples/first_order_general_adapt_new/first_order_general_adapt_new: hermes_common/libhermes_common.so
examples/first_order_general_adapt_new/first_order_general_adapt_new: /usr/bin/python2/library
examples/first_order_general_adapt_new/first_order_general_adapt_new: hermes_common/sparselib/mv/libmv.a
examples/first_order_general_adapt_new/first_order_general_adapt_new: hermes_common/sparselib/libsparse.a
examples/first_order_general_adapt_new/first_order_general_adapt_new: hermes_common/sparselib/spblas/libspblas.a
examples/first_order_general_adapt_new/first_order_general_adapt_new: examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/petr/Repositories/hermes-1d/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable first_order_general_adapt_new"
	cd /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/first_order_general_adapt_new.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/build: examples/first_order_general_adapt_new/first_order_general_adapt_new
.PHONY : examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/build

examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/clean:
	cd /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new && $(CMAKE_COMMAND) -P CMakeFiles/first_order_general_adapt_new.dir/cmake_clean.cmake
.PHONY : examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/clean

examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/depend:
	cd /home/petr/Repositories/hermes-1d && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new /home/petr/Repositories/hermes-1d /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new /home/petr/Repositories/hermes-1d/examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/first_order_general_adapt_new/CMakeFiles/first_order_general_adapt_new.dir/depend

