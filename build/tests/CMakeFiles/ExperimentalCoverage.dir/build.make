# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/asharm1/Documents/lpm-v2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/asharm1/Documents/lpm-v2/build

# Utility rule file for ExperimentalCoverage.

# Include any custom commands dependencies for this target.
include tests/CMakeFiles/ExperimentalCoverage.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/ExperimentalCoverage.dir/progress.make

tests/CMakeFiles/ExperimentalCoverage:
	cd /Users/asharm1/Documents/lpm-v2/build/tests && /Applications/CMake.app/Contents/bin/ctest -D ExperimentalCoverage

ExperimentalCoverage: tests/CMakeFiles/ExperimentalCoverage
ExperimentalCoverage: tests/CMakeFiles/ExperimentalCoverage.dir/build.make
.PHONY : ExperimentalCoverage

# Rule to build all files generated by this target.
tests/CMakeFiles/ExperimentalCoverage.dir/build: ExperimentalCoverage
.PHONY : tests/CMakeFiles/ExperimentalCoverage.dir/build

tests/CMakeFiles/ExperimentalCoverage.dir/clean:
	cd /Users/asharm1/Documents/lpm-v2/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/ExperimentalCoverage.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/ExperimentalCoverage.dir/clean

tests/CMakeFiles/ExperimentalCoverage.dir/depend:
	cd /Users/asharm1/Documents/lpm-v2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/asharm1/Documents/lpm-v2 /Users/asharm1/Documents/lpm-v2/tests /Users/asharm1/Documents/lpm-v2/build /Users/asharm1/Documents/lpm-v2/build/tests /Users/asharm1/Documents/lpm-v2/build/tests/CMakeFiles/ExperimentalCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/ExperimentalCoverage.dir/depend
