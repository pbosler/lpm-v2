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

# Include any dependencies generated for this target.
include tests/CMakeFiles/BivarConvergenceTest.exe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include tests/CMakeFiles/BivarConvergenceTest.exe.dir/compiler_depend.make

# Include the progress variables for this target.
include tests/CMakeFiles/BivarConvergenceTest.exe.dir/progress.make

# Include the compile flags for this target's objects.
include tests/CMakeFiles/BivarConvergenceTest.exe.dir/flags.make

tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o: tests/CMakeFiles/BivarConvergenceTest.exe.dir/flags.make
tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o: /Users/asharm1/Documents/lpm-v2/tests/BivarConvergenceTests.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/asharm1/Documents/lpm-v2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o"
	cd /Users/asharm1/Documents/lpm-v2/build/tests && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -c /Users/asharm1/Documents/lpm-v2/tests/BivarConvergenceTests.f90 -o CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o

tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.i"
	cd /Users/asharm1/Documents/lpm-v2/build/tests && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -E /Users/asharm1/Documents/lpm-v2/tests/BivarConvergenceTests.f90 > CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.i

tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.s"
	cd /Users/asharm1/Documents/lpm-v2/build/tests && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -S /Users/asharm1/Documents/lpm-v2/tests/BivarConvergenceTests.f90 -o CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.s

# Object files for target BivarConvergenceTest.exe
BivarConvergenceTest_exe_OBJECTS = \
"CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o"

# External object files for target BivarConvergenceTest.exe
BivarConvergenceTest_exe_EXTERNAL_OBJECTS =

tests/BivarConvergenceTest.exe: tests/CMakeFiles/BivarConvergenceTest.exe.dir/BivarConvergenceTests.f90.o
tests/BivarConvergenceTest.exe: tests/CMakeFiles/BivarConvergenceTest.exe.dir/build.make
tests/BivarConvergenceTest.exe: src/liblpmFortran.a
tests/BivarConvergenceTest.exe: tests/CMakeFiles/BivarConvergenceTest.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/asharm1/Documents/lpm-v2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable BivarConvergenceTest.exe"
	cd /Users/asharm1/Documents/lpm-v2/build/tests && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/BivarConvergenceTest.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
tests/CMakeFiles/BivarConvergenceTest.exe.dir/build: tests/BivarConvergenceTest.exe
.PHONY : tests/CMakeFiles/BivarConvergenceTest.exe.dir/build

tests/CMakeFiles/BivarConvergenceTest.exe.dir/clean:
	cd /Users/asharm1/Documents/lpm-v2/build/tests && $(CMAKE_COMMAND) -P CMakeFiles/BivarConvergenceTest.exe.dir/cmake_clean.cmake
.PHONY : tests/CMakeFiles/BivarConvergenceTest.exe.dir/clean

tests/CMakeFiles/BivarConvergenceTest.exe.dir/depend:
	cd /Users/asharm1/Documents/lpm-v2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/asharm1/Documents/lpm-v2 /Users/asharm1/Documents/lpm-v2/tests /Users/asharm1/Documents/lpm-v2/build /Users/asharm1/Documents/lpm-v2/build/tests /Users/asharm1/Documents/lpm-v2/build/tests/CMakeFiles/BivarConvergenceTest.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : tests/CMakeFiles/BivarConvergenceTest.exe.dir/depend
