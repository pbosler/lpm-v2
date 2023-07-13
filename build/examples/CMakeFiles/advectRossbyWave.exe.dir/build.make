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
include examples/CMakeFiles/advectRossbyWave.exe.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/CMakeFiles/advectRossbyWave.exe.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/advectRossbyWave.exe.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/advectRossbyWave.exe.dir/flags.make

examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o: examples/CMakeFiles/advectRossbyWave.exe.dir/flags.make
examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o: /Users/asharm1/Documents/lpm-v2/examples/AdvectRHWave.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/asharm1/Documents/lpm-v2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o"
	cd /Users/asharm1/Documents/lpm-v2/build/examples && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -c /Users/asharm1/Documents/lpm-v2/examples/AdvectRHWave.f90 -o CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o

examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.i"
	cd /Users/asharm1/Documents/lpm-v2/build/examples && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -E /Users/asharm1/Documents/lpm-v2/examples/AdvectRHWave.f90 > CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.i

examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.s"
	cd /Users/asharm1/Documents/lpm-v2/build/examples && /opt/homebrew/bin/mpifort $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -ffree-form -S /Users/asharm1/Documents/lpm-v2/examples/AdvectRHWave.f90 -o CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.s

# Object files for target advectRossbyWave.exe
advectRossbyWave_exe_OBJECTS = \
"CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o"

# External object files for target advectRossbyWave.exe
advectRossbyWave_exe_EXTERNAL_OBJECTS =

examples/advectRossbyWave.exe: examples/CMakeFiles/advectRossbyWave.exe.dir/AdvectRHWave.f90.o
examples/advectRossbyWave.exe: examples/CMakeFiles/advectRossbyWave.exe.dir/build.make
examples/advectRossbyWave.exe: src/liblpmFortran.a
examples/advectRossbyWave.exe: examples/CMakeFiles/advectRossbyWave.exe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/asharm1/Documents/lpm-v2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking Fortran executable advectRossbyWave.exe"
	cd /Users/asharm1/Documents/lpm-v2/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/advectRossbyWave.exe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/advectRossbyWave.exe.dir/build: examples/advectRossbyWave.exe
.PHONY : examples/CMakeFiles/advectRossbyWave.exe.dir/build

examples/CMakeFiles/advectRossbyWave.exe.dir/clean:
	cd /Users/asharm1/Documents/lpm-v2/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/advectRossbyWave.exe.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/advectRossbyWave.exe.dir/clean

examples/CMakeFiles/advectRossbyWave.exe.dir/depend:
	cd /Users/asharm1/Documents/lpm-v2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/asharm1/Documents/lpm-v2 /Users/asharm1/Documents/lpm-v2/examples /Users/asharm1/Documents/lpm-v2/build /Users/asharm1/Documents/lpm-v2/build/examples /Users/asharm1/Documents/lpm-v2/build/examples/CMakeFiles/advectRossbyWave.exe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/advectRossbyWave.exe.dir/depend

