# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.6

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin/cmake

# The command to remove a file.
RM = /cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build

# Utility rule file for NightlyCoverage.

# Include the progress variables for this target.
include CMakeFiles/NightlyCoverage.dir/progress.make

CMakeFiles/NightlyCoverage:
	/cvmfs/ilc.desy.de/sw/x86_64_gcc49_sl6/CMake/3.6.3/bin/ctest -D NightlyCoverage

NightlyCoverage: CMakeFiles/NightlyCoverage
NightlyCoverage: CMakeFiles/NightlyCoverage.dir/build.make

.PHONY : NightlyCoverage

# Rule to build all files generated by this target.
CMakeFiles/NightlyCoverage.dir/build: NightlyCoverage

.PHONY : CMakeFiles/NightlyCoverage.dir/build

CMakeFiles/NightlyCoverage.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NightlyCoverage.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NightlyCoverage.dir/clean

CMakeFiles/NightlyCoverage.dir/depend:
	cd /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build /afs/desy.de/user/s/sdharani/flc/eeZHqqbb/Hbb/build/CMakeFiles/NightlyCoverage.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NightlyCoverage.dir/depend

