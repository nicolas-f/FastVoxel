# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel

# Utility rule file for fastvoxel_swig_compilation.

# Include the progress variables for this target.
include CMakeFiles/fastvoxel_swig_compilation.dir/progress.make

CMakeFiles/fastvoxel_swig_compilation: CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON.stamp


CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON.stamp: src/fastvoxel.i
CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON.stamp: src/fastvoxel.i
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Swig compile src/fastvoxel.i for python"
	/usr/bin/cmake -E make_directory /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/CMakeFiles/_fastvoxel.dir
	/usr/bin/cmake -E touch /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON.stamp
	/usr/bin/cmake -E env SWIG_LIB=/usr/local/share/swig/4.0.2 /usr/local/bin/swig -python -outdir /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel -c++ -interface _fastvoxel -I/home/qgoestch/.local/lib/python2.7/site-packages/numpy/core/include -I/usr/local/include/python2.7 -I/home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/src -o /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON_wrap.cxx /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/src/fastvoxel.i

fastvoxel_swig_compilation: CMakeFiles/fastvoxel_swig_compilation
fastvoxel_swig_compilation: CMakeFiles/_fastvoxel.dir/fastvoxelPYTHON.stamp
fastvoxel_swig_compilation: CMakeFiles/fastvoxel_swig_compilation.dir/build.make

.PHONY : fastvoxel_swig_compilation

# Rule to build all files generated by this target.
CMakeFiles/fastvoxel_swig_compilation.dir/build: fastvoxel_swig_compilation

.PHONY : CMakeFiles/fastvoxel_swig_compilation.dir/build

CMakeFiles/fastvoxel_swig_compilation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fastvoxel_swig_compilation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fastvoxel_swig_compilation.dir/clean

CMakeFiles/fastvoxel_swig_compilation.dir/depend:
	cd /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel /home/qgoestch/Documents/These/Scripts/TLM_python/FastVoxel/CMakeFiles/fastvoxel_swig_compilation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fastvoxel_swig_compilation.dir/depend
