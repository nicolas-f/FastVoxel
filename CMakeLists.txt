cmake_minimum_required(VERSION 3.15...3.30)

# Maps to a solution file (fastvoxel.sln). The solution will
# have all targets (exe, lib, dll) as projects (.vcproj)
project(fastvoxel VERSION "${SKBUILD_PROJECT_VERSION}")

# Set compiler flags and options.
#############################################################
# enable use of c++11 features where available
# full c++11 support in clang 3.3+: http://clang.llvm.org/cxx_status.html
# for Mac, this is probably Apple LLVM 4.2 (based on LLVM 3.2svn, in XCode 4.6+)
#   or definitely Apple LLVM 5.0 (based on LLVM 3.3svn, in Xcode 5+):
#   https://gist.github.com/yamaya/2924292

IF (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    EXECUTE_PROCESS(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    IF (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        SET(USE_CXX_11 TRUE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    ENDIF()
ELSEIF (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    IF ((NOT APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "3.2")
            OR (APPLE AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.1"))
        SET(USE_CXX_11 TRUE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-error=c++11-narrowing")
    ENDIF()
ELSEIF (MSVC AND MSVC_VERSION GREATER 1600)
    SET(USE_CXX_11 TRUE)
ELSE()
    SET(USE_CXX_11 FALSE)
ENDIF()

# FastVoxel

set(FASTVOXEL_SOURCES
    src/triangle_feeder.cpp
    src/tools/octree44_triangleElement.cpp
    src/std_tools.cpp
    src/spatial_discretization.cpp
    src/scalar_field_creator.cpp
    src/point_feeder.cpp
    src/main_remesh.cpp
    src/input_output/ply/rply.c
    src/input_output/ply/rply_interface.cpp
    src/en_numeric.cpp
    src/Core/mathlib.cpp
    )

# Find Python 3 and numpy
find_package(Python3 COMPONENTS Interpreter Development.Module NumPy REQUIRED)

find_package(
        SWIG
        COMPONENTS python
        REQUIRED)

include(UseSWIG)

list(APPEND CMAKE_SWIG_FLAGS "-py3" "-DPY3")

configure_file(src/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/generated/__init__.py )

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/src")

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(src/fastvoxel.i PROPERTIES CPLUSPLUS ON)

swig_add_library(
        fastvoxel
        LANGUAGE python
        OUTPUT_DIR "${SKBUILD_PLATLIB_DIR}/fastvoxel"
        OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}
        SOURCES src/fastvoxel.i  ${FASTVOXEL_SOURCES})
# OUTPUT_DIR "${SKBUILD_PLATLIB_DIR}"
# Force use of release Python library even in debug builds
if(WIN32 AND CMAKE_BUILD_TYPE MATCHES Debug)
    # Remove the _d suffix that CMake adds for debug builds
    set_target_properties(Python3::Module PROPERTIES
            IMPORTED_LOCATION_DEBUG "${Python3_LIBRARY_RELEASE}"
    )
    # Prevent automatic _d suffix addition
    add_compile_definitions(SWIG_PYTHON_INTERPRETER_NO_DEBUG)
endif()


if(WIN32)
    set_property(TARGET fastvoxel PROPERTY SUFFIX "${Python_SOABI}.pyd")
else()
    set_property(TARGET fastvoxel
            PROPERTY SUFFIX "${Python_SOABI}${CMAKE_SHARED_MODULE_SUFFIX}")
endif()

target_link_libraries(fastvoxel PRIVATE Python3::Module Python3::NumPy)

#--------------#
#    INSTALL
#--------------#

set(FASTVOXEL_RESOURCES ${PROJECT_SOURCE_DIR}/resources/fastvoxel)

# Adds logic to INSTALL.vcproj to copy _fastvoxel.pyd to the destination directory
install (TARGETS fastvoxel DESTINATION fastvoxel)

install (FILES
        ${CMAKE_CURRENT_BINARY_DIR}/generated/__init__.py
        src/np_voxel.py
        DESTINATION fastvoxel)

