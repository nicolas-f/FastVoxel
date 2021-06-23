CMAKE_MINIMUM_REQUIRED(VERSION 3.6.0)

# Maps to a solution file (fastvoxel.sln). The solution will
# have all targets (exe, lib, dll) as projects (.vcproj)
project (fastvoxel VERSION 1.0.1)

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

message(STATUS " SWIG_EXECUTABLE     - $ENV{SWIG_EXECUTABLE}")
set(CMAKE_INSTALL_MFC_LIBRARIES ON) 
include( InstallRequiredSystemLibraries )
#---------------------------------------#
#    DEPENDENCY & EXECUTABLE (OR LIB)
#---------------------------------------#

# require swig > v3
FIND_PACKAGE(SWIG REQUIRED)


INCLUDE(${SWIG_USE_FILE})

# Find package for building
FIND_PACKAGE(PythonLibs 3.8 REQUIRED)

#Find numpy include directory

find_package(PythonInterp)
set(__numpy_out 1)

if (PYTHON_EXECUTABLE)
  # Find out the include path
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.get_include(), end='')\nexcept:pass\n"
            OUTPUT_VARIABLE __numpy_path)
  # And the version
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c
            "from __future__ import print_function\ntry: import numpy; print(numpy.__version__, end='')\nexcept:pass\n"
    OUTPUT_VARIABLE __numpy_version)
elseif(__numpy_out)
  message(STATUS "Python executable not found.")
endif(PYTHON_EXECUTABLE)


find_path(PYTHON_NUMPY_INCLUDE_DIR numpy/arrayobject.h
  HINTS "${__numpy_path}" "${PYTHON_INCLUDE_PATH}" NO_DEFAULT_PATH)

if(PYTHON_NUMPY_INCLUDE_DIR)
  set(PYTHON_NUMPY_FOUND 1 CACHE INTERNAL "Python numpy found")
endif(PYTHON_NUMPY_INCLUDE_DIR)


INCLUDE_DIRECTORIES(${PYTHON_NUMPY_INCLUDE_DIR})

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/src")

SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(src/fastvoxel.i PROPERTIES CPLUSPLUS ON)


IF ((${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 3.8) OR (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} EQUAL 3.8))
  SWIG_ADD_LIBRARY(fastvoxel LANGUAGE python SOURCES src/fastvoxel.i ${FASTVOXEL_SOURCES})
ELSE()
  SWIG_ADD_MODULE(fastvoxel python src/fastvoxel.i ${FASTVOXEL_SOURCES})
ENDIF ()



SWIG_LINK_LIBRARIES(fastvoxel ${PYTHON_LIBRARIES})

# Creates a folder "libraries" and adds target project (lib_interface.vcproj) under it
set_property(TARGET ${SWIG_MODULE_fastvoxel_REAL_NAME} PROPERTY FOLDER "libraries")

if(WIN32) # Check if we are on Windows
 if(MSVC) # Check if we are using the Visual Studio compiler
   # A 'Debug' build of the wrappers requires a debug build of the Python interpreter.
   # http://www.swig.org/Doc2.0/Python.html#Python_nn12
   # Fix linking with non existing python27_d.lib
   add_definitions(-DSWIG_PYTHON_INTERPRETER_NO_DEBUG)
 endif()
endif()

#--------------#
#    INSTALL
#--------------#

set(FASTVOXEL_RESSOURCES ${PROJECT_SOURCE_DIR}/fastvoxel)

# Adds logic to INSTALL.vcproj to copy _fastvoxel.pyd to the destination directory
install (TARGETS ${SWIG_MODULE_fastvoxel_REAL_NAME}
      LIBRARY DESTINATION fastvoxel/
)

if( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS )
    install( PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} DESTINATION fastvoxel COMPONENT System )
endif( CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS )

install (FILES
        ${FASTVOXEL_RESSOURCES}/__init__.py
		${CMAKE_CURRENT_BINARY_DIR}/fastvoxel.py
        DESTINATION fastvoxel/)
        
        

if(WIN32)
    #install(CODE "
    #include(BundleUtilities)   
    # Rename to .exe as cmake check if the file ends with .exe in windows before gathering dependencies..
    #file(RENAME \"\${CMAKE_INSTALL_PREFIX}/fastvoxel/_fastvoxel.pyd\" \"\${CMAKE_INSTALL_PREFIX}/fastvoxel/_fastvoxel.exe\")    
    #fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/fastvoxel/_fastvoxel.exe\" \"\"   \"${DIRS}\")
    #file(RENAME \"\${CMAKE_INSTALL_PREFIX}/fastvoxel/_fastvoxel.exe\" \"\${CMAKE_INSTALL_PREFIX}/fastvoxel/_fastvoxel.pyd\") 
    #")
endif()
