# Specify the minimum required version of CMake to build the project.
cmake_minimum_required( VERSION 3.4 )
project( facedeformer )

list( APPEND CMAKE_PREFIX_PATH "$ENV{HT}/cmake")
find_package( Houdini REQUIRED )

# Alglib finder with ALGLIB_PREFIX_PATH var:
set( ALGLIB_SRC_PATH "$ENV{ALGLIB_PREFIX_PATH}/src")

if (NOT EXISTS ${ALGLIB_SRC_PATH})
	message(FATAL_ERROR "Specify ALGLIB location with $ALGLIB_SRC_PATH. CMake will exit now.")
else()
	message(STATUS "ALGLIB_SRC_PATH: ${ALGLIB_SRC_PATH}")
endif()

# Collect files.
include_directories( ${ALGLIB_SRC_PATH} )
file( GLOB ALGLIB_SOURCES ${ALGLIB_SRC_PATH}/*.cpp )

# Eigen time
set( EIGEN_INCLUDE_DIR "$ENV{EIGEN_INCLUDE_DIR}")
if (NOT EXISTS ${EIGEN_INCLUDE_DIR})
	message(FATAL_ERROR "Specify EIGEN location with $EIGEN_INCLUDE_DIR. CMake will exit now.")
else()
	message(STATUS "EIGEN_INCLUDE_DIR: ${EIGEN_INCLUDE_DIR}")
endif()
include_directories(${EIGEN_INCLUDE_DIR})

set( CMAKE_BUILD_TYPE Debug)

# Add a library with two source files.
set( library_name SOP_FaceDeform )

add_library( ${library_name} SHARED
	${ALGLIB_SOURCES}
    src/SOP_FaceDeform.cpp
    src/SOP_FaceDeform.hpp
)

# Link against the Houdini libraries, and add required include directories and compile definitions.
target_link_libraries( ${library_name} Houdini )

# Configure several common target properties, such as its output directory.
houdini_configure_target( ${library_name} )