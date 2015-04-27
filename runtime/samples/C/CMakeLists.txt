cmake_minimum_required(VERSION 2.8)
PROJECT( pico )


#opencv should be built first!
set(OpenCV_DIR $ENV{OPENCV_OUT})
message("++++OpenCV_DIR:${OpenCV_DIR}")
FIND_PACKAGE( OpenCV REQUIRED)
if(OpenCV_FOUND)
	message("found")
	message( "OpenCV_INCLUDE_DIRS:${OpenCV_INCLUDE_DIRS}" )
	message( "OpenCV_INCLUDES:${OpenCV_INCLUDES}")
	message( "OpenCV_LIBRARIES:${OpenCV_LIBRARIES}")
	message( "OpenCV_LIBS:${OpenCV_LIBS}")
	message( "OpenCV_DEFINITIONS:${OpenCV_DEFINITIONS}")
else(OpenCV_FOUND)
	message(FATAL_ERROR "openCV not foud")
endif(OpenCV_FOUND)

#openMP to avoid link error : omp_xxx error!
# undefined reference to `GOMP_loop_dynamic_start
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#message("V4L2_LIB_DIR is ${V4L2_LIB_DIR}")
#message("V4L2_LIB_OUT is ${V4L2_LIB_OUT}")
#finding lib"v4l2-lib".a, only "v4l2-lib" is used to search!
#find_library(V4L2_LIBS NAMES v4l2-lib
#			HINTS ${V4L2_LIB_OUT})
#message( "xxx LIB::${V4L2_LIBS}" )

set(V4L2_LIB_DIR $ENV{V4L2_LIB_DIR})
message("V4L2_LIB_DIR is ${V4L2_LIB_DIR}")

find_path(V4L2_INC capture-v4l2.h
			HINTS ${V4L2_LIB_DIR} )
message("INC:${V4L2_INC}" )
include_directories("${V4L2_INC}")

#finding libv4l2.a
find_library (V4L2-LIB
    NAMES v4l2-x86
    PATHS ${V4L2_LIB_DIR}/build_x86/lib
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH)
message("V4L2-LIB=${V4L2-LIB} ")

#There are three kinds of target files that may be built: archive, library, and runtime.
#Executables are always treated as runtime targets.
#Static libraries are always treated as archive targets.
#Module libraries are always treated as library targets.
#For non-DLL platforms shared libraries are treated as library targets.
#For DLL platforms the DLL part of a shared library is treated as a runtime target and
#the corresponding import library is treated as an archive target.
#All Windows-based systems including Cygwin are DLL platforms.
#This property is initialized by the value of the variable CMAKE_RUNTIME_OUTPUT_DIRECTORY
#if it is set when a target is created.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${pico_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${pico_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${pico_BINARY_DIR}/lib)

ADD_EXECUTABLE( pico
sample.c
)

#dynamic or static link
TARGET_LINK_LIBRARIES( pico
				${OpenCV_LIBS}
				${V4L2-LIB})

#
ADD_EXECUTABLE( pico-speed
procspeed.c
)

TARGET_LINK_LIBRARIES( pico
    ${OpenCV_LIBS}
    ${V4L2-LIB} )
