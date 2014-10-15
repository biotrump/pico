cmake_minimum_required(VERSION 2.8)
PROJECT( pico )

#message("V4L2_LIB_DIR is ${V4L2_LIB_DIR}")
#message("V4L2_LIB_OUT is ${V4L2_LIB_OUT}")
#finding lib"v4l2-lib".a, only "v4l2-lib" is used to search!
#find_library(V4L2_LIBS NAMES v4l2-lib
#			HINTS ${V4L2_LIB_OUT})
#message( "xxx LIB::${V4L2_LIBS}" )

find_path(V4L2_INC capture-v4l2.h
			HINTS ${V4L2_LIB_DIR} )
message("INC:""${V4L2_INC}" )
include_directories("${V4L2_INC}")

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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PICO_OUT}/bin)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${V4L2_CAPTURE_OUT})
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${V4L2_LIB_OUT})

ADD_EXECUTABLE( pico
sample.c
)

#dynamic or static link
TARGET_LINK_LIBRARIES( pico ${OpenCV_LIBS} ${V4L2_LIB_OUT}/libv4l2-lib.a )
