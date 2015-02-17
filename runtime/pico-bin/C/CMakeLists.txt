cmake_minimum_required(VERSION 2.8)
PROJECT( cam-play)

#message("V4L2_LIB_DIR is ${V4L2_LIB_DIR}")
#message("V4L2_LIB_OUT is ${V4L2_LIB_OUT}")
#finding lib"v4l2-lib".a, only "v4l2-lib" is used to search!
#find_library(V4L2_LIBS NAMES libv4l2-lib
#			HINTS ${V4L2_LIB_OUT})
#if(V4L2_LIBS-NOTFOUND)
#message(FA )
#endif()
#message( ">>>>LIB::${V4L2_LIBS}" )

find_path(V4L2_INC capture-v4l2.h
			HINTS ${V4L2_LIB_DIR} )
message(">>> INC::${V4L2_INC}" )
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
#message("DSP_ICA_DIR ::${DSP_ICA_DIR}")
find_path(DSP_ICA_INCLUDES dsp-ica.h
			HINTS ${DSP_ICA_DIR} )
message(">>> INC::${DSP_ICA_INCLUDES}" )
include_directories("${DSP_ICA_INCLUDES}")

ADD_EXECUTABLE( cam-play
	cam-play.c
	)
#message(">>> DSP_ICA_LIBS::${DSP_ICA_LIBS}" )
#dynamic or static link
TARGET_LINK_LIBRARIES( cam-play ${OpenCV_LIBS}
	${V4L2_LIB_OUT}/libv4l2-lib.a
	${DSP_ICA_LIBS}/libdsp-ica.a)

#http://www.cmake.org/Wiki/CMake/Tutorials/Exporting_and_Importing_Targets
#Importing Targets
#add_library(ffts-lib STATIC IMPORTED)
#set_property(TARGET ffts-lib PROPERTY IMPORTED_LOCATION ffts/lib/libffts.x86.a)
#add_executable(myexe src1.c src2.c)
#target_link_libraries(myexe ffts-lib)
#
message("FFTS_DIR=${FFTS_DIR}")
#message("CMAKE_FIND_LIBRARY_PREFIXES=${CMAKE_FIND_LIBRARY_PREFIXES}")
#message("CMAKE_FIND_LIBRARY_SUFFIXES=${CMAKE_FIND_LIBRARY_SUFFIXES}")
#PATH_SUFFIXES : the subfolder to find the lib
#UNSET is quite important to clear the CACHE VAR in find_library
#otherwise the CACHE VAR is kept forever no matter libffts.a is killed forever!
UNSET(FFTSLIB CACHE)
find_library (FFTSLIB
			NAMES ffts
			PATHS ${FFTS_DIR}
			PATH_SUFFIXES lib
			NO_DEFAULT_PATH)
message("FFTSLIB=${FFTSLIB}")
ADD_LIBRARY(ffts-lib STATIC IMPORTED)
SET_PROPERTY(TARGET ffts-lib PROPERTY IMPORTED_LOCATION ${FFTSLIB})

include_directories("${FFTS_DIR}/include")

ADD_EXECUTABLE( pulse-v4l2
	pulse-v4l2.c
	)
#message(">>> DSP_ICA_LIBS::${DSP_ICA_LIBS}" )
#dynamic or static link
#	${V4L2_LIB_OUT}/libv4l2-lib.a
TARGET_LINK_LIBRARIES( pulse-v4l2
	${OpenCV_LIBS}
	${DSP_ICA_LIBS}/libdsp-ica.a
	ffts-lib)
