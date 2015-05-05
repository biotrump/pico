cmake_minimum_required(VERSION 2.8)
PROJECT( pico-rt )

message("CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
message("PROJECT_NAME : ${PROJECT_NAME}")
message("pico-rt_SOURCE_DIR : ${pico-rt_SOURCE_DIR}")
message("pico-rt_BINARY_DIR : ${pico-rt_BINARY_DIR}")
# --------------------------------------------------
if("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  add_definitions(-DDEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fPIC")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -g -O0 -fPIC -D_POSIX_C_SOURCE=200112L")
else("${CMAKE_BUILD_TYPE}" MATCHES "Debug")
  add_definitions(-DARMA_NO_DEBUG)
  add_definitions(-DNDEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -fPIC -D_POSIX_C_SOURCE=200112L")
endif("${CMAKE_BUILD_TYPE}" MATCHES "Debug")

# --------------------------------------------------
if(NOT OPENMP_FOUND)
	find_package(OpenMP)
	if(OPENMP_FOUND)
		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	endif(OPENMP_FOUND)
endif(NOT OPENMP_FOUND)
# --------------------------------------------------
message("++++OPENCV_OUT:${OPENCV_OUT}")
message("++++OpenCV_FOUND:${OpenCV_FOUND}")
if(NOT OpenCV_FOUND)
	#if("${OpenCV_DIR}" STREQUAL "")
	set(OpenCV_DIR $ENV{OPENCV_OUT})
	#endif()
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
endif(NOT OpenCV_FOUND)
# --------------------------------------------------
if("${V4L2_LIB_DIR}" STREQUAL "")
	set(V4L2_LIB_DIR $ENV{V4L2_LIB_DIR})
	find_path(V4L2_INC capture-v4l2.h
				HINTS ${V4L2_LIB_DIR} )
	message("V4L2_INC:${V4L2_INC}" )
	if(V4L2_INC)
		include_directories("${V4L2_INC}")
	endif(V4L2_INC)
	#finding libv4l2.a
	if(ANDROID_NDK)
		find_library (V4L2-LIB
			NAMES v4l2-arm
			PATHS ${V4L2_LIB_DIR}/build_and/lib
			PATH_SUFFIXES lib
			NO_DEFAULT_PATH)
	else(ANDROID_NDK)
	#lib"v4l2-x86_64".a
	if("${TARGET_ARCH}" STREQUAL "")
		set(TARGET_ARCH, "x86_64")
	endif()
		find_library (V4L2-LIB
			NAMES v4l2-${TARGET_ARCH}
			PATHS ${V4L2_LIB_OUT}
			PATH_SUFFIXES lib
			NO_DEFAULT_PATH)
	endif(ANDROID_NDK)
	message("V4L2-LIB=${V4L2-LIB} ")
endif("${V4L2_LIB_DIR}" STREQUAL "")

# --------------------------------------------------
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
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${pico-rt_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${pico-rt_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${pico-rt_BINARY_DIR}/lib)
# --------------------------------------------------
#dynamic or static link
#TARGET_LINK_LIBRARIES( pico-rt ${OpenCV_LIBS} ${V4L2_LIB_DIR}/build_x86/lib/libv4l2-x86.a )
if(ANDROID_NDK)

else(ANDROID_NDK)

	ADD_LIBRARY(pico-${TARGET_ARCH} STATIC
		sample.c
		../picornt.c
		)

	ADD_EXECUTABLE( pico-rt
		main.c
		)

	TARGET_LINK_LIBRARIES(pico-rt
		pico-${TARGET_ARCH}
		${OpenCV_LIBS}
		)

	ADD_EXECUTABLE( pico-v4l2
		v4l2.c
		)

	TARGET_LINK_LIBRARIES(pico-v4l2
		pico-${TARGET_ARCH}
		${OpenCV_LIBS}
#	${OpenCV_LIBS} ${V4L2_LIB_DIR}/build_x86/lib/libv4l2-x86.a
		)
	add_dependencies(pico-v4l2
					${OpenCV_LIBS}
					)
	add_dependencies(pico-rt
					${OpenCV_LIBS}
					)
endif(ANDROID_NDK)

