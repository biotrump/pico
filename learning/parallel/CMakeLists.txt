cmake_minimum_required(VERSION 2.8.5)
PROJECT( openmp-test)

#-DCMAKE_BUILD_TYPE:STRING=Debug
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE "Debug")
endif()
message("CMAKE_BUILD_TYPE:" "${CMAKE_BUILD_TYPE}")

#opencv should be built first!
#if(NOT OPENCV_BUILD)
#	set(OpenCV_DIR "~/build/opencv/build/3.0.x-static/")
#else()
#	set(OpenCV_DIR $ENV{OPENCV_BUILD})
#endif()


#message("++++OpenCV_DIR:""${OpenCV_DIR}")
#FIND_PACKAGE(OpenCV REQUIRED)
#if(OpenCV_FOUND)
#	message("found")
#	message( ">>>OpenCV_INCLUDE_DIRS:""${OpenCV_INCLUDE_DIRS}" )
#	message( ">>>OpenCV_INCLUDES:""${OpenCV_INCLUDES}")
#	message( ">>>OpenCV_LIBRARIES:""${OpenCV_LIBRARIES}")
#	message( ">>>OpenCV_LIBS:""${OpenCV_LIBS}")
#	message( ">>>OpenCV_LIBRARY_DIRS:""${OpenCV_LIBRARY_DIRS}")
#	message( ">>>OpenCV_DEFINITIONS:""${OpenCV_DEFINITIONS}")
#else()
#	message(FATAL_ERROR "openCV not foud")
#endif()
#include_directories("${OpenCV_INCLUDE_DIRS}")
#set(OpenCV_LIBRARY_DIRS "${OpenCV_DIR}/lib" CACHE PATH "Full path of opencv library directory")
#link_directories(${OpenCV_LIBRARY_DIRS})
#message( ">>>OpenCV_LIBRARY_DIRS:""${OpenCV_LIBRARY_DIRS}")

#SET(OPENCV_CPP_SAMPLES_REQUIRED_DEPS opencv_core opencv_imgproc opencv_flann
#    opencv_imgcodecs opencv_videoio opencv_highgui opencv_ml opencv_video
#    opencv_objdetect opencv_photo opencv_features2d opencv_calib3d
#    opencv_stitching opencv_videostab opencv_shape)

#openMP to avoid link error : omp_xxx error!
# undefined reference to `GOMP_loop_dynamic_start
find_package(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

ADD_EXECUTABLE( test1
	test1.c
	)

ADD_EXECUTABLE( test2
	test2.c
	)

ADD_EXECUTABLE( test3
	test3.c
	)
TARGET_LINK_LIBRARIES( test3 ${OpenCV_LIBS} )

ADD_EXECUTABLE( genki
	genki.c
	pico-nonface.c
	main.c
	)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PICO_OUT}/bin)

TARGET_LINK_LIBRARIES( genki ${OpenCV_LIBS} )
