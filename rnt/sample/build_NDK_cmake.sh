#!/bin/bash
if [ $# -eq 0 ]; then
	echo "$0 [arm|x86]"
	exit
fi
#export CMAKE_BUILD_TYPE "Debug"
export CMAKE_BUILD_TYPE="Release"

# Make sure you have NDK_ROOT defined in .bashrc or .bash_profile
# Modify INSTALL_DIR to suit your situation
#Lollipop	5.0 - 5.1	API level 21, 22
#KitKat	4.4 - 4.4.4	API level 19
#Jelly Bean	4.3.x	API level 18
#Jelly Bean	4.2.x	API level 17
#Jelly Bean	4.1.x	API level 16
#Ice Cream Sandwich	4.0.3 - 4.0.4	API level 15, NDK 8
#Ice Cream Sandwich	4.0.1 - 4.0.2	API level 14, NDK 7
#Honeycomb	3.2.x	API level 13
#Honeycomb	3.1	API level 12, NDK 6
#Honeycomb	3.0	API level 11
#Gingerbread	2.3.3 - 2.3.7	API level 10
#Gingerbread	2.3 - 2.3.2	API level 9, NDK 5
#Froyo	2.2.x	API level 8, NDK 4

if [ -z "${NDK_ROOT}"  ]; then
	export NDK_ROOT=${HOME}/NDK/android-ndk-r10d
#	export NDK_ROOT=${HOME}/NDK/android-ndk-r9
fi
export ANDROID_NDK=${NDK_ROOT}

if [[ ${NDK_ROOT} =~ .*"-r9".* ]]
then
#ANDROID_APIVER=android-8
#ANDROID_APIVER=android-9
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
#TOOL_VER="4.6"
#gfortran is in r9d V4.8.0
TOOL_VER="4.8.0"
else
#android 4.0.1 ICS and above
ANDROID_APIVER=android-14
TOOL_VER="4.9"
fi

case $(uname -s) in
  Darwin)
    CONFBUILD=i386-apple-darwin`uname -r`
    HOSTPLAT=darwin-x86
    CORE_COUNT=`sysctl -n hw.ncpu`
  ;;
  Linux)
    CONFBUILD=x86-unknown-linux
    HOSTPLAT=linux-`uname -m`
    CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
  ;;
CYGWIN*)
	CORE_COUNT=`grep processor /proc/cpuinfo | wc -l`
	;;
  *) echo $0: Unknown platform; exit
esac

#NDK cross compile toolchains
#default is arm
export ARCHI=$1
echo ARCHI=$ARCHI
case $ARCHI in
  arm)
    TARGPLAT=arm-linux-androideabi
    CONFTARG=arm-eabi
	echo "Using: $NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
	#export PATH="$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:\
	#$NDK_ROOT/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/${TARGPLAT}/bin/:$PATH"
	export PATH="${NDK_ROOT}/toolchains/${TARGPLAT}-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"
  ;;
  x86)
    TARGPLAT=i686-linux-android
    CONFTARG=x86
	echo "Using: $NDK_ROOT/toolchains/x86-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin"
	export PATH="${NDK_ROOT}/toolchains/x86-${TOOL_VER}/prebuilt/${HOSTPLAT}/bin/:$PATH"
#specify assembler for x86 SSE3, but ffts's sse.s needs 64bit x86.
#intel atom z2xxx and the old atoms are 32bit, so 64bit x86 in android can't work in
#most atom devices.
#http://forum.cvapp.org/viewtopic.php?f=13&t=423&sid=4c47343b1de899f9e1b0d157d04d0af1
#	export  CCAS="${TARGPLAT}-as"
#	export  CCASFLAGS="--64 -march=i686+sse3"
#	export  CCASFLAGS="--64"

  ;;
  mips)
  ## probably wrong
    TARGPLAT=mipsel-linux-android
    CONFTARG=mips
  ;;
  *) echo $0: Unknown target; exit
esac
#: ${NDK_ROOT:?}
echo $PATH

export SYS_ROOT="${NDK_ROOT}/platforms/${ANDROID_APIVER}/arch-${ARCHI}/"
export CC="${TARGPLAT}-gcc --sysroot=$SYS_ROOT"
export LD="${TARGPLAT}-ld"
export AR="${TARGPLAT}-ar"
export ARCH=${AR}
export RANLIB="${TARGPLAT}-ranlib"
export STRIP="${TARGPLAT}-strip"
#export CFLAGS="-Os -fPIE"
export CFLAGS="-Os -fPIE --sysroot=$SYS_ROOT"
export CXXFLAGS="-fPIE --sysroot=$SYS_ROOT"
export FORTRAN="${TARGPLAT}-gfortran --sysroot=$SYS_ROOT"

#!!! quite importnat for cmake to define the NDK's fortran compiler.!!!
#Don't let cmake decide it.
#export FC=/home/thomas/aosp/NDK/android-ndk-r9/toolchains/arm-linux-androideabi-4.8.0/prebuilt/linux-x86/bin/arm-linux-androideabi-gfortran
export FC=${FORTRAN}

#include path :
#platforms/android-21/arch-arm/usr/include/

if [ -z "$PICO_DIR" ]; then
	export DSP_HOME=${DSP_HOME:-`pwd`/../../../../dsp}
	export PICO_DIR=${PICO_DIR:-`pwd`/../..}
	export BLIS_DIR=${DSP_HOME}/blis
	export BLIS_DIR=${DSP_HOME}/blis
	export BLISLIB_DIR=${BLISLIB_DIR:-${BLIS_DIR}/lib}
	export FFTS_DIR=${DSP_HOME}/ffts
	export LAPACK_SRC=${LAPACK_SRC:-${DSP_HOME}/LAPACK}
else
	export BLISLIB_DIR=${BLISLIB_DIR:-$BLIS_OUT/lib}
fi
export LAPACKE_SRC=${LAPACKE_SRC:-${LAPACK_SRC}/LAPACKE}
export CBLAS_SRC=${CBLAS_SRC:-${LAPACK_SRC}/CBLAS}
export PICORT_HOME=${PICORT_HOME:-`pwd`}

case $ARCHI in
	arm)
		#armv7a
		#$BLIS_OUT/lib/libblis-NDK-arm.a
		export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-arm}
		#export FFTS_LIB_DIR=${FFTS_OUT}/lib
		#libffts-NDK-arm.a
		FFTS_LIB_NAME=ffts-NDK-${TARGET_ARCH}
		#export LAPACK_BUILD=${LAPACK_BUILD:-${LAPACK_SRC}/build_NDK_arm}
		;;
	x86)
		#x86
		#$BLIS_OUT/lib/libblis-NDK-x86.a
		export BLIS_LIB_NAME=${BLIS_LIB_NAME:-blis-NDK-x86}
		#export FFTS_LIB_DIR=${FFTS_OUT}/lib
		#libffts-NDK-x86.a
		FFTS_LIB_NAME=ffts-NDK-x86
		#export LAPACK_BUILD=${LAPACK_BUILD:-${LAPACK_SRC}/build_NDK_x86}
	;;
	mips)
	;;
	*) echo $0: Unknown target; exit
esac
if [ ! -d ${PICO_OUT} ]; then
	mkdir -p ${PICO_OUT}
else
	rm -rf ${PICO_OUT}/*
fi
pushd ${PICO_OUT}

case $ARCHI in
  arm)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_HOME}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=${TARGPLAT}-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="armeabi-v7a with VFPV3" \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DTARGET_ARCH=${TARGET_ARCH} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_HOME}

	;;
  x86)
	cmake -DCMAKE_TOOLCHAIN_FILE=${PICORT_HOME}/android.toolchain.cmake \
	-DANDROID_NATIVE_API_LEVEL=${ANDROID_APIVER} \
	-DANDROID_NDK=${NDK_ROOT} -DANDROID_TOOLCHAIN_NAME=x86-${TOOL_VER} \
	-DCMAKE_BUILD_TYPE=Release -DANDROID_ABI="x86" \
	-DFFTS_DIR:FILEPATH=${FFTS_DIR} -DFFTS_LIB_NAME=${FFTS_LIB_NAME} \
	-DFFTS_OUT:FILEPATH=${FFTS_OUT} \
	-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
	-DLAPACK_OUT:FILEPATH=${LAPACK_OUT} \
	-DUSE_BLIS=1 -DBLIS_DIR:FILEPATH=${BLIS_DIR} -DBLISLIB_DIR:FILEPATH=${BLISLIB_DIR} \
	-DBLIS_LIB_NAME=${BLIS_LIB_NAME} \
	-DTARGET_ARCH=${TARGET_ARCH} \
	-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ${PICORT_HOME}

  ;;
  mips)
  ;;
  *) echo $0: Unknown target; exit
esac

make -j${CORE_COUNT}

popd

