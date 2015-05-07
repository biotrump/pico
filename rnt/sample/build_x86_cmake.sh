#!/bin/bash
#Thomas Tsai <thomas@biotrump.com>
#temp for my workspace

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

if [ ! -d ${PICO_OUT} ]; then
	mkdir -p ${PICO_OUT}
else
	rm -rf ${PICO_OUT}/*
fi
pushd ${PICO_OUT}

cmake -DHAVE_SSE=1 -DTARGET_ARCH=${TARGET_ARCH} \
${PICORT_HOME}

#-DLAPACK_SRC:FILEPATH=${LAPACK_SRC} -DLAPACK_BUILD:FILEPATH=${LAPACK_BUILD} \
#-DLAPACK_LIB:FILEPATH=${LAPACK_LIB} -DLAPACKE_SRC:FILEPATH=${LAPACKE_SRC} \
#-DCBLAS_SRC:FILEPATH=${CBLAS_SRC} ..

make ${MAKE_FLAGS}

popd
