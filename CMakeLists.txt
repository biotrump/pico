# ----------------------------------------------------------------------------
# Process subdirectories
# ----------------------------------------------------------------------------
#$ENV{VAR} to get shell env
set(V4L2_LIB_DIR $ENV{V4L2_LIB_DIR})
set(V4L2_LIB_OUT $ENV{V4L2_LIB_OUT} )
add_subdirectory(runtime/samples/C)
add_subdirectory(runtime/pico-bin/C)
add_subdirectory(learning)
