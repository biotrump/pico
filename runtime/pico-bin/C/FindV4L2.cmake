# - Find V4L2
# This module finds if V4L2 Software package is installed
# and determines where the binaries and header files are.
# This code sets the following variables:
#
#  V4L2_FOUND          - True if V4L2 API found
#  V4L2_PATH:          - Path to the V4L2 API folder
#  V4L2_LIBRARY_DIR    - V4L2 libraries folder
#
# Created: 5 Aug 2011 by Marian Zajko (marian.zajko@ximea.com)
# Updated: 25 June 2012 by Igor Kuzmin (parafin@ximea.com)
# Updated: 22 October 2012 by Marian Zajko (marian.zajko@ximea.com)
#

set(V4L2_FOUND)
set(V4L2_PATH)
set(V4L2_LIBRARY_DIR)

if(WIN32)
  # Try to find the V4L2 API path in registry.
  GET_FILENAME_COMPONENT(V4L2_PATH "[HKEY_CURRENT_USER\\Software\\V4L2\\CamSupport\\API;Path]" ABSOLUTE)

  if(EXISTS ${V4L2_PATH})
    set(V4L2_FOUND 1)
    # set LIB folders
    if(X86_64)
      set(V4L2_LIBRARY_DIR "${V4L2_PATH}/x64")
    else()
      set(V4L2_LIBRARY_DIR "${V4L2_PATH}/x86")
    endif()
  else()
    set(V4L2_FOUND 0)
  endif()
elseif(APPLE)
  if(EXISTS /Library/Frameworks/m3api.framework)
    set(V4L2_FOUND 1)
  else()
    set(V4L2_FOUND 0)
  endif()
else()
  if(EXISTS /home/thomas/build/misc/v4l2-lib)
    set(V4L2_FOUND 1)
    # set folders
    set(V4L2_PATH /home/thomas/build/misc/v4l2-lib)
    set(V4L2_LIBRARY_DIR "${V4L2_PATH}/build")
  else()
    set(V4L2_FOUND 0)
  endif()
endif()

mark_as_advanced(FORCE V4L2_FOUND)
mark_as_advanced(FORCE V4L2_PATH)
mark_as_advanced(FORCE V4L2_LIBRARY_DIR)
