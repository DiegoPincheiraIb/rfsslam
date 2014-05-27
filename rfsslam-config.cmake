# - Config file for the FooBar package
# It defines the following variables
#  FOOBAR_INCLUDE_DIRS - include directories for FooBar
#  FOOBAR_LIBRARIES    - libraries to link against
#  FOOBAR_EXECUTABLE   - the bar executable
 
FIND_PACKAGE(Eigen REQUIRED)
FIND_PACKAGE(Boost COMPONENTS math_c99 timer system thread filesystem graph REQUIRED)

# Compute paths
get_filename_component(RFSSLAM_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

find_path(RFSSLAM_INCLUDE_DIR "rfsslam" PATHS "/usr/local/include")
find_path(RFSSLAM_LIBRARY "rfsslam/librfsslam.a" PATHS "/usr/local/lib")
set(RFSSLAM_INCLUDE_DIRS ${RFS_INCLIDE_DIR} ${Boost_INCLUDE_DIRS} ${Eigen_INCLUDE_DIRS})


# Our library dependencies (contains definitions for IMPORTED targets)
if(NOT TARGET rfsslam AND NOT RFSSLAM_BINARY_DIR)
  include("${RFSSLAM_CMAKE_DIR}/rfsslam-targets.cmake")
endif()
 
# These are IMPORTED targets created by rfsslam-targets.cmake
set(RFSSLAM_LIBRARIES rfsslam ${Boost_LIBRARIES} ${Eigen_LIBRARIES})
#set(RFSSLAM_EXECUTABLE foobar)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(RFSSLAM DEFAULT_MSG RFSSLAM_INCLUDE_DIR RFSSLAM_LIBRARY)
mark_as_advanced(RFSSLAM_INCLUDE_DIR RFSSLAM_LIBRARY)

if(RFSSLAM_FOUND)
  message(STATUS "RFSSLAM found (include: ${RFSSLAM_INCLUDE_DIR})")
endif(RFSSLAM_FOUND)