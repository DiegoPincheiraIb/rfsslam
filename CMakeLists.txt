
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#
# If the user specifies -DCMAKE_BUILD_TYPE on the command line, take their definition
# and dump it in the cache along with proper documentation, otherwise set CMAKE_BUILD_TYPE
# to Debug prior to calling PROJECT()
#
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
ENDIF()

PROJECT(RFSSLAM)
SET(RFSSLAM_MAJOR_VERSION 1)
SET(RFSSLAM_MINOR_VERSION 1)
SET(RFSSLAM_PATCH_VERSION 0)
SET(RFSSLAM_VERSION ${RFSSLAM_MAJOR_VERSION}.${RFSSLAM_MINOR_VERSION}.${RFSSLAM_PATCH_VERSION})
message(STATUS "RFSSLAM Version ${RFSSLAM_VERSION}")

# Compiler settings
# Note: to change compiler, use for example...
# cmake -D CMAKE_C_COMPILER=gcc CMAKE_CXX_COMPILER=c++ .
#SET(CMAKE_C_COMPILER "gcc")
#SET(CMAKE_CXX_COMPILER "g++")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
   add_compile_options (-fdiagnostics-color=always)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   add_compile_options (-fcolor-diagnostics)
endif ()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#MESSAGE("CMAKE_CXX_COMPILTER_ID = ${CMAKE_CXX_COMPILER_ID}")
IF(APPLE)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-qualifiers ") # For Eigen const warnings
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    MESSAGE("NOTE: Libraries being linked to RFSSLAM must also be compiled using GCC")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-enum-compare ") # For Eigen enum compare warnings
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q ") #required for march=native errors. This forces clang assembler to be used, and will produce some argument unused warnings. 
  ENDIF()
ENDIF()

# OpenMP
IF(APPLE AND ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang"))
  SET(USE_OPENMP OFF) # OpenMP does not yet work with Clang, but hopefully soon
ELSE()
  SET(USE_OPENMP ON CACHE BOOL "Use OpenMP to compile multithreaded versions of SLAM algorithms.")
ENDIF()
IF(USE_OPENMP)
  IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp ")
  ELSE()
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
  ENDIF()
ENDIF()

#Threads

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# Boost
FIND_PACKAGE(Boost COMPONENTS timer chrono system filesystem graph program_options REQUIRED)
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME ON)

# For cmake custom library searches
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Eigen
FIND_PACKAGE(Eigen REQUIRED)

# Google Ceres

FIND_PACKAGE(Ceres REQUIRED)

# G2O 

FIND_PACKAGE(g2o REQUIRED)

# Google perftools
FIND_PACKAGE(Perftools)
IF(Perftools_FOUND)
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(USE_CPU_PROFILER OFF CACHE BOOL "Use Google Perftools for CPU profiling.")
    SET(USE_HEAP_PROFILER OFF CACHE BOOL "Use Google Perftools for heap profiling.")
  ELSE()
    SET(USE_CPU_PROFILER OFF) # Does not play nice with Clang
    SET(USE_HEAP_PROFILER OFF)
  ENDIF()
  IF(USE_CPU_PROFILER)
    ADD_DEFINITIONS(-D_PERFTOOLS_CPU)
    IF(APPLE)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_pie ") # Required due to ASLR
    ENDIF()
  ENDIF()
  IF(USE_HEAP_PROFILER)
    ADD_DEFINITIONS(-D_PERFTOOLS_HEAP)
    IF(APPLE)
      SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,-no_pie ") # Required due to ASLR
    ENDIF()
  ENDIF()
ENDIF()

# sanitizer

if (USE_SANITIZER STREQUAL "Address")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address  -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address  -fno-omit-frame-pointer")
    set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}  -fsanitize=address  -fno-omit-frame-pointer")
elseif (USE_SANITIZER STREQUAL "Thread")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=thread")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread")
    set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}  -fsanitize=thread")
elseif (USE_SANITIZER STREQUAL "Undefined")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=undefined")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=undefined")
    set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}  -fsanitize=undefined")
elseif (USE_SANITIZER STREQUAL "Leak")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=leak")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=leak")
    set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}  -fsanitize=leak")
elseif (USE_SANITIZER STREQUAL "Memory")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=memory")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=memory")
    set (CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS}  -fsanitize=memory")
else()
    SET(USE_SANITIZER "disabled" CACHE STRING "set sanitizer type.")
endif()
# VTK visualization
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# GTest
# If using Ubuntu apt-get to install, go to /usr/src/gtest,
# run cmake . and make  . After compiling , 
# move libgtest_main.a and libgtest.a to /usr/lib
enable_testing()
find_package(GTest REQUIRED)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
#MESSAGE("CMAKE_RUNTIME_OUTPUT_DIRECTORY = " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ./lib)
#MESSAGE("CMAKE_LIBRARY_OUTPUT_DIRECTORY = " ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
#MESSAGE("CMAKE_ARCHIVE_OUTPUT_DIRECTORY = " ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin/examples")
FILE(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/bin/test")

#ADD_SUBDIRECTORY()
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(test)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${EIGEN_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GPERFTOOLS_INCLUDE_DIRS})
IF(GTEST_FOUND)
  INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
ENDIF()
#LINK_DIRECTORIES()

ADD_LIBRARY(rfsslam 
  src/TimeStamp.cpp 
  src/Timer.cpp
  src/Frame.cpp
  src/ProcessModel_Odometry1D.cpp 
  src/ProcessModel_Odometry2D.cpp
  src/ProcessModel_Odometry6D.cpp
  src/ProcessModel_Ackerman2D.cpp
  src/MeasurementModel_Rng1D.cpp
  src/MeasurementModel_RngBrg.cpp
  src/MeasurementModel_RngBrg_amplitude.cpp
  src/MeasurementModel_XY.cpp
  src/MeasurementModel_6D.cpp
  src/Visualizer6D.cpp
  src/MeasurementModel_VictoriaPark.cpp
  src/KalmanFilter_RngBrg_amplitude.cpp
  src/KalmanFilter.cpp
  src/HungarianMethod.cpp
  src/MurtyAlgorithm.cpp
  src/BruteForceAssignment.cpp
  src/CostMatrix.cpp
  src/PermutationLexicographic.cpp
  src/MatrixPermanent.cpp
  src/misc/memProfile.cpp
  src/GaussianGenerators.cpp
  )
TARGET_LINK_LIBRARIES(rfsslam ${Boost_LIBRARIES})

ADD_EXECUTABLE(rbphdslam2dSim src/rbphdslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rbphdslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rbphdslam2dCityTrees10000 src/rbphdslamCityTrees10000.cpp)
TARGET_LINK_LIBRARIES(rbphdslam2dCityTrees10000 ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})


ADD_EXECUTABLE(rblmbslam2dCityTrees10000 src/rblmbslamCityTrees10000.cpp)
TARGET_LINK_LIBRARIES(rblmbslam2dCityTrees10000 ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})


ADD_EXECUTABLE(rbphdslam6dSim src/rbphdslam6dSim.cpp)
TARGET_LINK_LIBRARIES(rbphdslam6dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY} ${VTK_LIBRARIES} Threads::Threads)

ADD_EXECUTABLE(rbphdslam2dSim_amplitude src/rbphdslam2dSim_amplitude.cpp)
TARGET_LINK_LIBRARIES(rbphdslam2dSim_amplitude ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(rbcbmemberslam2dSim src/rbcbmemberslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rbcbmemberslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rblmbslam2dSim src/rblmbslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rblmbslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rfspsoslam1dSim src/rfspsoslam1dSim.cpp)
TARGET_LINK_LIBRARIES(rfspsoslam1dSim ${Boost_LIBRARIES} rfsslam ${CERES_LIBRARIES} ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rfspsoslam2dSim src/rfspsoslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rfspsoslam2dSim ${Boost_LIBRARIES} rfsslam ${CERES_LIBRARIES} ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rfsceresslam1dSim src/rfsceresslam1dSim.cpp)
TARGET_LINK_LIBRARIES(rfsceresslam1dSim ${Boost_LIBRARIES} rfsslam ${CERES_LIBRARIES} ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(rfsceresslam2dSim src/rfsceresslam2dSim.cpp)
TARGET_LINK_LIBRARIES(rfsceresslam2dSim ${Boost_LIBRARIES} rfsslam ${CERES_LIBRARIES} ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(fastslam2dSim src/fastslam2dSim.cpp)
TARGET_LINK_LIBRARIES(fastslam2dSim ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY}) 

ADD_EXECUTABLE(rbphdslam_VictoriaPark src/rbphdslam_VictoriaPark.cpp)
TARGET_LINK_LIBRARIES(rbphdslam_VictoriaPark ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(fastslam_VictoriaPark src/fastslam_VictoriaPark.cpp)
TARGET_LINK_LIBRARIES(fastslam_VictoriaPark ${Boost_LIBRARIES} rfsslam ${Perftools_PROFILER_LIBRARY} ${Perftools_TCMALLOC_LIBRARY})

ADD_EXECUTABLE(analysis2dSim src/analysis2dSim.cpp)
TARGET_LINK_LIBRARIES(analysis2dSim ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(spatialIndexTree src/examples/spatialIndexTree.cpp)
TARGET_LINK_LIBRARIES(spatialIndexTree ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(convertLogFiles src/convertLogFiles.cpp)
TARGET_LINK_LIBRARIES(convertLogFiles ${Boost_LIBRARIES})

IF(GTEST_FOUND)

  #ADD_EXECUTABLE(unitTest test/test.cpp src/Measurement.cpp src/Pose.cpp)
  #TARGET_LINK_LIBRARIES(unitTest ${GTEST_LIBRARIES} ${Boost_LIBRARIES} pthread TimeStamp MeasurementModel ProcessModel LinearAssignment)
  #ADD_TEST(unitTest bin/unitTest)
  
  ADD_EXECUTABLE(spatialIndexTreeTest test/SpatialIndexTreeTest.cpp)
  TARGET_LINK_LIBRARIES(spatialIndexTreeTest ${GTEST_LIBRARIES} ${Boost_LIBRARIES} rfsslam pthread)
  ADD_TEST(unitTest spatialIndexTreeTest)

  ADD_EXECUTABLE(matrixPermanentTest test/MatrixPermanentTest.cpp)
  TARGET_LINK_LIBRARIES(matrixPermanentTest ${GTEST_LIBRARIES} rfsslam pthread)
  ADD_TEST(unitTest matrixPermanentTest)

ENDIF()

# Examples

ADD_EXECUTABLE(linearAssignment_CostMatrixPartitioning src/examples/linearAssignment_CostMatrixPartitioning.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_CostMatrixPartitioning
${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(linearAssignment_LexicographicOrdering src/examples/linearAssignment_LexicographicOrdering.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_LexicographicOrdering ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(linearAssignment_MurtyAlgorithm src/examples/linearAssignment_MurtyAlgorithm.cpp)
TARGET_LINK_LIBRARIES(linearAssignment_MurtyAlgorithm ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(ospaError src/examples/ospaError.cpp)
TARGET_LINK_LIBRARIES(ospaError ${Boost_LIBRARIES} rfsslam)


ADD_EXECUTABLE(ospaErrorFromFile src/examples/ospaErrorFromFile.cpp)
TARGET_LINK_LIBRARIES(ospaErrorFromFile ${Boost_LIBRARIES} rfsslam)

ADD_EXECUTABLE(cpu_timer_test src/tests/cpu_timer_info.cpp)
TARGET_LINK_LIBRARIES(cpu_timer_test ${Boost_LIBRARIES})

ADD_EXECUTABLE(jcbb_test src/tests/jcbb_test.cpp)
TARGET_LINK_LIBRARIES(jcbb_test ${Boost_LIBRARIES} rfsslam)

set_target_properties(linearAssignment_CostMatrixPartitioning
	linearAssignment_LexicographicOrdering
	linearAssignment_MurtyAlgorithm 
	ospaError
	spatialIndexTree
	PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/examples"
)

IF(GTEST_FOUND)
  set_target_properties(cpu_timer_test jcbb_test
    spatialIndexTreeTest
    matrixPermanentTest
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/test"
    )
ENDIF()


# Installation

INSTALL(TARGETS rfsslam
  DESTINATION /usr/local/lib/rfsslam
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ
  EXPORT rfsslam-targets)
INSTALL(EXPORT rfsslam-targets DESTINATION /usr/local/lib/rfsslam)
CONFIGURE_FILE(rfsslam-config.cmake.in "rfsslam-config.cmake" @ONLY)
CONFIGURE_FILE(rfsslam-config-version.cmake.in "rfsslam-config-version.cmake" @ONLY)
INSTALL(FILES ${CMAKE_BINARY_DIR}/rfsslam-config.cmake ${CMAKE_BINARY_DIR}/rfsslam-config-version.cmake
  DESTINATION /usr/local/lib/rfsslam)

INSTALL(DIRECTORY include/ DESTINATION /usr/local/include/rfsslam
  FILES_MATCHING PATTERN "*.hpp"
  PERMISSIONS OWNER_READ OWNER_WRITE GROUP_READ WORLD_READ)
