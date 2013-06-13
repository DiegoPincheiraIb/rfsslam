# Keith Leung 2013

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
SET(CMAKE_BUILD_TYPE Debug)
#CMAKE_BUILD_TYPE(CMAKE_CXX_FLAGS_RELEASE)
PROJECT(RB_PHD_FILTER)

FILE(MAKE_DIRECTORY "bin")
FILE(MAKE_DIRECTORY "lib")

FIND_PACKAGE(Eigen REQUIRED)

# Change the following include directories depending on where gtest is installed
# If using Ubuntu apt-get to install, after compiling (this has to be done manually), 
# move libgtest_main.a and libgtest.a to /usr/lib
enable_testing()
find_package(GTest REQUIRED)


SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
MESSAGE("CMAKE_RUNTIME_OUTPUT_DIRECTORY = " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
MESSAGE("CMAKE_ARCHIVE_OUTPUT_DIRECTORY = " ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#ADD_SUBDIRECTORY()
INCLUDE_DIRECTORIES(include)
INCLUDE_DIRECTORIES(src)
INCLUDE_DIRECTORIES(test)
INCLUDE_DIRECTORIES(${Eigen_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})
#LINK_DIRECTORIES()

ADD_LIBRARY(Pose src/Pose.cpp)
ADD_LIBRARY(ProcessModel src/ProcessModel.cpp)
ADD_LIBRARY(Landmark src/Landmark.cpp)
ADD_LIBRARY(Measurement src/Measurement.cpp)
ADD_LIBRARY(MeasurementModel src/MeasurementModel.cpp)

# Is there a point to making this into a library since it's a template?
ADD_LIBRARY(Particle src/Particle.hpp)
SET_TARGET_PROPERTIES(Particle PROPERTIES LINKER_LANGUAGE CXX)

# Is there a point to making this into a library since it's a template?
ADD_LIBRARY(ParticleFilter src/ParticleFilter.hpp)
SET_TARGET_PROPERTIES(ParticleFilter PROPERTIES LINKER_LANGUAGE CXX)

# Is there a point to making this into a library since it's a template?
ADD_LIBRARY(GaussianMixture src/GaussianMixture.hpp)
SET_TARGET_PROPERTIES(GaussianMixture PROPERTIES LINKER_LANGUAGE CXX)

ADD_EXECUTABLE(sandbox src/test.cpp)
TARGET_LINK_LIBRARIES(sandbox Pose Particle ProcessModel Landmark Measurement MeasurementModel GaussianMixture)

ADD_EXECUTABLE(unitTest test/test.cpp)
TARGET_LINK_LIBRARIES(unitTest ${GTEST_LIBRARIES} pthread Landmark Pose Measurement MeasurementModel)
ADD_TEST(unitTest bin/unitTest)
