cmake_minimum_required( VERSION 3.1 )
project( featureExtraction )

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS"-std=c++17 -Wall -Wextra")
set(CMAKE_CXX_FLAGS_RELEASE"-O3")

find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
# find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
find_package(PDAL REQUIRED CONFIG)


include_directories(

	${PCL_INCLUDE_DIRS}
  ${CERES_INCLUDE_DIRS}
  # ${OpenCV_INCLUDE_DIRS}
  )

include_directories( "/usr/include/eigen3" )


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS filesystem)
message(STATUS "Boost version: ${Boost_VERSION}")
# This is needed if your Boost version is newer than your CMake version
# or if you have an old version of CMake (<3.5)
if(NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
    set_property(TARGET Boost::filesystem PROPERTY
        INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
endif()


add_executable( featureExtraction src/CTrajectory.cpp src/CMapTree.cpp src/IntegratedScan.cpp src/utility.cpp src/featureExtraction.cpp
src/cylinder_estimator.cpp src/cylinder_model_fit.cpp src/Mapping.cpp src/main.cpp )
target_link_libraries( featureExtraction PUBLIC Boost::filesystem ${PCL_LIBRARIES} ${CERES_LIBRARIES} ${PDAL_LIBRARIES})
target_include_directories(featureExtraction PRIVATE
    ${PDAL_INCLUDE_DIRS}
    ${PDAL_INCLUDE_DIRS}/pdal)