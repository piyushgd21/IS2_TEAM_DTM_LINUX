cmake_minimum_required(VERSION 3.10)
project(featureExtraction)

# ====================
# C++ Standard
# ====================
set(CMAKE_CXX_STANDARD 14 CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ====================
# Compiler Flags
# ====================
#f (MSVC)
 #  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4 /std:c++14 / bigobj/ permissive")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
 #  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
#  else()
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#   set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
# endif()


if (MSVC)
    add_compile_options(
        /std:c++14
        /bigobj
        /permissive
        /W4
    )
endif()
# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Ensure /Zi is used to generate .pdb (debug info in Release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /O2 /Zi")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG")

# Avoid _ITERATOR_DEBUG_LEVEL mismatch (Boost is built in release mode)
add_compile_definitions(_ITERATOR_DEBUG_LEVEL=0)

# ====================
# Manual Include Paths (Avoiding Conda Eigen)
# ====================

# Eigen (manually installed with PCL)
include_directories("C:/Program Files/PCL 1.8.1/3rdParty/Eigen/include/eigen3")

# Gflags + Glog
set(GFLAGS_INCLUDE_DIR "C:/ThirdParty/install2/include")
set(GFLAGS_LIBRARY "C:/ThirdParty/install2/lib/gflags_nothreads_static.lib")
set(GLOG_INCLUDE_DIR "C:/ThirdParty/install2/include")
set(GLOG_LIBRARY "C:/ThirdParty/install2/lib/glog.lib")

# Boost
set(Boost_INCLUDE_DIR "C:/Program Files/PCL 1.8.1/3rdParty/Boost/include/boost-1_64")
set(Boost_LIBRARY_DIR "C:/Program Files/PCL 1.8.1/3rdParty/Boost/lib")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

# PCL
set(PCL_DIR "C:/Program Files/PCL 1.8.1/cmake")

# PDAL (OSGeo4W)
set(PDAL_INCLUDE_DIR "C:/OSGeo4W/include")
set(PDAL_LIBRARY "C:/OSGeo4W/lib/pdalcpp.lib")
set(PDAL_DIR "C:/OSGeo4W/lib/cmake/PDAL")

# Prefix path to help CMake locate PDALConfig.cmake
set(CMAKE_PREFIX_PATH "${PDAL_DIR};C:/Program Files/PCL 1.8.1;C:/ThirdParty/install2")

# ====================
# Find Packages
# ====================
# Do NOT find_package(Eigen3) — we use manually installed Eigen
find_package(PCL REQUIRED)
find_package(Ceres REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(PDAL REQUIRED CONFIG)

# ====================
# Include Directories
# ====================
include_directories(
    ${PCL_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${GFLAGS_INCLUDE_DIR}
    ${GLOG_INCLUDE_DIR}
    ${PDAL_INCLUDE_DIR}
)

# ====================
# Add Executable
# ====================
add_executable(featureExtraction
    src/CTrajectory.cpp
    src/CMapTree.cpp
    src/IntegratedScan.cpp
    src/utility.cpp
    src/featureExtraction.cpp
    src/cylinder_estimator.cpp
    src/cylinder_model_fit.cpp
    src/Mapping.cpp
    src/main.cpp
)

# ====================
# Link Libraries
# ====================
target_link_libraries(featureExtraction
    ${Boost_LIBRARIES}
    ${PCL_LIBRARIES}
    ${CERES_LIBRARIES}
    ${GLOG_LIBRARY}
    ${GFLAGS_LIBRARY}
    ${PDAL_LIBRARY}
)

# ====================
# Force compiler standard (safety net)
# ====================
target_compile_options(featureExtraction PRIVATE "/std:c++14")
