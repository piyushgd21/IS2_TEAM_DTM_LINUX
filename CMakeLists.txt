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
if (MSVC)
    add_compile_options(
        /std:c++14
        /bigobj
        /permissive
        /W4
    )
else()
    # GCC/Clang flags for Linux
    add_compile_options(
        -Wall
        -Wextra
        -Wno-unused-parameter
        -Wno-sign-compare
    )
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (MSVC)
    # Ensure /Zi is used to generate .pdb (debug info in Release)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /O2 /Zi")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG")

    # Avoid _ITERATOR_DEBUG_LEVEL mismatch (Boost is built in release mode)
    add_compile_definitions(_ITERATOR_DEBUG_LEVEL=0)
else()
    # Linux optimization flags
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
endif()

# ====================
# Find Packages
# ====================
find_package(PCL REQUIRED)
find_package(Ceres REQUIRED)
find_package(Boost REQUIRED COMPONENTS filesystem)
find_package(PDAL REQUIRED)
find_package(gflags REQUIRED)
find_package(glog REQUIRED)

# ====================
# Include Directories
# ====================
include_directories(
    ${PCL_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIRS}
    ${PDAL_INCLUDE_DIRS}
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
    glog::glog
    gflags
    ${PDAL_LIBRARIES}
)

# ====================
# Force compiler standard (safety net)
# ====================
if (MSVC)
    target_compile_options(featureExtraction PRIVATE "/std:c++14")
else()
    target_compile_features(featureExtraction PRIVATE cxx_std_14)
endif()
