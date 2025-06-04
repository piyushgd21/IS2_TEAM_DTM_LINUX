cmake_minimum_required(VERSION 3.5...3.27)
project(featureExtraction)
# --- Manually set Boost (installed via vcpkg) ---
set(BOOST_ROOT "C:/Users/deshpa74/vcpkg/installed/x64-windows")
set(Boost_INCLUDE_DIR "${BOOST_ROOT}/include")
set(Boost_LIBRARY_DIR "${BOOST_ROOT}/lib")
set(PCL_DIR "C:/Users/deshpa74/vcpkg/installed/x64-windows/share/pcl")
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(MSVC)
    add_compile_options(/W4)  # Use high warning level for MSVC
else()
    add_compile_options(-Wall -Wextra)
endif()

# Toolchain file for vcpkg (set it only ONCE, preferably via CLI when running CMake)
# set(CMAKE_TOOLCHAIN_FILE "C:/Users/deshpa74/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
set(PDAL_DIR "C:/Users/deshpa74/vcpkg/installed/x64-windows/share/pdal")
# ---- Dependencies ----
find_package(Eigen3 REQUIRED)
find_package(PDAL REQUIRED CONFIG)
find_package(Ceres REQUIRED CONFIG)
find_package(glog CONFIG REQUIRED)
# find_package(Boost 1.88 REQUIRED COMPONENTS filesystem)
find_package(PCL REQUIRED)
# find_package(OpenCV REQUIRED) # Optional if you're using OpenCV later

# ---- Add Executable ----
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
target_compile_definitions(featureExtraction PRIVATE GLOG_USE_GLOG_EXPORT)


target_compile_definitions(featureExtraction PRIVATE GLOG_NO_ABBREVIATED_SEVERITIES)

# ---- Include Directories ----
# Manually add if you need custom locations
# target_include_directories(featureExtraction PRIVATE
#     ${PCL_INCLUDE_DIRS}
#     ${CERES_INCLUDE_DIRS}
#     ${Boost_INCLUDE_DIRS}
#     ${PDAL_INCLUDE_DIRS}
#     C:/libs/eigen-3.4.0
#     C:/ceres-solver-2.2.0/include
#     C:/glog-install/include
# )

# ---- Link Libraries ----
target_link_libraries(featureExtraction PRIVATE
    pcl_common
    pcl_filters
    pcl_io
    pcl_io_ply
    pcl_kdtree
    pcl_octree
    pcl_search
    pcl_sample_consensus
    pcl_features
    pcl_ml
    pcl_segmentation
    pcl_surface
    pcl_registration
    pcl_keypoints
    pcl_tracking
    pcl_recognition
    pcl_stereo
    ceres
    pdalcpp
    ${Boost_LIBRARY_DIR}/boost_filesystem-vc143-mt-x64-1_88.lib
    ${Boost_LIBRARY_DIR}/boost_serialization-vc143-mt-x64-1_88.lib
    ${Boost_LIBRARY_DIR}/boost_iostreams-vc143-mt-x64-1_88.lib

    glog::glog
)

# ---- Optional: Boost Compatibility Patch ----
# Needed only if Boost::filesystem isn't registered
# if(NOT TARGET Boost::filesystem)
#     add_library(Boost::filesystem IMPORTED INTERFACE)
#     set_property(TARGET Boost::filesystem PROPERTY
#         INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
#     set_property(TARGET Boost::filesystem PROPERTY
#         INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
# endif()

# ---- PCL Compatibility (optional cleanup if needed) ----
# If not using vcpkg version of PCL, manually set root:
# set(PCL_ROOT "C:/Program Files/PCL 1.15.0")
