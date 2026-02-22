## System Requirements

- Ubuntu 18.04 or later (or equivalent Linux distribution)
- At least 8GB RAM recommended
- C++14 compatible compiler (GCC 7+ or Clang 5+)

## Installation Steps

### 1. Install System Dependencies

```bash
sudo apt update
sudo apt install -y build-essential cmake git
```

### 2. Install Required Libraries

```bash
# Point Cloud Library (PCL)
sudo apt install -y libpcl-dev

# Ceres Solver (optimization library)
sudo apt install -y libceres-dev

# Boost libraries
sudo apt install -y libboost-all-dev

# PDAL (Point Data Abstraction Library)
sudo apt install -y libpdal-dev

# Google libraries
sudo apt install -y libgflags-dev libgoogle-glog-dev
```

### 3. Clone and Build the Project

```bash
# Clone the repository
git clone https://github.com/piyushgd21/IS2_TEAM_DTM_LINUX.git
cd IS2_TEAM_DTM_LINUX

# Create build directory
mkdir build
cd build

# Configure with CMake (explicitly set Release build type)
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build the project using modern CMake build command
cmake --build . -j
```

**Alternative traditional build method:**
```bash
# After cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
```

### 4. Verify Build Success

After building, you should have an executable called `featureExtraction` in the build directory:

```bash
ls -la featureExtraction
```

## Data Preparation

The program expects input data in a specific directory structure. Create your input folder with this structure:

```
your_dataset_folder/
├── sequence/
│   ├── setting.txt          # Configuration file (copy from project root)
│   ├── times.txt            # Scan timing information
│   ├── lidar_para.txt       # LiDAR calibration parameters
│   ├── trajectory.txt       # Raw trajectory data
│   ├── vlp/                 # LiDAR scan data (.bin files)
│   │   ├── 7500.bin
│   │   ├── 7501.bin
│   │   └── ...
│   ├── tree_trunk/         # Tree feature files (.txt)
│   │   ├── 7500.txt
│   │   ├── 7501.txt
│   │   └── ...
│   └── ground/             # Ground feature files (.txt)
│       ├── 7500.txt
│       ├── 7501.txt
│       └── ...
```

## Configuration

1. Copy the `setting.txt` from the project root to your `sequence/` folder
2. Modify `setting.txt` parameters as needed for your dataset:
   - `initScanIndex` and `endScanIndex`: Range of scan files to process
   - `channel`: Number of LiDAR channels (32 for Velodyne VLP-32)
   - Other parameters for feature extraction and mapping

## Running the Program

### Create a Batch Process File

Create a text file (e.g., `batch_process.txt`) with the following format:

```
-folder /path/to/your_dataset_folder
-setting setting.txt
```

### Execute the Program

```bash
cd build
./featureExtraction /path/to/batch_process.txt
```

## Output

The program will create an output folder `your_dataset_folder/setting_result/` containing:

- `odometry/` - Odometry results and trajectory information
- `loopclosure/` - Loop closure data
- `TrajectoryGlobal.txt`
- Various log files and debug information
- Mapped point cloud data


## Troubleshooting

### Common Issues:

1. **Missing libraries**: If CMake fails to find libraries, ensure all packages are installed correctly
2. **Compilation errors**: Make sure you have C++14 support (`g++ --version` should show 7.0+)
3. **Memory issues**: Large datasets may require more RAM; try processing smaller scan ranges first
4. **Data format errors**: Ensure your .bin files and .txt feature files match the expected format

### Build Issues:

If you encounter build problems, try:

```bash
# Clean and rebuild
cd build
rm -rf *
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j
```

### Runtime Issues:

- Check that all input files exist and paths are correct
- Verify scan indices in `setting.txt` match your data files
- Ensure sufficient disk space for output files

## Additional Notes

- The program processes LiDAR scans sequentially for odometry and mapping
- Processing time depends on dataset size and hardware
- The system extracts features from trees and ground for SLAM
- Results include trajectory estimation and mapped point clouds

This setup should allow you to build and run the IS2-TEAM LiDAR SLAM system without issues.