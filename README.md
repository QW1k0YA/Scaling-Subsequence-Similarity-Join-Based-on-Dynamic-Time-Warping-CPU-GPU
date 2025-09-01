# Scaling-Subsequence-Similarity-Join-Based-on-Dynamic-Time-Warping-CPU-GPU

## INTRODUCTION
This repository provides the CPU and GPU implementations of the top-1 subsequence similarity join algorithm based on Dynamic Time Warping (DTW).

## Environment and Tools
- **Operating System**: Windows 11  
- **Hardware**: NVIDIA GeForce RTX 3060 Laptop GPU with 3072 CUDA cores  
- **Development and Debugging Tools**:
  - **CPU**: MinGW, CLion 2023.2.2, CMake  
  - **GPU**: Microsoft Visual Studio (MSVS), CLion 2023.2.2, CMake, CUDA 12.2  

## CPU Version
1. **Parameters**:
   - `argv[1]`: Path to the input data file  
   - `argv[2]`: Output path for execution time  
   - `argv[3]`: Length of the subsequence  
   - `argv[4]`: Ratio of band length to subsequence length  

2. The executable can be run from the command line as:
   ```bash
   CPU_dtw.exe "path/to/data" "path/to/time_output" 256 0.01
   ```

## GPU Version
1. **Parameters**:
   - `argv[1]`: Path to the input data file  
   - `argv[2]`: Output path for execution time  
   - `argv[3]`: Length of the subsequence  
   - `argv[4]`: Ratio of band length to subsequence length  

2. The executable can be run from the command line as:
   ```bash
   GPU_dtw.exe "path/to/data" "path/to/time_output" 256 0.01
   ```

3. **Important Note**:  
   If the subsequence length or band length exceeds 31, you must adjust the macro definition `REGISTER_NUM` at the top of `src/new_dtw_motifGUI_malloc.cu`.  
   Set `REGISTER_NUM` to a value greater than or equal to `ceil(w / 31)`, where `w` is the band length.

   **Example**:  
   For the command:
   ```bash
   GPU_dtw.exe "path/to/data" "path/to/time_output" 1024 0.05
   ```
   The required `REGISTER_NUM` is `ceil(1024 * 0.05 / 31) = ceil(51.2 / 31) = ceil(1.65) = 2`, so `#define REGISTER_NUM 2` is sufficient.

   However, for:
   ```bash
   GPU_dtw.exe "path/to/data" "path/to/time_output" 1024 0.1
   ```
   The band length is `1024 * 0.1 = 102.4`, so `ceil(102.4 / 31) = ceil(3.3) = 4`.  
   Therefore, `REGISTER_NUM` must be set to **4** in this case.
```
