# Scaling-Subsequence-Similarity-Join-Based-on-Dynamic-Time-Warping-CPU-GPU
## INTRODUCTION
   The repository privide CPU and GPU code of the top-1 subsequence similarity join based on Dynamic Time Warping.
## Environment and tools
+ Software environment: windows 11
+ Hardware: NVIDIA GeForce RTX 3060 laptop graphic card, which has 30 Stream Processing Units.
+ Development and debugging tools
   + CPU:MinGW;Clion 2023.2.2;Cmake
   + GPU:MSVS;Clion 2023.2.2;Cmake;cuda 12.2
## CPU_version
1. **parameters** :
argv[1]:path to the data file
argv[2]:output path of execute time
argv[3]:len of subsequence
argv[4]:The ratio of band length to subsequence length  
2. The exe program can be executed by the cmd prompt of
```
CPU_dtw.exe "path of data" "path of time" 256 0.01
```
## GPU_version
1. **parameters** :
argv[1]:path to the data file
argv[2]:output path of execute time
argv[3]:len of subsequence
argv[4]:The ratio of band length to subsequence length  
2. The exe program can be executed by the cmd prompt of
```
GPU_dtw.exe "path of data" "path of time" 256 0.01
```
3. If the length of subsequence and the length of the band is over 31,it is needed to change the Macro definition of `REGISTER_NUM` in the top of `src/new_dtw_motifGUI_malloc.cu` to the proper size, which must be over or equel to `ceil(w/31)`[w:the length of band].  
for example: 
```
GPU_dtw.exe "path of data" "path of time" 1024 0.05
```
is ok for `#define REGISTER_NUM 2` because ceil(1024*0.05/31) = 2 ,but is not suitable for the following command:
```
GPU_dtw.exe "path of data" "path of time" 1024 0.1
```
In this situation, `REGISTER_NUM` is needed to be set to 4;
