cmake_minimum_required(VERSION 3.24)
project(GPU_dtw LANGUAGES C CXX CUDA)


set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CUDA_COMPILER "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.2/bin/nvcc.exe")

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O2")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}  -rdc=true -arch=sm_86")
set(CMAKE_BUILD_PARALLEL_LEVEL 4)

set(SOURCES_KERNEL_DIAG lib/GLOBAL_DIAG.cu lib/diag_mask_global.cu lib/diag_fast.cu)
add_library(obj_kernel_DIAG OBJECT ${SOURCES_KERNEL_DIAG})
target_compile_options(obj_kernel_DIAG PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:-O2>
        $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>
        $<$<COMPILE_LANGUAGE:CUDA>:-g>
        $<$<COMPILE_LANGUAGE:CUDA>:-maxrregcount=96>
        $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas;-v>

)


file(GLOB_RECURSE CUDA_SOURCES
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu"
        "${CMAKE_CURRENT_SOURCE_DIR}/allmatrixop/*.cu"
)



add_executable(GPU_dtw  ${CUDA_SOURCES}
)

set_target_properties(GPU_dtw PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
)

target_link_libraries(GPU_dtw PRIVATE cudart cudadevrt ${CUDAToolkit_LIBRARIES} obj_kernel_DIAG)
target_include_directories(GPU_dtw PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/allunder
        ${CMAKE_CURRENT_SOURCE_DIR}/alldef
        ${CMAKE_CURRENT_SOURCE_DIR}/alldtw
)
target_include_directories(obj_kernel_DIAG PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/allunder
        ${CMAKE_CURRENT_SOURCE_DIR}/alldef
        ${CMAKE_CURRENT_SOURCE_DIR}/alldtw
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
        # 如果 GPU_parameters.h 在其他目录，也加上
)
