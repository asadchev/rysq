#message(${CUDA_INCLUDE_DIRS})

list(APPEND CUDA_NVCC_FLAGS -std=c++14)

CUDA_ADD_LIBRARY(rysq-cuda kernel.cc kernel.cu)
