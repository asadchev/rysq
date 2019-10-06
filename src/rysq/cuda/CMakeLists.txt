#message(${CUDA_INCLUDE_DIRS})

set(maxrregcount 32)

list(APPEND CUDA_NVCC_FLAGS -std=c++14 --ptxas-options=--verbose -maxrregcount=${maxrregcount} -keep --source-in-ptx -lineinfo)

CUDA_ADD_LIBRARY(rysq-cuda kernel.cu)
