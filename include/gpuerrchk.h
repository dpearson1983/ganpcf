#ifndef _GPUERRCHK_H_
#define _GPUERRCHK_H_

#include <cuda.h>
#include <cublas_v2.h>
#include <sstream>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
    if (code != cudaSuccess) {
        std::stringstream err;
        err << "GPUassert -- " << file << "(" << line << "): " << cudaGetErrorString(code) << std::endl;
        throw std::runtime_error(err.str());
    }
} 

// inline void gpuAssert(cublasStatus_t code, const char *file, int line) {
//     if (code != cudaSuccess) {
//         std::stringstream err;
//         err << "GPUassert -- " << file << "(" << line << "): " << cudaGetErrorString(code) << std::endl;
//         throw std::runtime_error(err.str());
//     }
// } 

#endif
