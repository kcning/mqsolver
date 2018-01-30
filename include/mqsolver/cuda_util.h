/* cuda_util.h: define common utilities for CUDA */

#ifndef __CUDA_UTIL_H__
#define __CUDA_UTIL_H__

#include <stdio.h>
#include <stdbool.h>
#include <cuda_runtime.h>


#include "util.h"


/* function: cuda_check
 * usage: check the return value of a CUDA API call and print error msg
 *      if sth is wrong
 * arguments:
 *      1) err: the cudaError_t error
 *      2) file: filename of the source code
 *      3) line: line number of the source code
 *      4) fatal: a flag indicating if the program should abort on error
 * return: void
 */
__host__ void
cuda_check(cudaError_t err, const char* file, const int line, bool fatal);

/* a macro that fills in source file and line number */
#define CUDA_CHECK(call) do { \
    cuda_check((call), __FILE__, __LINE__, 1); \
} while(0)


#ifndef __cplusplus

typedef struct cudaDeviceProp cudaDeviceProp;

#endif

/* check and print info of a GPU device */
__host__ void
check_gpu(const cudaDeviceProp* const dev_prop);

#endif  /* __CUDA_UTIL_H__ */
