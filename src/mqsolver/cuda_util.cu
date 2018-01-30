/* cuda_util.cu: implementation of cuda_util.h */

extern "C" {

#include "cuda_util.h"

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
cuda_check(cudaError_t err, const char* file, const int line, bool fatal) {
    if(cudaSuccess != err) {
        PRINTF_ERR_STAMP("[!] GPU error: %s:%d, code: %d, reason: %s\n",
                         file, line, err, cudaGetErrorString(err));

        if(fatal) {
            PRINTF_ERR_STAMP("[!] aborting...\n");
            exit(err);
        }
    }
}

/* check and print info of a GPU device */
__host__ void
check_gpu(const cudaDeviceProp* const dev_prop) {

    PRINTF_STAMP("[+] CUDA info:\n");
    int cuda_rtv;
    CUDA_CHECK(cudaRuntimeGetVersion(&cuda_rtv));
    PRINTF_STAMP("\t\tRuntime version: %d.%d\n", cuda_rtv / 1000,
                 (cuda_rtv % 100) / 10);

    PRINTF_STAMP("\t\tCapability: %d.%d\n", dev_prop->major, dev_prop->minor);

    PRINTF_STAMP("\t\tGlobal memory size: %.2fMB (%" PRIu64 " bytes)\n",
           dev_prop->totalGlobalMem / MBFLOAT,
           (uint64_t) dev_prop->totalGlobalMem);

    PRINTF_STAMP("\t\tGPU clock rate: %0.2f MHz\n", dev_prop->clockRate / 1000.0f);
    PRINTF_STAMP("\t\tMemory clock rate: %0.2f MHz\n",
                 dev_prop->memoryClockRate / 1000.0f);
    PRINTF_STAMP("\t\tMemory bus width: %d bits\n", dev_prop->memoryBusWidth);
    PRINTF_STAMP("\t\tMax memory pitch: %0.2f MB\n",
                 dev_prop->memPitch / MBFLOAT);
    PRINTF_STAMP("\t\tL1 support global cache? %s\n",
                 dev_prop->globalL1CacheSupported ? "yes" : "no");
    PRINTF_STAMP("\t\tL1 support local cache? %s\n",
                 dev_prop->localL1CacheSupported ? "yes" : "no");

    PRINTF_STAMP("\t\tL2 cache size: %d bytes\n", dev_prop->l2CacheSize);
    PRINTF_STAMP("\t\tConstant memory size: %lu bytes\n", dev_prop->totalConstMem);
    PRINTF_STAMP("\t\tShared memory size per block: %lu bytes\n",
                 dev_prop->sharedMemPerBlock);

    PRINTF_STAMP("\t\tNumber of registers available per block: %d\n",
                 dev_prop->regsPerBlock);
    PRINTF_STAMP("\t\tMax number of threads per block: %d\n",
                 dev_prop->maxThreadsPerBlock);
    PRINTF_STAMP("\t\tNumber of registers available per thread: %d\n",
                 dev_prop->regsPerBlock / dev_prop->maxThreadsPerBlock);

    PRINTF_STAMP("\t\tWarp size: %d\n", dev_prop->warpSize);
    PRINTF_STAMP("\t\tMax number of threads per multiprocessor: %d\n",
                 dev_prop->maxThreadsPerMultiProcessor);
    PRINTF_STAMP("\t\tNumber of multiprocessors: %d\n", dev_prop->multiProcessorCount);

    PRINTF_STAMP("\t\tMax sizes of each dimension of a block: (%d x %d x %d)\n",
                 dev_prop->maxThreadsDim[0],
                 dev_prop->maxThreadsDim[1],
                 dev_prop->maxThreadsDim[2]);
    PRINTF_STAMP("\t\tMax sizes of each dimension of a grid: (%d x %d x %d)\n",
                 dev_prop->maxGridSize[0],
                 dev_prop->maxGridSize[1],
                 dev_prop->maxGridSize[2]);
    PRINTF_STAMP("\t\tConcurrent copy and execution? %s\n",
                 dev_prop->deviceOverlap ? "yes" : "no");
    PRINTF_STAMP("\t\tLaunch concurrent kernels? %s\n",
                 dev_prop->concurrentKernels ? "yes" : "no");
    //PRINTF_STAMP("\t\tSingle/Double precision performance ratio: %d\n",
    //             dev_prop->singleToDoublePrecisionPerfRatio);
    PRINTF_STAMP("\t\tNumber of asynchronous engines: %d\n",
                 dev_prop->asyncEngineCount);
    //PRINTF_STAMP("\t\tNative atomic operations between host and device?: %s\n",
    //             dev_prop->hostNativeAtomicSupported ? "yes" : "no");
    
    PRINTF_STAMP("[+] end CUDA info\n");
}

}
