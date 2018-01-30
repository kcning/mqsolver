cmake_minimum_required(VERSION 3.8 FATAL_ERROR)  # CUDA support

project(mqsolver LANGUAGES C CUDA)

# set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

# ========================================================================== #
# language standard requirements
# ========================================================================== #

if(NOT DEFINED CMAKE_C_STANDARD)
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED true)
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED true)
endif()

# ========================================================================== #
# compilation flags
# ========================================================================== #

include(CheckCCompilerFlag)

check_c_compiler_flag("-Wall" COMPILER_C_ALL_WARN)
if(COMPILER_C_ALL_WARN)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall"
        CACHE STRING "Show all warning messages" FORCE
    )
endif()

check_c_compiler_flag("-Wextra" COMPILER_C_EXTRA_WARN)
if(COMPILER_C_EXTRA_WARN)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wextra"
        CACHE STRING "Show extra warning messages" FORCE
    )
endif()

check_c_compiler_flag("-pipe" COMPILER_C_PIPE)
if(COMPILER_C_PIPE)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pipe"
        CACHE STRING "Speed up compilation by piping" FORCE
    )
endif()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler='${CMAKE_C_FLAGS}'"
    CACHE STEING "Pass the same C flags to back-end compiler of nvcc" FORCE
)

check_c_compiler_flag("-mavx" COMPILER_C_AVX)
if(COMPILER_C_AVX)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx"
        CACHE STRING "Optimize with AVX instructions" FORCE
    )
endif()

check_c_compiler_flag("-mavx2" COMPILER_C_AVX2)
if(COMPILER_C_AVX2)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx2"
        CACHE STRING "Optimize with AVX2 instructions" FORCE
    )
endif()

check_c_compiler_flag("-march=native" COMPILER_C_ARCH)
if(COMPILER_C_ARCH)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native"
        CACHE STRING "Optimize based on the architecture" FORCE
    )
endif()

check_c_compiler_flag("-mtune=native" COMPILER_C_TUNE)
if(COMPILER_C_TUNE)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mtune=native"
        CACHE STRING "Fine-tune the program based on the CPU" FORCE
    )
endif()

check_c_compiler_flag("-fomit-frame-pointer" COMPILER_C_NO_FP)
if(COMPILER_C_NO_FP)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fomit-frame-pointer"
        CACHE STRING "Omit frame pointers" FORCE
    )
endif()

check_c_compiler_flag("-ffast-math" COMPILER_C_FAST_MATH)
if(COMPILER_C_FAST_MATH)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math"
        CACHE STRING "Trade floating point precision for speed" FORCE
    )
endif()

set(CMAKE_CUDA_FLAGS_RELEASE
    "${CMAKE_CUDA_FLAGS_RELEASE} -Xcompiler='${CMAKE_C_FLAGS_RELEASE}'"
    CACHE STEING "Pass the same C release flags to back-end compiler of nvcc"
    FORCE
)

# with an nvcc that supports CUDA standard 11, this option should be supported
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas '-dlcm=ca'"
    CACHE STRING "Enable L1 cache on GPU" FORCE
)

# ========================================================================== #
# dynamically generate C code
# ========================================================================== #

set(GEN_SRCS
    gc_decl_lsys.def
    gc_copy_lsys.def
    gc_check_lsys.def
    gc_solve_lsys.def
    gc_extract_sol.def
    gc_dep_lsys.def
    gc_gauss.def
)

if(NOT DEFINED KEEP_VAR_NUM)
    message(FATAL_ERROR "The number of variables to keep is not defined")
endif()

add_custom_command(OUTPUT ${GEN_SRCS}
    COMMAND python3 ${PROJECT_SOURCE_DIR}/bin/meta.py ${KEEP_VAR_NUM}
    DEPENDS ${PROJECT_SOURCE_DIR}/bin/meta.py
    COMMENT "Generate C source code based on the choice of k"
)

add_custom_target(gen_code DEPENDS ${GEN_SRCS}
    COMMENT "Check if re-generation is required"
)

# ========================================================================== #
# source code
# ========================================================================== #

set(C_SRCS
    src/mqsolver/util.c
    src/mqsolver/options.c
    src/mqsolver/algor.c
    src/mqsolver/mq_math.c
    src/mqsolver/debug.c
    src/mqsolver/mqfix.c
    src/mqsolver/drow.c
    src/mqsolver/macaulay.c
    src/threadpool/threadpool.c
)

set(CUDA_SRCS
    src/mqsolver/cuda_util.cu
    src/mqsolver/mqsolver.cu
    src/mqsolver/graycode.cu
    src/mqsolver/fix.cu
    src/mqsolver/rmac.cu
)

# ========================================================================== #
# main program
# ========================================================================== #

add_executable(${PROJECT_NAME}
    ${C_SRCS}
    ${CUDA_SRCS}
    src/main.c
)

target_include_directories(${PROJECT_NAME}
    PRIVATE ${PROJECT_SOURCE_DIR}/include/mqsolver
    PRIVATE ${PROJECT_SOURCE_DIR}/include/threadpool
    PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    PRIVATE ${PROJECT_BINARY_DIR}  # include generated C code
)

target_link_libraries(${PROJECT_NAME}
    PRIVATE m
    PRIVATE pthread
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    LINKER_LANGUAGE C
    CUDA_SEPARABLE_COMPILATION ON
)

target_compile_definitions(${PROJECT_NAME} PRIVATE
    KNUM=${KEEP_VAR_NUM}
)

add_dependencies(${PROJECT_NAME} gen_code)
