/*
 * SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 * SPDX-License-Identifier: LicenseRef-NvidiaProprietary
 *
 * NVIDIA CORPORATION, its affiliates and licensors retain all intellectual
 * property and proprietary rights in and to this material, related
 * documentation and any modifications thereto. Any use, reproduction,
 * disclosure or distribution of this material and related documentation
 * without an express license agreement from NVIDIA CORPORATION or
 * its affiliates is strictly prohibited.
 */

#ifndef CUDSS_DATA_TYPES_PUBLIC_HEADER_H
#define CUDSS_DATA_TYPES_PUBLIC_HEADER_H

#include <stddef.h>        // size_t
#include <library_types.h> // for cudaDataType_t
#include <cuda_runtime.h>  // for cudaStream_t

#ifdef __cplusplus
extern "C" {
#endif

struct cudssContext;
typedef struct cudssContext *cudssHandle_t; // the library context handle

struct cudssMatrix;
typedef struct cudssMatrix* cudssMatrix_t; // opaque generic matrix struct (for both dense and sparse matrices)

struct cudssData;
typedef struct cudssData* cudssData_t; // opaque object type which stores internal data (like LU factors) as well as some user-provided pointers

struct cudssConfig;
typedef struct cudssConfig* cudssConfig_t; // opaque object type which stores solver settings (e.g., algorithmic knobs)

typedef enum cudssDataType_t {
    CUDSS_DATA_TYPE_UNSET = 1024,
    CUDSS_R_32F = CUDA_R_32F,
    CUDSS_R_64F = CUDA_R_64F,
    CUDSS_C_32F = CUDA_C_32F,
    CUDSS_C_64F = CUDA_C_64F,
    CUDSS_R_64F_64F = 1025 + CUDA_R_64F,
    CUDSS_R_32I = CUDA_R_32I,
    CUDSS_R_64I = CUDA_R_64I
} cudssDataType_t;

typedef enum cudssConfigParam_t {
    CUDSS_CONFIG_REORDERING_ALG,
    CUDSS_CONFIG_FACTORIZATION_ALG,
    CUDSS_CONFIG_SOLVE_ALG,
    CUDSS_CONFIG_MATCHING_ALG,      // not used by default
    CUDSS_CONFIG_SOLVE_MODE,        // not supported right now (for transpose, conj transpose solves) (only off)
    CUDSS_CONFIG_IR_N_STEPS,
    CUDSS_CONFIG_IR_TOL,
    CUDSS_CONFIG_PIVOT_TYPE,
    CUDSS_CONFIG_PIVOT_THRESHOLD,
    CUDSS_CONFIG_PIVOT_EPSILON,
    CUDSS_CONFIG_MAX_LU_NNZ,              // only for CUDSS_REORDERING_ALG_BTF_COLAMD and CUDSS_REORDERING_ALG_COLAMD reordering algorithms
    CUDSS_CONFIG_HYBRID_MEMORY_MODE,             // by default: disabled
    CUDSS_CONFIG_HYBRID_DEVICE_MEMORY_LIMIT,
    CUDSS_CONFIG_USE_CUDA_REGISTER_MEMORY,// by default: enabled
    CUDSS_CONFIG_HOST_NTHREADS,
    CUDSS_CONFIG_HYBRID_EXECUTE_MODE,     // default: 0 - disabled
    CUDSS_CONFIG_PIVOT_EPSILON_ALG,
    CUDSS_CONFIG_ND_NLEVELS,
    CUDSS_CONFIG_UBATCH_SIZE, // "U" - stands for Uniform
    CUDSS_CONFIG_UBATCH_INDEX,
    CUDSS_CONFIG_USE_SUPERPANELS, // by default: enabled
    CUDSS_CONFIG_DEVICE_COUNT,
    CUDSS_CONFIG_DEVICE_INDICES,
    CUDSS_CONFIG_SCHUR_MODE,
    CUDSS_CONFIG_DETERMINISTIC_MODE,
    CUDSS_CONFIG_ND_UBFACTOR  //unbalance factor (in percent) for nested dissection
} cudssConfigParam_t;

typedef enum cudssDataParam_t {
    CUDSS_DATA_INFO,                 // (out)
    CUDSS_DATA_LU_NNZ,               // (out)
    CUDSS_DATA_NPIVOTS,              // (out)
    CUDSS_DATA_INERTIA,              // (out, non-trivial for non-positive-definite matrices)
    CUDSS_DATA_PERM_REORDER_ROW,     // (out)
    CUDSS_DATA_PERM_REORDER_COL,     // (out)
    CUDSS_DATA_PERM_ROW,             // (out, supported only for CUDSS_REORDERING_ALG_BTF_COLAMD and CUDSS_REORDERING_ALG_COLAMD reordering algorithms)
    CUDSS_DATA_PERM_COL,             // (out, supported only for CUDSS_REORDERING_ALG_BTF_COLAMD and CUDSS_REORDERING_ALG_COLAMD reordering algorithms)
    CUDSS_DATA_DIAG,                 // (out)
    CUDSS_DATA_USER_PERM,            // (in, out) for the user to provide a permutation or retrieve the provided permutation
    CUDSS_DATA_HYBRID_DEVICE_MEMORY_MIN,
    CUDSS_DATA_COMM_DEVICE,          // (in) communicator for device buffers
    CUDSS_DATA_COMM_HOST,            // (in) communicator for host buffers
    CUDSS_DATA_MEMORY_ESTIMATES,
    CUDSS_DATA_PERM_MATCHING,        // (out, supported only when matching is enabled)
    CUDSS_DATA_SCALE_ROW,            // (out, supported only when matching is enabled)
    CUDSS_DATA_SCALE_COL,            // (out, supported only when matching with scaling is enabled)
    CUDSS_DATA_NSUPERPANELS,         // (out)
    CUDSS_DATA_USER_SCHUR_INDICES,   // (in) indices[n] (1s for Schur complement, 0s for the rest)
    CUDSS_DATA_SCHUR_SHAPE,          // (out) int64_t shape[3]
    CUDSS_DATA_SCHUR_MATRIX,         // (out) cudssMatrix_t
    // (in, out) for the user to provide or retrieve the auxiliary partition tree information
    CUDSS_DATA_USER_ND_PARTITION_TREE,
    CUDSS_DATA_ND_PARTITION_TREE,     // (out) retrieves the partition tree information
    CUDSS_DATA_USER_HOST_INTERRUPT,
    CUDSS_DATA_IR_N_STEPS,           // (out) number of iterative refinement steps performed
    CUDSS_DATA_UBATCH_MASK,          // (in, out) integer array (0/1 values) to selectively process matrices in uniform batch
    CUDSS_DATA_FLOPS                 // (out)
} cudssDataParam_t;

typedef enum cudssPhase_t {
    CUDSS_PHASE_REORDERING             = 1 << 0,
    CUDSS_PHASE_SYMBOLIC_FACTORIZATION = 1 << 1,
    CUDSS_PHASE_ANALYSIS               = CUDSS_PHASE_REORDERING | CUDSS_PHASE_SYMBOLIC_FACTORIZATION,
    CUDSS_PHASE_FACTORIZATION          = 1 << 2,
    CUDSS_PHASE_REFACTORIZATION        = 1 << 3,
    CUDSS_PHASE_SOLVE_FWD_PERM         = 1 << 4,
    CUDSS_PHASE_SOLVE_FWD              = 1 << 5,
    CUDSS_PHASE_SOLVE_DIAG             = 1 << 6,
    CUDSS_PHASE_SOLVE_BWD              = 1 << 7,
    CUDSS_PHASE_SOLVE_BWD_PERM         = 1 << 8,
    CUDSS_PHASE_SOLVE_REFINEMENT       = 1 << 9,
    CUDSS_PHASE_SOLVE                  = CUDSS_PHASE_SOLVE_FWD_PERM | CUDSS_PHASE_SOLVE_FWD | CUDSS_PHASE_SOLVE_DIAG | CUDSS_PHASE_SOLVE_BWD | CUDSS_PHASE_SOLVE_BWD_PERM | CUDSS_PHASE_SOLVE_REFINEMENT
} cudssPhase_t;

typedef enum cudssStatus_t {
    CUDSS_STATUS_SUCCESS,
    CUDSS_STATUS_NOT_INITIALIZED,
    CUDSS_STATUS_ALLOC_FAILED,
    CUDSS_STATUS_INVALID_VALUE,
    CUDSS_STATUS_NOT_SUPPORTED,
    CUDSS_STATUS_EXECUTION_FAILED,
    CUDSS_STATUS_INTERNAL_ERROR,
    CUDSS_STATUS_IR_FAILED
} cudssStatus_t;

typedef enum cudssMatrixType_t {
    CUDSS_MTYPE_GENERAL,
    CUDSS_MTYPE_SYMMETRIC,
    CUDSS_MTYPE_HERMITIAN,
    CUDSS_MTYPE_SPD,
    CUDSS_MTYPE_HPD
} cudssMatrixType_t;

typedef enum cudssMatrixViewType_t {
    CUDSS_MVIEW_FULL,
    CUDSS_MVIEW_LOWER,
    CUDSS_MVIEW_UPPER
} cudssMatrixViewType_t;

typedef enum cudssIndexBase_t {
    CUDSS_BASE_ZERO,
    CUDSS_BASE_ONE
} cudssIndexBase_t;

typedef enum cudssLayout_t {
    CUDSS_LAYOUT_COL_MAJOR,
    CUDSS_LAYOUT_ROW_MAJOR
} cudssLayout_t;

typedef enum cudssReorderingAlg_t {
    CUDSS_REORDERING_ALG_DEFAULT,
    CUDSS_REORDERING_ALG_BTF_COLAMD,
    CUDSS_REORDERING_ALG_COLAMD,
    CUDSS_REORDERING_ALG_AMD,
    CUDSS_REORDERING_ALG_NESTED_DISSECTION,
    CUDSS_REORDERING_ALG_NONE   // Uses natural (identity) order for the internal ordering when no user permutation is supplied.
} cudssReorderingAlg_t;

typedef enum cudssFactorizationAlg_t {
    CUDSS_FACTORIZATION_ALG_DEFAULT,
    CUDSS_FACTORIZATION_ALG_MULTIBLOCK,
    CUDSS_FACTORIZATION_ALG_GENERAL
} cudssFactorizationAlg_t;

typedef enum cudssPivotEpsilonAlg_t {
    CUDSS_PIVOT_EPSILON_ALG_DEFAULT,
    CUDSS_PIVOT_EPSILON_ALG_SCALED,
    CUDSS_PIVOT_EPSILON_ALG_STATIC,
} cudssPivotEpsilonAlg_t;

typedef enum cudssSolveAlg_t {
    CUDSS_SOLVE_ALG_DEFAULT,
    CUDSS_SOLVE_ALG_GENERAL
} cudssSolveAlg_t;

typedef enum cudssMatchingAlg_t {
    CUDSS_MATCHING_ALG_NONE = 0,        /// Matching is disabled
    CUDSS_MATCHING_ALG_MAX_DIAG_COUNT = 1,
    CUDSS_MATCHING_ALG_MAX_MIN_DIAG = 2,
    CUDSS_MATCHING_ALG_MAX_MIN_DIAG_ALT = 3,
    CUDSS_MATCHING_ALG_MAX_DIAG_SUM = 4,
    CUDSS_MATCHING_ALG_MAX_DIAG_PRODUCT = 5,
    CUDSS_MATCHING_ALG_AUTO = 6,        /// The library selects the alg to apply
} cudssMatchingAlg_t;

typedef enum cudssPivotType_t {
    CUDSS_PIVOT_AUTO,          // Automatic pivoting strategy [default] (resolves based on algorithm and matrix type)
    CUDSS_PIVOT_NONE,          // Disable pivoting search (valid for all algorithms) 
    CUDSS_PIVOT_GLOBAL_COL,    // For CUDSS_REORDERING_ALG_BTF_COLAMD and CUDSS_REORDERING_ALG_COLAMD - all matrix types
    CUDSS_PIVOT_GLOBAL_ROW,    // For CUDSS_REORDERING_ALG_BTF_COLAMD and CUDSS_REORDERING_ALG_COLAMD - all matrix types
    CUDSS_PIVOT_DIAGONAL,      // For symmetric indefinite matrices, CUDSS_REORDERING_ALG_AMD and CUDSS_REORDERING_ALG_NESTED_DISSECTION
    CUDSS_PIVOT_LOCAL_BLOCK,   // For general matrices, CUDSS_REORDERING_ALG_AMD and CUDSS_REORDERING_ALG_NESTED_DISSECTION
    CUDSS_PIVOT_BUNCH_KAUFMAN  // Reserved for future (not supported yet)
} cudssPivotType_t;

typedef enum cudssMatrixFormat_t {
  CUDSS_MFORMAT_DENSE       = 1,
  CUDSS_MFORMAT_CSR         = 2,
  CUDSS_MFORMAT_BATCH       = 4,
  CUDSS_MFORMAT_DISTRIBUTED = 8
} cudssMatrixFormat_t;

#define CUDSS_ALLOCATOR_NAME_LEN 64

typedef struct {
  void* ctx;
  int (*device_alloc)(void* ctx, void** ptr, size_t size, cudaStream_t stream);
  int (*device_free)(void* ctx, void* ptr, size_t size, cudaStream_t stream);
  char name[CUDSS_ALLOCATOR_NAME_LEN];
} cudssDeviceMemHandler_t;


#ifdef __cplusplus

// Type corresponding to CUDSS_R_64F_64F
struct alignas(2 * sizeof(double)) cudss_fp64mp2_t {
    double hi;
    double lo;
};

#else

// Type corresponding to CUDSS_R_64F_64F
typedef struct {
    _Alignas(2 * sizeof(double)) double hi;
    double lo;
} cudss_fp64mp2_t;

#endif // __cplusplus


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CUDSS_DATA_TYPES_PUBLIC_HEADER_H */
