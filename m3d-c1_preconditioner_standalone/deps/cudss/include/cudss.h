/*
 * Copyright 2023-2026 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO LICENSEE:
 *
 * This source code and/or documentation ("Licensed Deliverables") are
 * subject to NVIDIA intellectual property rights under U.S. and
 * international Copyright laws.
 *
 * These Licensed Deliverables contained herein is PROPRIETARY and
 * CONFIDENTIAL to NVIDIA and is being provided under the terms and
 * conditions of a form of NVIDIA software license agreement by and
 * between NVIDIA and Licensee ("License Agreement") or electronically
 * accepted by Licensee.  Notwithstanding any terms or conditions to
 * the contrary in the License Agreement, reproduction or disclosure
 * of the Licensed Deliverables to any third party without the express
 * written consent of NVIDIA is prohibited.
 *
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, NVIDIA MAKES NO REPRESENTATION ABOUT THE
 * SUITABILITY OF THESE LICENSED DELIVERABLES FOR ANY PURPOSE.  IT IS
 * PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND.
 * NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THESE LICENSED
 * DELIVERABLES, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
 * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY
 * SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY
 * DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THESE LICENSED DELIVERABLES.
 *
 * U.S. Government End Users.  These Licensed Deliverables are a
 * "commercial item" as that term is defined at 48 C.F.R. 2.101 (OCT
 * 1995), consisting of "commercial computer software" and "commercial
 * computer software documentation" as such terms are used in 48
 * C.F.R. 12.212 (SEPT 1995) and is provided to the U.S. Government
 * only as a commercial end item.  Consistent with 48 C.F.R.12.212 and
 * 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the Licensed Deliverables with
 * only those rights set forth herein.
 *
 * Any use of the Licensed Deliverables in individual and commercial
 * software must include, in the user documentation and internal
 * comments to the code, the above Disclaimer and U.S. Government End
 * Users Notice.
 */

#ifndef CUDSS_PUBLIC_HEADER_H
#define CUDSS_PUBLIC_HEADER_H

#include <cuda_runtime.h>  // for cudaStream_t
#include <stddef.h>        // size_t
#include <stdint.h>        // int64_t
#include <stdio.h>         // FILE

#include "cudss_data_types.h"
#include "cudss_distributed_interface.h"
#include "cudss_threading_interface.h"

#define CUDSS_VERSION_MAJOR 0
#define CUDSS_VERSION_MINOR 8
#define CUDSS_VERSION_PATCH 0
#define CUDSS_VERSION (CUDSS_VERSION_MAJOR * 10000 + \
                       CUDSS_VERSION_MINOR *  100 +  \
                       CUDSS_VERSION_PATCH)

#if !defined(CUDSSAPI)
#    if defined(_WIN32)
#        define CUDSSAPI __stdcall
#    else
#        define CUDSSAPI
#    endif
#endif

#ifdef __cplusplus
extern "C" {
#endif


// Set/Get APIs for cudssConfig_t and cudssData_t

cudssStatus_t CUDSSAPI cudssConfigSet(cudssConfig_t config, cudssConfigParam_t param, const void *value, size_t sizeInBytes);

cudssStatus_t CUDSSAPI cudssConfigGet(const cudssConfig_t config, cudssConfigParam_t param,  void *value, size_t sizeInBytes, size_t *sizeWritten);

cudssStatus_t CUDSSAPI cudssDataSet(const cudssHandle_t handle, cudssData_t data, cudssDataParam_t param, const void *value, size_t sizeInBytes);

cudssStatus_t CUDSSAPI cudssDataGet(const cudssHandle_t handle, const cudssData_t data, cudssDataParam_t param, void *value, size_t sizeInBytes, size_t *sizeWritten);

// Main cuDSS routine

cudssStatus_t CUDSSAPI cudssExecute(cudssHandle_t handle, int phase, const cudssConfig_t solverConfig, cudssData_t solverData, const cudssMatrix_t inputMatrix, cudssMatrix_t solution, const cudssMatrix_t rhs);

// Setting the stream (in the library handle)

cudssStatus_t CUDSSAPI cudssSetStream(cudssHandle_t handle, cudaStream_t stream);

// Setting per-device streams for multi-GPU mode (in the library handle)
// streams[i] will be used for device_indices[i]
cudssStatus_t CUDSSAPI cudssSetMgStreams(cudssHandle_t handle, const cudaStream_t* streams, int stream_count);

// Setting the communication layer library name (in the library handle)

cudssStatus_t CUDSSAPI cudssSetCommLayer(cudssHandle_t handle, const char* commLibFileName);

// Setting the threading layer library name (in the library handle)

cudssStatus_t CUDSSAPI cudssSetThreadingLayer(cudssHandle_t handle, const char* thrLibFileName);

// Create/Destroy APIs (allocating structures + set defaults)

cudssStatus_t CUDSSAPI cudssConfigCreate(cudssConfig_t *solverConfig);
cudssStatus_t CUDSSAPI cudssConfigDestroy(cudssConfig_t solverConfig);

cudssStatus_t CUDSSAPI cudssDataCreate(const cudssHandle_t handle, cudssData_t *solverData);
cudssStatus_t CUDSSAPI cudssDataDestroy(cudssHandle_t handle, cudssData_t solverData);

cudssStatus_t CUDSSAPI cudssCreate(cudssHandle_t *handle);
cudssStatus_t CUDSSAPI cudssCreateMg(cudssHandle_t *handle_pt, int device_count, const int *device_indices);
cudssStatus_t CUDSSAPI cudssDestroy(cudssHandle_t handle);

// Versioning

cudssStatus_t CUDSSAPI cudssGetProperty(libraryPropertyType propertyType, int* value);

// Create/Destroy API helpers for matrix wrappers

cudssStatus_t CUDSSAPI cudssMatrixCreateDn(cudssMatrix_t *matrix, int64_t nrows, int64_t ncols, int64_t ld, const void *values, cudssDataType_t valueType,  cudssLayout_t layout);

cudssStatus_t CUDSSAPI cudssMatrixCreateCsr(cudssMatrix_t *matrix, int64_t nrows, int64_t ncols, int64_t nnz, const void *rowStart, const void *rowEnd, const void *colIndices, const void *values, cudssDataType_t offsetType, cudssDataType_t indexType, cudssDataType_t valueType, cudssMatrixType_t mtype, cudssMatrixViewType_t mview, cudssIndexBase_t indexBase);

cudssStatus_t CUDSSAPI cudssMatrixCreateBatchDn(cudssMatrix_t *matrix, int64_t batchCount, const void *nrows, const void *ncols, const void *ld, const void *const *values, cudssDataType_t integerType, cudssDataType_t valueType,  cudssLayout_t layout);

cudssStatus_t CUDSSAPI cudssMatrixCreateBatchCsr(cudssMatrix_t *matrix, int64_t batchCount, const void *nrows, const void *ncols, const void *nnz, const void *const *rowStart, const void *const *rowEnd, const void *const *colIndices, const void *const *values, cudssDataType_t offsetType, cudssDataType_t indexType, cudssDataType_t valueType, cudssMatrixType_t mtype, cudssMatrixViewType_t mview, cudssIndexBase_t indexBase);

cudssStatus_t CUDSSAPI cudssMatrixDestroy(cudssMatrix_t matrix);

// Setters/Getters API helpers for matrix wrappers

cudssStatus_t CUDSSAPI cudssMatrixGetDn(const cudssMatrix_t matrix,  int64_t* nrows, int64_t* ncols, int64_t* ld, void **values, cudssDataType_t* type, cudssLayout_t* layout);

cudssStatus_t CUDSSAPI cudssMatrixGetCsr(const cudssMatrix_t matrix, int64_t* nrows, int64_t* ncols, int64_t* nnz, void **rowStart, void **rowEnd, void **colIndices, void **values, cudssDataType_t* offsetType, cudssDataType_t* indexType, cudssDataType_t* valueType, cudssMatrixType_t* mtype, cudssMatrixViewType_t* mview, cudssIndexBase_t* indexBase);

cudssStatus_t CUDSSAPI cudssMatrixSetValues(cudssMatrix_t matrix, const void *values);

cudssStatus_t CUDSSAPI cudssMatrixSetCsrPointers(cudssMatrix_t matrix, const void *rowOffsets, const void *rowEnd, const void *colIndices, const void *values);

cudssStatus_t CUDSSAPI cudssMatrixGetBatchDn(const cudssMatrix_t matrix, int64_t *batchCount, void **nrows, void **ncols, void **ld, void ***values, cudssDataType_t *indexType, cudssDataType_t *valueType, cudssLayout_t* layout);

cudssStatus_t CUDSSAPI cudssMatrixGetBatchCsr(const cudssMatrix_t matrix, int64_t *batchCount, void **nrows, void **ncols, void **nnz, void ***rowStart, void ***rowEnd, void ***colIndices, void ***values, cudssDataType_t* offsetType, cudssDataType_t* indexType, cudssDataType_t* valueType, cudssMatrixType_t* mtype, cudssMatrixViewType_t* mview, cudssIndexBase_t* indexBase);

cudssStatus_t CUDSSAPI cudssMatrixSetBatchValues(cudssMatrix_t matrix, const void *const *values);

cudssStatus_t CUDSSAPI cudssMatrixSetBatchCsrPointers(cudssMatrix_t matrix, const void *const *rowOffsets, const void *const *rowEnd, const void *const *colIndices, const void *const *values);

cudssStatus_t CUDSSAPI cudssMatrixGetFormat(const cudssMatrix_t matrix, int* format);

cudssStatus_t CUDSSAPI cudssMatrixSetDistributionRow1d(cudssMatrix_t matrix, int64_t first_row, int64_t last_row);

cudssStatus_t CUDSSAPI cudssMatrixGetDistributionRow1d(const cudssMatrix_t matrix, int64_t *first_row, int64_t *last_row);

// Memory allocator API

cudssStatus_t CUDSSAPI cudssGetDeviceMemHandler(const cudssHandle_t handle, cudssDeviceMemHandler_t* handler);

cudssStatus_t CUDSSAPI cudssSetDeviceMemHandler(cudssHandle_t handle, const cudssDeviceMemHandler_t* handler);

// Logger API

typedef void(*cudssLoggerCallback_t)(int logLevel, const char* functionName, const char* message);

cudssStatus_t CUDSSAPI cudssLoggerSetCallback(cudssLoggerCallback_t callback);

cudssStatus_t CUDSSAPI cudssLoggerSetFile(FILE* file);

cudssStatus_t CUDSSAPI cudssLoggerOpenFile(const char* logFile);

cudssStatus_t CUDSSAPI cudssLoggerSetLevel(int level);

cudssStatus_t CUDSSAPI cudssLoggerSetMask(int mask);

cudssStatus_t CUDSSAPI cudssLoggerForceDisable(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CUDSS_PUBLIC_HEADER_H */
