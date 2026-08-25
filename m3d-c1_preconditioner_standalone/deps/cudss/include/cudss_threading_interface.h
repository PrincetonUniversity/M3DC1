/*
 * SPDX-FileCopyrightText: Copyright (c) 2025-2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
 * SPDX-License-Identifier: LicenseRef-NvidiaProprietary
 *
 * NVIDIA CORPORATION, its affiliates and licensors retain all intellectual
 * property and proprietary rights in and to this material, related
 * documentation and any modifications thereto. Any use, reproduction,
 * disclosure or distribution of this material and related documentation
 * without an express license agreement from NVIDIA CORPORATION or
 * its affiliates is strictly prohibited.
*/

#ifndef CUDSS_THR_IFACE_PUBLIC_HEADER_H
#define CUDSS_THR_IFACE_PUBLIC_HEADER_H

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*cudss_thr_func_t)(int task, void *ctx);

/* Interface to threading layer */
typedef struct {
    int (*cudssGetMaxThreads)();
    void (*cudssParallelFor)(int nthr_requested, int ntasks, void *ctx, cudss_thr_func_t f);
    int (*cudssThreadingGetProperty)(libraryPropertyType propertyType, int* value);
} cudssThreadingInterface_t;

#ifdef __cplusplus
} /* extern "C" */
#endif
#endif /* CUDSS_THR_IFACE_PUBLIC_HEADER_H */
