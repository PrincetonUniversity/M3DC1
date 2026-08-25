// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * cudss_commlayer_craympi.c -- custom cuDSS MGMN communication layer on
 * cray-mpich (GPU-aware MPI), replacing the prebuilt NCCL / OpenMPI shims.
 *
 * Ported to cuDSS 0.8.0: cudssDistributedInterface_t now has 21 function
 * pointers -- each of the 10 original collectives/p2p split into a `…Device`
 * and a `…Host` variant, plus `cudssDistributedGetProperty`. All datatype
 * arguments are now `cudssDataType_t` (was `cudaDataType_t`).
 *
 * cuDSS loads this .so via cudssSetCommLayer(handle, "<path>"); it resolves the
 * global symbol `cudssDistributedInterface` and calls through it. The
 * `void *comm` passed to every function is exactly what the driver supplied via
 *   cudssDataSet(handle, data, CUDSS_DATA_COMM_DEVICE/_HOST, &comm, sizeof(comm));
 * For this shim the driver stores an MPI_Comm in a stable location and passes
 * its address, so every `comm` argument here is an `MPI_Comm *`.
 *
 * GPU-awareness: the Device variants operate on DEVICE pointers. cray-mpich is
 * GPU-aware (MPICH_GPU_SUPPORT_ENABLED=1) but MPI is stream-unaware, so the
 * Device variants cudaStreamSynchronize(stream) at the start of every
 * collective/p2p to ensure (a) the GPU has finished producing the send buffer
 * before MPI reads it and (b) results are visible once the (synchronous) MPI
 * call returns. The Host variants operate on HOST pointers and need no stream
 * sync.
 *
 * Build: cc -O2 -fPIC -shared -I$CUDSS_INC cudss_commlayer_craympi.c -o \
 *        libcudss_commlayer_craympi.so
 * The Cray `cc` wrapper provides GPU-aware cray-mpich + MPI includes.
 */
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>
#include <cuda_runtime.h>
#include <library_types.h>            /* libraryPropertyType */

#include "cudss_distributed_interface.h"
#include "cudss.h"                    /* CUDSS_VERSION_* for GetProperty */

/* All interface functions return 0 on success (matching the NCCL/OpenMPI
 * layers: a nonzero return signals a comm-layer error to cuDSS). */
#define CUDSS_COMM_OK   0
#define CUDSS_COMM_ERR  1

/* Optional call trace (CUDSS_SHIM_TRACE=1): prints each comm call cuDSS makes,
 * so we can see the SOLVE-phase collective sequence and pin down data bugs. */
static int trace_on(void)
{
    static int t = -1;
    if (t < 0) { const char *e = getenv("CUDSS_SHIM_TRACE"); t = (e && atoi(e)) ? 1 : 0; }
    return t;
}
#define TRACE(fmt, ...) do { if (trace_on()) { \
    fprintf(stderr, "[shim] " fmt "\n", ##__VA_ARGS__); fflush(stderr); } } while (0)

/* Synchronize the CUDA stream so the device has finished producing the buffer
 * before stream-unaware MPI touches it. stream may be 0 (the default stream). */
static int sync_stream(cudaStream_t stream)
{
    cudaError_t e = cudaStreamSynchronize(stream);
    if (e != cudaSuccess) {
        fprintf(stderr,
                "[cudss_commlayer_craympi] cudaStreamSynchronize failed: %s\n",
                cudaGetErrorString(e));
        return CUDSS_COMM_ERR;
    }
    return CUDSS_COMM_OK;
}

/* Map cuDSS's cudssDataType_t to an MPI_Datatype. Abort with a clear message on
 * an unmapped type rather than silently corrupting a collective. */
static MPI_Datatype map_dtype(cudssDataType_t dt)
{
    switch (dt) {
        case CUDSS_R_32I: return MPI_INT;        /* 32-bit signed int   */
        case CUDSS_R_64F: return MPI_DOUBLE;     /* 64-bit real         */
        case CUDSS_R_32F: return MPI_FLOAT;      /* 32-bit real         */
        case CUDSS_R_64I: return MPI_LONG_LONG;  /* 64-bit signed int   */
        default:
            fprintf(stderr,
                    "[cudss_commlayer_craympi] FATAL: unmapped cudssDataType_t "
                    "%d -- add it to map_dtype().\n", (int)dt);
            MPI_Abort(MPI_COMM_WORLD, 2);
            return MPI_DATATYPE_NULL; /* unreachable */
    }
}

/* Map cuDSS's cudssOpType_t to an MPI_Op. */
static MPI_Op map_op(cudssOpType_t op)
{
    switch (op) {
        case CUDSS_SUM: return MPI_SUM;
        case CUDSS_MAX: return MPI_MAX;
        case CUDSS_MIN: return MPI_MIN;
        default:
            fprintf(stderr,
                    "[cudss_commlayer_craympi] FATAL: unmapped cudssOpType_t "
                    "%d -- add it to map_op().\n", (int)op);
            MPI_Abort(MPI_COMM_WORLD, 3);
            return MPI_OP_NULL; /* unreachable */
    }
}

/* ===================================================================== *
 *  Core MPI implementations (datatype/op-mapped). `do_sync` selects the
 *  Device (GPU-aware, stream-synced) vs Host (no sync) behavior; the MPI
 *  calls are otherwise identical on device- or host-resident buffers.
 * ===================================================================== */

static int impl_rank(void *comm, int *rank)
{
    return (MPI_Comm_rank(*(MPI_Comm *)comm, rank) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_size(void *comm, int *size)
{
    return (MPI_Comm_size(*(MPI_Comm *)comm, size) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_send(const void *buffer, int count, cudssDataType_t datatype,
                     int dest, int tag, void *comm, cudaStream_t stream,
                     int do_sync)
{
    TRACE("Send  count=%d dt=%d dest=%d tag=%d host=%d stream=%p", count, (int)datatype, dest, tag, !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    return (MPI_Send(buffer, count, map_dtype(datatype), dest, tag,
                     *(MPI_Comm *)comm) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_recv(void *buffer, int count, cudssDataType_t datatype,
                     int root, int tag, void *comm, cudaStream_t stream,
                     int do_sync)
{
    TRACE("Recv  count=%d dt=%d src=%d tag=%d host=%d stream=%p", count, (int)datatype, root, tag, !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    return (MPI_Recv(buffer, count, map_dtype(datatype), root, tag,
                     *(MPI_Comm *)comm, MPI_STATUS_IGNORE) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_bcast(void *buffer, int count, cudssDataType_t datatype,
                      int root, void *comm, cudaStream_t stream, int do_sync)
{
    TRACE("Bcast count=%d dt=%d root=%d host=%d stream=%p", count, (int)datatype, root, !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    return (MPI_Bcast(buffer, count, map_dtype(datatype), root,
                      *(MPI_Comm *)comm) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_reduce(const void *sendbuf, void *recvbuf, int count,
                       cudssDataType_t datatype, cudssOpType_t op, int root,
                       void *comm, cudaStream_t stream, int do_sync)
{
    TRACE("Reduce count=%d dt=%d op=%d root=%d inplace=%d host=%d stream=%p", count, (int)datatype, (int)op, root, (int)(sendbuf==recvbuf), !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    MPI_Comm c = *(MPI_Comm *)comm;
    /* cuDSS (like NCCL, which allows aliased reduce) may pass sendbuf==recvbuf
     * for an in-place reduction. MPI forbids aliasing and requires MPI_IN_PLACE,
     * which for MPI_Reduce is significant ONLY at the root (the root takes its
     * input from recvbuf and overwrites it with the result). cray-mpich's alias
     * check fires only at the root, so non-root ranks pass the buffers as-is. */
    const void *sb = sendbuf;
    if (sendbuf == recvbuf) {
        int myrank;
        if (MPI_Comm_rank(c, &myrank) != MPI_SUCCESS) return CUDSS_COMM_ERR;
        if (myrank == root) sb = MPI_IN_PLACE;
    }
    return (MPI_Reduce(sb, recvbuf, count, map_dtype(datatype),
                       map_op(op), root, c) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_allreduce(const void *sendbuf, void *recvbuf, int count,
                          cudssDataType_t datatype, cudssOpType_t op,
                          void *comm, cudaStream_t stream, int do_sync)
{
    TRACE("Allreduce count=%d dt=%d op=%d inplace=%d host=%d stream=%p", count, (int)datatype, (int)op, (int)(sendbuf==recvbuf), !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    /* cuDSS may pass sendbuf == recvbuf for an in-place reduction; MPI requires
     * MPI_IN_PLACE in that case (passing the same pointer is undefined). */
    const void *sb = (sendbuf == recvbuf) ? MPI_IN_PLACE : sendbuf;
    return (MPI_Allreduce(sb, recvbuf, count, map_dtype(datatype),
                          map_op(op), *(MPI_Comm *)comm) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

static int impl_scatterv(const void *sendbuf, const int *sendcounts,
                         const int *displs, cudssDataType_t sendtype,
                         void *recvbuf, int recvcount, cudssDataType_t recvtype,
                         int root, void *comm, cudaStream_t stream, int do_sync)
{
    TRACE("Scatterv recvcount=%d sdt=%d rdt=%d root=%d host=%d stream=%p", recvcount, (int)sendtype, (int)recvtype, root, !do_sync, (void *)stream);
    if (do_sync && sync_stream(stream)) return CUDSS_COMM_ERR;
    return (MPI_Scatterv(sendbuf, sendcounts, displs, map_dtype(sendtype),
                         recvbuf, recvcount, map_dtype(recvtype), root,
                         *(MPI_Comm *)comm) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

/* cuDSS allocated the storage that `newcomm` points at, sizing it from the
 * sizeof passed at CUDSS_DATA_COMM_* (sizeof(MPI_Comm*) = 8 in the driver).
 * Under cray-mpich (MPICH) MPI_Comm is a 4-byte int, so writing it into that
 * 8-byte slot is safe. */
static int impl_split(const void *comm, int color, int key, void *newcomm)
{
    TRACE("CommSplit color=%d key=%d", color, key);
    return (MPI_Comm_split(*(MPI_Comm *)comm, color, key,
                           (MPI_Comm *)newcomm) == MPI_SUCCESS)
               ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

/* Free a communicator previously produced by comm_split. Guard against freeing
 * predefined / null communicators (MPI forbids freeing those). */
static int impl_free(void *comm)
{
    MPI_Comm *c = (MPI_Comm *)comm;
    if (c == NULL || *c == MPI_COMM_NULL ||
        *c == MPI_COMM_WORLD || *c == MPI_COMM_SELF) {
        return CUDSS_COMM_OK;
    }
    return (MPI_Comm_free(c) == MPI_SUCCESS) ? CUDSS_COMM_OK : CUDSS_COMM_ERR;
}

/* ===================================================================== *
 *  Device variants: GPU-resident buffers -> stream-synced GPU-aware MPI.
 * ===================================================================== */

static int comm_rank_dev(void *comm, int *rank) { return impl_rank(comm, rank); }
static int comm_size_dev(void *comm, int *size) { return impl_size(comm, size); }

static int comm_send_dev(const void *buffer, int count, cudssDataType_t dt,
                         int dest, int tag, void *comm, cudaStream_t stream)
{ return impl_send(buffer, count, dt, dest, tag, comm, stream, 1); }

static int comm_recv_dev(void *buffer, int count, cudssDataType_t dt,
                         int root, int tag, void *comm, cudaStream_t stream)
{ return impl_recv(buffer, count, dt, root, tag, comm, stream, 1); }

static int comm_bcast_dev(void *buffer, int count, cudssDataType_t dt,
                          int root, void *comm, cudaStream_t stream)
{ return impl_bcast(buffer, count, dt, root, comm, stream, 1); }

static int comm_reduce_dev(const void *sb, void *rb, int count,
                           cudssDataType_t dt, cudssOpType_t op, int root,
                           void *comm, cudaStream_t stream)
{ return impl_reduce(sb, rb, count, dt, op, root, comm, stream, 1); }

static int comm_allreduce_dev(const void *sb, void *rb, int count,
                              cudssDataType_t dt, cudssOpType_t op,
                              void *comm, cudaStream_t stream)
{ return impl_allreduce(sb, rb, count, dt, op, comm, stream, 1); }

static int comm_scatterv_dev(const void *sb, const int *sc, const int *displs,
                             cudssDataType_t st, void *rb, int rc,
                             cudssDataType_t rt, int root, void *comm,
                             cudaStream_t stream)
{ return impl_scatterv(sb, sc, displs, st, rb, rc, rt, root, comm, stream, 1); }

static int comm_split_dev(const void *comm, int color, int key, void *newcomm)
{ return impl_split(comm, color, key, newcomm); }

static int comm_free_dev(void *comm) { return impl_free(comm); }

/* ===================================================================== *
 *  Host variants: host-resident buffers -> plain MPI, no stream sync.
 * ===================================================================== */

static int comm_rank_host(void *comm, int *rank) { return impl_rank(comm, rank); }
static int comm_size_host(void *comm, int *size) { return impl_size(comm, size); }

static int comm_send_host(const void *buffer, int count, cudssDataType_t dt,
                          int dest, int tag, void *comm, cudaStream_t stream)
{ return impl_send(buffer, count, dt, dest, tag, comm, stream, 0); }

static int comm_recv_host(void *buffer, int count, cudssDataType_t dt,
                          int root, int tag, void *comm, cudaStream_t stream)
{ return impl_recv(buffer, count, dt, root, tag, comm, stream, 0); }

static int comm_bcast_host(void *buffer, int count, cudssDataType_t dt,
                           int root, void *comm, cudaStream_t stream)
{ return impl_bcast(buffer, count, dt, root, comm, stream, 0); }

static int comm_reduce_host(const void *sb, void *rb, int count,
                            cudssDataType_t dt, cudssOpType_t op, int root,
                            void *comm, cudaStream_t stream)
{ return impl_reduce(sb, rb, count, dt, op, root, comm, stream, 0); }

static int comm_allreduce_host(const void *sb, void *rb, int count,
                               cudssDataType_t dt, cudssOpType_t op,
                               void *comm, cudaStream_t stream)
{ return impl_allreduce(sb, rb, count, dt, op, comm, stream, 0); }

static int comm_scatterv_host(const void *sb, const int *sc, const int *displs,
                              cudssDataType_t st, void *rb, int rc,
                              cudssDataType_t rt, int root, void *comm,
                              cudaStream_t stream)
{ return impl_scatterv(sb, sc, displs, st, rb, rc, rt, root, comm, stream, 0); }

static int comm_split_host(const void *comm, int color, int key, void *newcomm)
{ return impl_split(comm, color, key, newcomm); }

static int comm_free_host(void *comm) { return impl_free(comm); }

/* Report the comm-layer version. cuDSS requires the layer version be
 * >= the cuDSS library version, so report the headers we built against
 * (0.8.0). Return 0 on success. */
static int comm_get_property(libraryPropertyType type, int *value)
{
    if (!value) return CUDSS_COMM_ERR;
    switch (type) {
        case MAJOR_VERSION: *value = CUDSS_VERSION_MAJOR; break;
        case MINOR_VERSION: *value = CUDSS_VERSION_MINOR; break;
        case PATCH_LEVEL:   *value = CUDSS_VERSION_PATCH; break;
        default:            *value = 0; break;
    }
    return CUDSS_COMM_OK;
}

/* The symbol cuDSS looks up by name. Must be named exactly
 * `cudssDistributedInterface` and fully populated (21 fields in 0.8.0).
 * Designated initializers keep the mapping robust to field reordering. */
cudssDistributedInterface_t cudssDistributedInterface = {
    .cudssCommRankDevice      = comm_rank_dev,
    .cudssCommSizeDevice      = comm_size_dev,
    .cudssSendDevice          = comm_send_dev,
    .cudssRecvDevice          = comm_recv_dev,
    .cudssBcastDevice         = comm_bcast_dev,
    .cudssReduceDevice        = comm_reduce_dev,
    .cudssAllreduceDevice     = comm_allreduce_dev,
    .cudssScattervDevice      = comm_scatterv_dev,
    .cudssCommSplitDevice     = comm_split_dev,
    .cudssCommFreeDevice      = comm_free_dev,
    .cudssCommRankHost        = comm_rank_host,
    .cudssCommSizeHost        = comm_size_host,
    .cudssSendHost            = comm_send_host,
    .cudssRecvHost            = comm_recv_host,
    .cudssBcastHost           = comm_bcast_host,
    .cudssReduceHost          = comm_reduce_host,
    .cudssAllreduceHost       = comm_allreduce_host,
    .cudssScattervHost        = comm_scatterv_host,
    .cudssCommSplitHost       = comm_split_host,
    .cudssCommFreeHost        = comm_free_host,
    .cudssDistributedGetProperty = comm_get_property
};
