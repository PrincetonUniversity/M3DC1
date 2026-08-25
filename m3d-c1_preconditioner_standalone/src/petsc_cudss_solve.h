// Copyright (c) 2026 NVIDIA Corporation & Affiliates. All rights reserved.

/*
 * petsc_cudss_solve.h -- reusable entry point for the cuDSS block-Jacobi
 * PCShell solve implemented in petsc_cudss_solve.c.
 *
 * Solves A x = b with the KSP configured from the PETSc options database
 * (e.g. -ksp_type lgmres), preconditioned by a PCShell whose Apply runs a
 * cuDSS LU solve on each plane diagonal block (block-Jacobi across planes).
 *
 * Requirements:
 *   - plane-aligned distribution: nplanes > 0, nprocs % nplanes == 0, each
 *     rank's ownership range lies entirely inside one plane's row block;
 *   - A/b/x use GPU types (-mat_type aijcusparse -vec_type cuda);
 *   - real double scalars (cuDSS path).
 *
 * setKspType_cudss is the cuDSS analogue of matrix_solve::setKspType():
 * one-time setup that extracts the plane diagonal block from A, runs cuDSS
 * ANALYSIS+FACTORIZATION, creates the outer KSP (operators + tolerances +
 * KSPSetFromOptions) and installs the cuDSS block-Jacobi PCShell on it.
 * The returned KSP owns the cuDSS block and plane communicator; both are
 * freed by KSPDestroy.  m3dc1_scorec calls it from
 * matrix_solve::solve_cudss / solve_cudss_with_guess when _kspSet == 0 and
 * reuses the KSP across solves (the factorization is NOT refreshed when A's
 * values change); the standalone driver calls it per rep.
 *
 * petsc_cudss_solve solves A x = b with a KSP previously set up by
 * setKspType_cudss.  It only handles the initial guess, KSPSolve and
 * reporting; nothing is destroyed.
 *
 * use_initial_guess: if PETSC_TRUE, the incoming x is used as the nonzero
 * initial guess; otherwise x is zeroed.
 * out_its: if non-NULL, receives the KSP iteration count.
 */
#ifndef PETSC_CUDSS_SOLVE_H
#define PETSC_CUDSS_SOLVE_H

#include <petsc.h>

#ifdef __cplusplus
extern "C" {
#endif

PetscErrorCode setKspType_cudss(Mat A, PetscInt nplanes, KSP *out_ksp);

PetscErrorCode petsc_cudss_solve(KSP ksp, Mat A, Vec b, Vec x,
                                 PetscInt nplanes,
                                 PetscBool use_initial_guess,
                                 PetscInt *out_its);

#ifdef __cplusplus
}
#endif

#endif /* PETSC_CUDSS_SOLVE_H */
