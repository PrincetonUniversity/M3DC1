// kokkos_gemm.cc
// Kokkos replacement for the BLAS dgemm call in define_basis().
//
// Computes  C(M x N) = A(M x K) * B^T(K x N)
// where A, B, C are Fortran column-major (LayoutLeft) double arrays.
// This is equivalent to:
//   dgemm('N','T', M, N, K, 1.0, A, M, B, N, 0.0, C, M)
//
// Callable from Fortran via iso_c_binding:
//   interface
//     subroutine kokkos_dgemm_nt(M,N,K,A,B,C) bind(C,name='kokkos_dgemm_nt')
//       use iso_c_binding
//       integer(c_int), value, intent(in) :: M, N, K
//       real(c_double), intent(in)        :: A(*), B(*)
//       real(c_double), intent(out)       :: C(*)
//     end subroutine
//   end interface
//
// Build note: requires -DUSEKOKKOS in FOPTS and linking against libkokkos.
// On CPU, DefaultExecutionSpace maps to OpenMP (with -DKokkos_ENABLE_OPENMP)
// or Serial; on GPU it maps to CUDA/HIP.  The A/B/C pointers must reside in
// Kokkos::HostSpace (consistent with the USEBLAS CPU path).

#include <Kokkos_Core.hpp>

extern "C" {

// kokkos_dgemm_nt: compute C = A * B^T  (N-T variant)
// Arguments are passed by value (Fortran VALUE attribute in interface).
void kokkos_dgemm_nt(int M, int N, int K,
                     const double* A,   // M x K, column-major
                     const double* B,   // N x K, column-major
                     double*       C)   // M x N, column-major (output)
{
  using HostExec = Kokkos::DefaultHostExecutionSpace;

  // Use raw pointer arithmetic to avoid Kokkos::View::operator() overload
  // issues inside KOKKOS_LAMBDA (lambda captures are const, causing int&
  // mismatch with MemoryTraits<Unmanaged> views).
  // Column-major (LayoutLeft / Fortran) indexing:
  //   A(i,l) = A[i + M*l],  B(j,l) = B[j + N*l],  C(i,j) = C[i + M*j]

  // C(i,j) = sum_l  A(i,l) * B(j,l)
  // Outer (i,j) iterations are parallelised; inner K-loop is sequential.
  Kokkos::parallel_for(
    "kokkos_dgemm_nt",
    Kokkos::MDRangePolicy<HostExec, Kokkos::Rank<2>>({0, 0}, {M, N}),
    KOKKOS_LAMBDA(int i, int j) {
      double sum = 0.0;
      for (int l = 0; l < K; ++l)
        sum += A[i + M*l] * B[j + N*l];
      C[i + M*j] = sum;
    });

  Kokkos::fence();
}

} // extern "C"
