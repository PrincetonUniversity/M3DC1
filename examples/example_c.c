//==========================================================================
// fio_example
// ~~~~~~~~~~~
// In this example, several files from linear calculations are loaded.
// The fields are evaluated at several points, and the results are summed.
//==========================================================================

#include <fusion_io_defs.h>
#include <stdio.h>

//const int nfiles = 4;
const int nfiles = 1;
int isrc[nfiles];
int ipres[nfiles], ine[nfiles], ini[nfiles];
int imag[nfiles], ij[nfiles];

void free_sources();

int main()
{
  // factor by which to multiply perturbed (linear) part of solution
  const double factor = 1.;

  // time slice to read
  const int timeslice = 1;

  char filename[nfiles][256];

  double p, ne, ni, b[3], x[3], x2[3], curr[3], db[9], db_n[9];
  double t1, t3[3], t32[3], t9[9];
  const double epsilon = 1e-6;

  const int npts = 10;
  double R0, R1, Z0, Z1, phi0, phi1;
  int i, j, ierr;

  double psi_axis, psi_lcfs;
  int ipsi_axis, ipsi_lcfs;

  strcpy(filename[0], "/p/tsc/nferraro/data/DIII-D/146626/2250/mesh22aa_kap6_amu6_d0_extsub/C1.h5");
/*
  strcpy(filename[0], "/p/tsc/nferraro/data/DIII-D/126006/3600_efit06/orlov/mesh21a_kap6_amu6_n=1/C1.h5");
  strcpy(filename[1], "/p/tsc/nferraro/data/DIII-D/126006/3600_efit06/orlov/mesh21a_kap6_amu6_n=2/C1.h5");
  strcpy(filename[2], "/p/tsc/nferraro/data/DIII-D/126006/3600_efit06/orlov/mesh21a_kap6_amu6_n=3/C1.h5");
  strcpy(filename[3], "/p/tsc/nferraro/data/DIII-D/126006/3600_efit06/orlov/mesh21a_kap6_amu6_n=4/C1.h5");
*/
  // read files and fields
  for(i=0; i<nfiles; i++) {
    printf("Reading %s\n", filename[i]);
    ierr =  fio_open_source(FIO_M3DC1_SOURCE, filename[i], &(isrc[i]));
    if(ierr != 0) {
      printf("Error reading %s\n", filename[i]);
      free_sources();
      return 1;
    }
    
    // Set options appropriate to this source
    ierr = fio_get_options(isrc[i]);
    ierr = fio_set_int_option(FIO_TIMESLICE, timeslice);
    ierr = fio_set_real_option(FIO_LINEAR_SCALE, factor);

    // For first file, read full fields (equilibrium + perturbed)
    // For subsequent file, read only perturbed parts
    // if(i > 0) ierr = fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);
    ierr = fio_set_int_option(FIO_PART, FIO_PERTURBED_ONLY);

    // read fields
    // magnetic field and total pressure are species-independent
    ierr = fio_get_field(isrc[i], FIO_TOTAL_PRESSURE, &(ipres[i]));
    ierr = fio_get_field(isrc[i], FIO_MAGNETIC_FIELD, &(imag[i]));
    ierr = fio_get_field(isrc[i], FIO_CURRENT_DENSITY, &(ij[i]));

     // density is species dependent; specify electrons
    ierr = fio_set_int_option(FIO_SPECIES, FIO_ELECTRON);
    ierr = fio_get_field(isrc[i], FIO_DENSITY, &(ine[i]));
    
    ierr = fio_set_int_option(FIO_SPECIES, FIO_MAIN_ION);
    ierr = fio_get_field(isrc[i], FIO_DENSITY, &(ini[i]));
  }
 
  ierr = fio_get_series(isrc[0], FIO_MAGAXIS_PSI, &ipsi_axis);
  ierr = fio_get_series(isrc[0], FIO_LCFS_PSI, &ipsi_lcfs);
  printf("ipsi_axis: %d\n", ipsi_axis);
  printf("ipsi_lcfs: %d\n", ipsi_lcfs);

  ierr = fio_eval_series(ipsi_axis, 0.d, &psi_axis);
  ierr = fio_eval_series(ipsi_lcfs, 0.d, &psi_lcfs);
  printf("Psi at magnetic axis: %g\n", psi_axis);
  printf("Psi at lcfs: %g\n", psi_lcfs);
    
  ierr = fio_close_series(ipsi_axis);
  ierr = fio_close_series(ipsi_lcfs, ierr);

  R0 = 1.6;
  R1 = 2.1;
  Z0 = 0.4;
  Z1 = 0.4;
  phi0 = 0.;
  phi1 = 0.;

  for(i=0; i<npts; i++) {
    x[0] = R0 + (R1-R0)*(double)i/(double)(npts-1);
    x[1] = phi0 + (phi1-phi0)*(double)i/(double)(npts-1);
    x[2] = Z0 + (Z1-Z0)*(double)i/(double)(npts-1);

    printf("(%g, %g, %g)\n", x[0], x[1], x[2]);

    p = 0.;
    ne = 0.;
    ni = 0.;
    for(j=0; j<3; j++) {
      b[j] = 0.;
      curr[j] = 0.;
    }
    for(j=0; j<9; j++) { 
      db[j] = 0.;
      db_n[j] = 0.;
    }

    for(j=0; j<nfiles; j++) {
      ierr = fio_eval_field(ipres[j], x, &t1);
      p += t1;
      ierr = fio_eval_field(ine[j], x, &t1);
      ne += t1;
      ierr = fio_eval_field(ini[j], x, &t1);
      ni += t1;
      ierr = fio_eval_field(imag[j], x, t3);
      b[0] += t3[0];
      b[1] += t3[1];
      b[2] += t3[2];

      // calculate numerical derivatives of B
      // dR
      x2[0] = x[0] + epsilon;  x2[1] = x[1];  x2[2] = x[2];
      ierr = fio_eval_field(imag[j], x2, t32);
      db_n[FIO_DR_R] += (t32[0] - t3[0]) / epsilon;
      db_n[FIO_DR_PHI] += (t32[1] - t3[1]) / epsilon;
      db_n[FIO_DR_Z] += (t32[2] - t3[2]) / epsilon;
      // dPHI
      x2[0] = x[0];  x2[1] = x[1] + epsilon;  x2[2] = x[2];
      ierr = fio_eval_field(imag[j], x2, t32);
      db_n[FIO_DPHI_R] += (t32[0] - t3[0]) / epsilon;
      db_n[FIO_DPHI_PHI] += (t32[1] - t3[1]) / epsilon;
      db_n[FIO_DPHI_Z] += (t32[2] - t3[2]) / epsilon;
      // dZ
      x2[0] = x[0];  x2[1] = x[1];  x2[2] = x[2] + epsilon;
      ierr = fio_eval_field(imag[j], x2, t32);
      db_n[FIO_DZ_R] += (t32[0] - t3[0]) / epsilon;
      db_n[FIO_DZ_PHI] += (t32[1] - t3[1]) / epsilon;
      db_n[FIO_DZ_Z] += (t32[2] - t3[2]) / epsilon;     

      ierr = fio_eval_field(ij[j], x, t3);
      curr[0] += t3[0];
      curr[1] += t3[1];
      curr[2] += t3[2];

      ierr = fio_eval_field_deriv(imag[j], x, t9);
      for(j=0; j<9; j++)
	db[j] += t9[j];

    }

    printf("\tpressure = %g\n", p);
    printf("\telectron density = %g\n", ne);
    printf("\tion density = %g\n", ni);
    printf("\ttotal B = (%g, %g, %g)\n", b[0], b[1], b[2]);
    printf("\ttotal J = (%g, %g, %g)\n", curr[0], curr[1], curr[2]);

    printf("Interpolated derivatives: \n");
    printf("\tdB/dR = (%g, %g, %g)\n", 
	   db[FIO_DR_R], db[FIO_DR_PHI], db[FIO_DR_Z]);
    printf("\tdB/dphi = (%g, %g, %g)\n", 
	   db[FIO_DPHI_R], db[FIO_DPHI_PHI], db[FIO_DPHI_Z]);
    printf("\tdB/dZ = (%g, %g, %g)\n", 
	   db[FIO_DZ_R], db[FIO_DZ_PHI], db[FIO_DZ_Z]);

    printf("Numerical derivatives: \n");
    printf("\tdB/dR = (%g, %g, %g)\n", 
	   db_n[FIO_DR_R], db_n[FIO_DR_PHI], db_n[FIO_DR_Z]);
    printf("\tdB/dphi = (%g, %g, %g)\n", 
	   db_n[FIO_DPHI_R], db_n[FIO_DPHI_PHI], db_n[FIO_DPHI_Z]);
    printf("\tdB/dZ = (%g, %g, %g)\n", 
	   db_n[FIO_DZ_R], db_n[FIO_DZ_PHI], db_n[FIO_DZ_Z]);
  }

  free_sources();

  return 0;
}

void free_sources()
{
  int i, ierr;
  for(i=0; i<nfiles; i++) {
    ierr = fio_close_field(ine[i]);
    ierr = fio_close_field(ini[i]);
    ierr = fio_close_field(ipres[i]);
    ierr = fio_close_field(imag[i]);
    ierr = fio_close_field(ij[i]);
    ierr = fio_close_source(isrc[i]);
  }
}
