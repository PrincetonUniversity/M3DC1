#include <m3dc1_source.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

int main(int argc, const char* argv[])
{
  const int MAX_SIZE = 100000000;       // Maximum allowed number of points
  const int precision = 5;              // Precision of output
  const char filename_in[] = "C1.h5";   // Input file
  const char filename_out[] = "B.out";  // Output file

  int NR, NZ, NPHI;

  if(argc < 4) {
    std::cout 
      << "Usage:  m3dc1_convert NR NZ NPHI\n"
      << " NR, NZ, NPHI: Number of points in R, Z, and PHI directions\n\n"
      << "Output (B.out) is in following format: \n"
      << " NR NZ NPHI\n"
      << " RMIN RMAX ZMIN ZMAX PHIMIN PHIMAX\n"
      << " BR1 BZ1 BPHI1\n"
      << " BR2 BZ2 BPHI2\n"
      << " BR3 BZ3 BPHI3\n"
      << " ...\n\n"
      << "where\n"
      << " RMIN, RMAX: minimum and maximum R coordinate\n"
      << " ZMIN, ZMAX: minimum and maximum Z coordinate\n"
      << " PHIMIN, PHIMAX: minimum and maximum PHI coordinate\n"
      << " BR1 BZ1 BPHI1: components of field at point 1\n\n"
      << "points are ordered in the following way:\n"
      << " i = ir + iz*NR + iphi*NR*NZ\n"
      << "coordinates of point i are\n"
      << " (R, Phi, Z) = \n"
      << " (ir  *(  RMAX-  RMIN)/(NR  -1) +   RMIN, \n"
      << "  iphi*(PHIMAX-PHIMIN)/(NPHI-1) + PHIMIN, \n"
      << "  iz  *(  ZMAX-  ZMIN)/(NZ  -1) +   ZMIN)"
      << std::endl;
      return 1;
  }
  
  NR = atoi(argv[1]);
  NZ = atoi(argv[2]);
  NPHI = atoi(argv[3]);

  if(NR < 2 || NZ < 2 || NPHI < 2) {
    std::cerr << "Error: invalid number of points" << std::endl;
    return 1;
  }
  if(NR*NZ*NPHI > MAX_SIZE) {
    std::cerr << "Error: number of points exceeds MAX_SIZE.\n"
	      << "       To increase limit, raise MAX_SIZE and recompile."
	      << std::endl;
    return 1;
  }

  m3dc1_source source0(filename_in, -1);
  m3dc1_source source1(filename_in,  1);

  std::cerr << "loading ..." << std::endl;
  if(!source0.load()) {
    return 1;
  }
  if(!source1.load()) {
    return 1;
  }

  std::cerr << "done." << std::endl;

  double rmin, rmax, phimin, phimax, zmin, zmax;

  // Read extent of domain
  source0.extent(&rmin, &rmax, &zmin, &zmax);
  phimin = 0.;
  phimax = 2.*M_PI*(double)(NPHI-1)/(double)NPHI;

  double dr   = (rmax - rmin)/(double)(NR-1);
  double dphi = (phimax - phimin)/(double)(NPHI-1);
  double dz   = (zmax - zmin)/(double)(NZ-1);

  // Open file for output
  std::fstream file_out;
  file_out.open(filename_out, std::fstream::out | std::fstream::trunc);
  file_out.setf(std::ios::scientific,std::ios::floatfield);
  file_out.precision(precision);

  // Write array dimensions
  file_out << std::setw(5) << NR  
	   << std::setw(5) << NZ 
	   << std::setw(5) << NPHI << std::endl;
  
  // Write physical dimensions
  const int width = precision + 8;
  file_out << std::setw(width) << rmin 
	   << std::setw(width) << rmax
	   << std::setw(width) << zmin 
	   << std::setw(width) << zmax 
	   << std::setw(width) << phimin 
	   << std::setw(width) << phimax << std::endl;

  // Write field data
  for(int i=0; i<NPHI; i++) {
    double Phi = dphi*i + phimin;
    std::cout << "Plane " << i+1 << " of " << NPHI << "..." << std::endl;
    for(int j=0; j<NZ; j++) {
      double Z = dz*j + zmin;
      for(int k=0; k<NR; k++) {
	double R = dr*k + rmin;
	double Br = 0;
	double Bz = 0;
	double Bphi = 0;
	source0.eval(R, Phi, Z, &Br, &Bphi, &Bz);
	source1.eval(R, Phi, Z, &Br, &Bphi, &Bz);
	file_out << std::setw(width) << Br 
		 << std::setw(width) << Bz 
		 << std::setw(width) << Bphi << '\n';
      }
    }
  }

  file_out.close();

  std::cout << "Done.  Field data output to " << filename_out << std::endl;

  return 0;
}
