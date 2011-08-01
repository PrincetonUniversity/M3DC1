#include <m3dc1_source.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>

double factor = 1;
int NR, NZ, NPHI;
double rmin, rmax, phimin, phimax, zmin, zmax;
const int precision = 5;              // Precision of output
const int width = precision + 8;

std::ostream& write_header(std::ostream&);

int main(int argc, const char* argv[])
{
  const int MAX_SIZE = 100000000;       // Maximum allowed number of points
  const char filename_in[] = "C1.h5";   // Input file
  const char filename_bout[] = "B.out"; // Magnetic field output file
  const char filename_pout[] = "p.out"; // Pressure output file
  const char filename_nout[] = "ne.out"; // Electron density output file

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
  if(argc>=5) factor = atof(argv[4]); 

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
  source1.factor = factor;

  std::cerr << "loading ..." << std::endl;
  if(!source0.load()) {
    return 1;
  }
  if(!source1.load()) {
    return 1;
  }
  std::cerr << "done." << std::endl;

  std::cout << "NR, NZ, PHI = " << NR << " " << NZ << " " << NPHI << std::endl;
  std::cout << "factor = " << source1.factor << std::endl;

  // Read extent of domain
  source0.extent(&rmin, &rmax, &zmin, &zmax);
  phimin = 0.;
  phimax = 2.*M_PI*(double)(NPHI-1)/(double)NPHI;

  double dr   = (rmax - rmin)/(double)(NR-1);
  double dphi = (phimax - phimin)/(double)(NPHI-1);
  double dz   = (zmax - zmin)/(double)(NZ-1);

  // Open file for output
  std::fstream file_bout, file_pout, file_nout;
  file_bout.open(filename_bout, std::fstream::out | std::fstream::trunc);
  file_pout.open(filename_pout, std::fstream::out | std::fstream::trunc);
  file_nout.open(filename_nout, std::fstream::out | std::fstream::trunc);

  write_header(file_bout);
  write_header(file_pout);
  write_header(file_nout);

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
	double p = 0;
	double n = 0;
	source0.eval(R, Phi, Z, &Br, &Bphi, &Bz);
	source1.eval(R, Phi, Z, &Br, &Bphi, &Bz);
	file_bout << std::setw(width) << Br 
		  << std::setw(width) << Bz 
		  << std::setw(width) << Bphi << '\n';
	source0.eval_pn(R, Phi, Z, &p, &n);
	source1.eval_pn(R, Phi, Z, &p, &n);
	file_pout << std::setw(width) << p << '\n';
	file_nout << std::setw(width) << n << '\n';
      }
    }
  }

  file_bout.close();
  file_pout.close();
  file_nout.close();

  std::cout << "Done.\n" 
	    << "Field data output to " << filename_bout << "\n"
	    << "Pressure data output to " << filename_pout << "\n"
	    << "Density data output to " << filename_nout << std::endl;

  return 0;
}

std::ostream& write_header(std::ostream& os)
{
  os.setf(std::ios::scientific,std::ios::floatfield);
  os.precision(precision);

  // Write array dimensions
  os << std::setw(5) << NR  
     << std::setw(5) << NZ 
     << std::setw(5) << NPHI << std::endl;
  
  // Write physical dimensions
  os << std::setw(width) << rmin 
     << std::setw(width) << rmax
     << std::setw(width) << zmin 
     << std::setw(width) << zmax 
     << std::setw(width) << phimin 
     << std::setw(width) << phimax << std::endl;

  return os;
}
