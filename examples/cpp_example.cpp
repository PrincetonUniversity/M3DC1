#include <iostream>
#include <m3dc1_file.h>

int main()
{
  m3dc1_file file;
  const char filename[] = "C1.h5";
  const int timeslice = 1;
  const double r0=1.3, z0=0., phi0=0.;
  const double dr=0.05, dz=0., dphi=0.;
  const int npts = 10;
  

  // Open m3dc1 file
  if(!file.open(filename)) {
    std::cerr << "Error opening " << filename << std::endl;
    return 1;
  }

  // read toroidal field
  double rzero, bzero;
  file.read_parameter("rzero", &rzero);
  file.read_parameter("bzero", &bzero);
  double fzero = rzero*bzero;
  std::cout << "fzero = " << fzero << std::endl;

  // load field data
  m3dc1_field* psi = file.load_field("psi",timeslice,m3dc1_file::M3DC1_TOTAL);
  m3dc1_field* f = file.load_field("f",timeslice,m3dc1_file::M3DC1_TOTAL);

  // evaluate fields
  double psival[m3dc1_field::OP_NUM];
  double fval[m3dc1_field::OP_NUM];

  const m3dc1_field::m3dc1_get_op psiop = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_VAL | m3dc1_field::GET_DVAL | m3dc1_field::GET_PVAL);
  const m3dc1_field::m3dc1_get_op fop = 
    (m3dc1_field::m3dc1_get_op) 
    (m3dc1_field::GET_DVAL | m3dc1_field::GET_DDVAL | 
     m3dc1_field::GET_PVAL);

  
  for(int i=0; i<npts; i++) {
    double r = r0 + i*dr;
    double z = z0 + i*dz;
    double phi = phi0 + i*dphi;

    std::cerr << "Evaluating at (" 
	      << r << ", " << phi << ", " << z << ")" << std::endl;

    if(!psi->eval(r, phi, z, psiop, psival)) {
      std::cerr << "Error reading psi at" << std::endl;
    }
    if(!f->eval(r, phi, z, fop, fval)) {
      std::cerr << "Error reading f"  << std::endl;
    }

    // evaluate vector potential
    // A = R^2 Grad(Phi) x Grad(f) + psi Grad(Phi) - F0 Ln(R) Grad(Z)

    double AR   =  r*fval[m3dc1_field::OP_DZ];
    double AZ   = -r*fval[m3dc1_field::OP_DR] - fzero*log(r);
    double APHI =  psival[m3dc1_field::OP_1]/r;
    std::cout << " A = (" << AR << ", " << APHI << ", " << AZ << ")" 
	      << std::endl;

    // evaluate derivatives of vector potential
    double AR_R     =  r*fval[m3dc1_field::OP_DRZ] + fval[m3dc1_field::OP_DZ];
    double AR_Z     =  r*fval[m3dc1_field::OP_DZZ];
    double AR_PHI   =  r*fval[m3dc1_field::OP_DZP];
    
    double AZ_R     = -r*fval[m3dc1_field::OP_DRR] - fval[m3dc1_field::OP_DR]
      - fzero/r;
    double AZ_Z     = -r*fval[m3dc1_field::OP_DRZ];
    double AZ_PHI   = -r*fval[m3dc1_field::OP_DRP];
   
    double APHI_R   =  psival[m3dc1_field::OP_DR]/r 
      - psival[m3dc1_field::OP_1]/(r*r);
    double APHI_Z   =  psival[m3dc1_field::OP_DZ]/r;
    double APHI_PHI =  psival[m3dc1_field::OP_DP]/r;
    
    std::cout << " dA/dR = (" 
	      << AR_R   << ", " << APHI_R   << ", " << AZ_R << ")\n"
	      << " dA/dZ = (" 
	      << AR_Z   << ", " << APHI_Z   << ", " << AZ_Z << ")\n"
	      << " dA/dPhi=(" 
	      << AR_PHI << ", " << APHI_PHI << ", " << AZ_PHI << ")"
	      << std::endl;
  }

  // Close m3dc1 file
  file.close();

  return 0;
}
