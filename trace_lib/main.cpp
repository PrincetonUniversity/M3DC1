#include "trace_integrator.h"
#include "m3dc1_source.h"
#include "geqdsk_source.h"
#include "diiid_coils.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <time.h>
#include <mpi.h>

trace_integrator tracer;
int transits = 100;
int steps_per_transit = 100;
int surfaces = 11;
double dR = 0.;
double dZ = 0.;
double dR0 = 0.;
double dZ0 = 0.;
double Phi = 0.;
double Phi0 = 0.;
double R0 = 0.;
double Z0 = 0.;
double angle = 0.;
bool qout = true;
bool pout = true;
double ds = 0.;

bool R0_set = false;
bool Z0_set = false;
double phase_set = false;


void print_help();
void print_parameters();
bool process_command_line(int argc, char* argv[]);
bool process_line(const std::string& opt, const int argc, 
		  const std::string argv[]);
void delete_sources();

int my_offset(const int rank, const int size);
int my_size(const int rank, const int size);

int main(int argc, char* argv[])
{
  int size, rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(!process_command_line(argc, argv)) {
    print_help();
    return 1;
  }

  tracer.plane = angle;

  if(tracer.sources.size() == 0) {
    std::cerr << "No sources specified.  Returning." << std::endl;
    print_help();
    return 1;
  }

  if(!tracer.load()) {
    std::cerr << " Error loading tracer" << std::endl;
    return 1;
  }

  double R_axis, Z_axis;
  tracer.center(&R_axis, &Z_axis);

  std::cerr << "Magnetic axis is at ( " << R_axis << ", "
	    << Z_axis << ")" << std::endl;

  if(!R0_set) R0 = R_axis;
  if(!Z0_set) Z0 = Z_axis;

  if(dR0<=0.) dR0 = dR;

  print_parameters();

  std::fstream plotfile, plotfile2;
  char outfile[256];

  if(pout) {
    plotfile.open("gplot", std::fstream::out | std::fstream::trunc);
    plotfile << "set size ratio -1" << '\n'
	     << "set xlabel 'R'" << '\n'
	     << "set ylabel 'Z'" << '\n'
	     << "set pointsize 0.1" << '\n'
	     << "plot ";

    for(int j=0; j<surfaces; j++) {
      sprintf(outfile, "out%d", j);     
      if(j>0) plotfile << ",\\\n";
      plotfile << "'" << outfile << "'" << " title ''";
    }
    plotfile << std::endl;
    plotfile.close();

    plotfile2.open("gplotpsi", std::fstream::out | std::fstream::trunc);
    plotfile2 << "set size ratio 1" << '\n'
	      << "set xlabel 'theta'" << '\n'
	      << "set ylabel 'psi'" << '\n'
	      << "set pointsize 0.1" << '\n'
	      << "plot [-180:180] ";

    for(int j=0; j<surfaces; j++) {
      sprintf(outfile, "out%d", j);     
      if(j>0) plotfile2 << ",\\\n";
      plotfile2 << "'" << outfile << "'" << " using 3:4 title ''";
    }
    plotfile2 << std::endl;
    plotfile2.close();
  }

  int* offsets = new int[size];
  int* local_surfaces = new int[size];

  for(int i=0; i<size; i++) {
    offsets[i] = my_offset(i, size);
    local_surfaces[i] = my_size(i, size);
  }

  double *my_q_min, *my_q_max, *my_q_mean, *my_r, *my_phase;
  my_r = new double[local_surfaces[rank]];
  my_q_min = new double[local_surfaces[rank]];
  my_q_max = new double[local_surfaces[rank]];
  my_q_mean = new double[local_surfaces[rank]];
  my_phase = new double[local_surfaces[rank]];

  time_t t1 = time(0);
  for(int j=0; j<local_surfaces[rank]; j++) {
    int s = j + offsets[rank];
    double R = R0 + dR0 + dR*s;
    double Z = Z0 + dZ0 + dZ*s;
    double *rr, *zz;
    bool result;

    int n;

    if(pout)  {
      sprintf(outfile, "out%d", s);
      tracer.open_file(outfile);
    }

    if(ds==0.) {
      n = 1;
      rr = new double[1];
      zz = new double[1];
      rr[0] = R;
      zz[0] = Z;
    } else { 
      if(!tracer.get_surface(R, Phi, Z, ds, &rr, &zz, &n))
	std::cerr << "Warning: not on closed surface" << std::endl;

      if(n==0) {
	std::cerr << "Error: no points found on surface" << std::endl;
	break;
      }
    }

    my_r[j] = sqrt((R - R_axis)*(R - R_axis) + (Z - Z_axis)*(Z - Z_axis));
    my_phase[j] = Phi;

    std::cerr << "Surface " << s+1 << " of " << surfaces
	      << " (" << n << " pts at this surface)" << std::endl;
    std::cerr << "(" << R << "," << Phi <<  ", "  << Z << ") ... ";

    // Cycle through surface points
    trace_integrator::integrator_data data;
    my_q_mean[j] = 0.;
    my_q_max[j] = 0.;
    my_q_min[j] = 0.;
    int tt = 0;
    result = false;
    for(int k=0; k<n; k++) {   
      if(phase_set) Phi = Phi0;
      else Phi = tracer.find_min_bn(rr[k],zz[k]);

      std::cerr << atan2(zz[k] - Z_axis, rr[k] - R_axis) << "\t";
      //      dfile << rr[k] << '\t' <<  zz[k] << '\t';
      tracer.set_pos(rr[k],Phi,zz[k]);

      // perform integration
      result = tracer.integrate(transits, steps_per_transit, &data);
      if(!result) {
	double RR, PP, ZZ;
	tracer.get_pos(&RR, &PP, &ZZ);
	std::cerr << "lost at (" << RR << ", " << ZZ <<  ")" << std::endl;
      }

      // write connection length
      //      if(result)
      //	dfile << 0. << std::endl;
      //      else
      //      dfile << data.distance << '\t' << log10(data.distance) << std::endl;

      if(k==0 || data.q < my_q_min[j]) my_q_min[j] = data.q;
      if(k==0 || data.q > my_q_max[j]) my_q_max[j] = data.q;
				   
      my_q_mean[j] += data.q*data.toroidal_transits;
      tt += data.toroidal_transits;
    }
    my_q_mean[j] /= (double)tt;
    
    std::cerr << tt << " tor. transits, q = " << my_q_mean[j] << std::endl;

    tracer.close_file();

    delete[] rr;
    delete[] zz;
  }
  time_t t2 = time(0);
  std::cerr << "Computation completed in " << t2 - t1 << " s." << std::endl;

  
  // assemble and write data
  // ~~~~~~~~~~~~~~~~~~~~~~~

  // allocate global data
  double *total_r, *total_q_min, *total_q_max, *total_q_mean, *total_phase;

  if(rank==0) {
    total_r = new double[surfaces];
    if(qout) {
      total_q_min = new double[surfaces];
      total_q_max = new double[surfaces];
      total_q_mean = new double[surfaces];
      total_phase = new double[surfaces];
    }
  }
  
  MPI_Gatherv(my_r, local_surfaces[rank], MPI_DOUBLE, 
	      total_r, local_surfaces, offsets, MPI_DOUBLE, 
	      0, MPI_COMM_WORLD);
  if(qout) {
    MPI_Gatherv(my_q_min, local_surfaces[rank], MPI_DOUBLE, 
		total_q_min, local_surfaces, offsets, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_q_max, local_surfaces[rank], MPI_DOUBLE, 
		total_q_max, local_surfaces, offsets, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_q_mean, local_surfaces[rank], MPI_DOUBLE, 
		total_q_mean, local_surfaces, offsets, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
    MPI_Gatherv(my_phase, local_surfaces[rank], MPI_DOUBLE, 
		total_phase, local_surfaces, offsets, MPI_DOUBLE, 
		0, MPI_COMM_WORLD);
  }

  if(rank==0) {
    // write data
    if(qout) {
      std::fstream qfile;     
      qfile.open("q.out", std::fstream::out | std::fstream::trunc);
      for(int i=0; i<surfaces; i++) {
	qfile << total_r[i] << '\t' 
	      << total_q_mean[i] << '\t'
	      << total_q_min[i] << '\t'
	      << total_q_max[i] << '\n';
      }
      qfile.close();
    }
    
    // delete global data
    delete[] total_r;
    if(qout) {
      delete[] total_q_min;
      delete[] total_q_max;
      delete[] total_q_mean;
      delete[] total_phase;
    }
  }

  // delete local data
  delete[] my_r;
  delete[] my_q_min;
  delete[] my_q_max;
  delete[] my_q_mean;
  delete[] my_phase;

  delete[] offsets;
  delete[] local_surfaces;

  std::cerr << "Deleting sources." << std::endl;
  delete_sources();

  std::cerr << "===============================\n" 
	    << "Poincare computation complete.\n"
	    << "To view Poincare plot, use 'gnuplot gplot'" << std::endl;

  MPI_Finalize();

  return 0;
}

void delete_sources()
{
  trace_source_list::iterator i;
  i = tracer.sources.begin();

  while(i != tracer.sources.end()) {
    delete *i;
    i++;
  }
}

bool process_command_line(int argc, char* argv[])
{
  const int max_args = 4;
  const int num_opts = 15;
  std::string arg_list[num_opts] = 
    { "-geqdsk", "-m3dc1", "-diiid-i",
      "-dR", "-dZ", "-dR0", "-dZ0", 
      "-ds", "-p", "-t", "-s", "-a",
      "-pout", "-qout", "-phi0"};
  std::string opt = "";
  std::string arg[max_args];
  int args = 0;
  bool is_opt;
  bool processed = true;

  for(int i=1; i<argc; i++) {
    // determine if current cl arg is an option
    is_opt = false;
    for(int j=0; j<num_opts; j++) {
      if(argv[i]==arg_list[j]) {
	is_opt = true;
	break;
      }
    }
    
    if(is_opt) {     // if so, process current option
      if(!processed)
	if(!process_line(opt, args, arg)) return false;

      opt = argv[i];
      args = 0;
      processed = false;
      
    } else {         // otherwise, add argument
      if(args >= max_args)
	std::cerr << "Too many arguments for " << opt << std::endl;
      else 
	arg[args++] = argv[i];
    }
  }

  if(!processed)
    if(!process_line(opt, args, arg)) return false;

  return true;

}

int my_offset(const int rank, const int size)
{
  return rank*(surfaces/size);
}

int my_size(const int rank, const int size)
{
  int min = my_offset(rank, size);
  int max = (rank==size-1) ? surfaces : my_offset(rank+1, size);
  return max-min;
}


bool process_line(const std::string& opt, const int argc, const std::string argv[])
{
  bool argc_err = false;

  if(opt=="-m3dc1") {
    m3dc1_source* s = new m3dc1_source();
    if(argc>=1) s->filename = argv[0];
    if(argc>=2) s->time = atoi(argv[1].c_str());
    if(argc>=3) s->factor = atof(argv[2].c_str());
    if(argc>=4) s->shift = atof(argv[3].c_str())*M_PI/180.;
    tracer.sources.push_back(s);
  } else if(opt=="-geqdsk") {
    geqdsk_source* s = new geqdsk_source();
    if(argc>=1) s->filename = argv[0];
    tracer.sources.push_back(s);
  } else if(opt=="-diiid-i") {
    coil_source* s = new coil_source();
    double current;
    int n;
    double phase;
    if(argc>=1) current = atof(argv[0].c_str()); else current = 1e3;
    if(argc>=2) n = atoi(argv[1].c_str());       else n = 1;
    if(argc>=3) phase = atof(argv[2].c_str());   else phase = 0.;
    diiid_icoils(s, current, n, 0., phase);
    tracer.sources.push_back(s);
  } else if(opt=="-dR") {
    if(argc==1) dR = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dZ") {
    if(argc==1) dZ = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dR0") {
    if(argc==1) dR0 = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-dZ0") {
    if(argc==1) dZ0 = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-ds") {
    if(argc==1) ds = atof(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-p") {
    if(argc==1) surfaces = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-t") {
    if(argc==1) transits = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-s") {
    if(argc==1) steps_per_transit = atoi(argv[0].c_str());
    else argc_err = true;
  } else if(opt=="-a") {
    if(argc==1) angle = atof(argv[0].c_str())*M_PI/180.;
    else argc_err = true;
  } else if(opt=="-phi0") {
    if(argc==1) {
      Phi0 = atof(argv[0].c_str())*M_PI/180.;
      phase_set = true;
    } else argc_err = true;
  } else if(opt=="-pout") {
    if(argc==1) pout = (atoi(argv[0].c_str()) == 1);
    else argc_err = true;
  } else if(opt=="-qout") {
    if(argc==1) qout = (atoi(argv[0].c_str()) == 1);
    else argc_err = true;

  } else {
    std::cerr << "Unrecognized option " << opt << std::endl;
    return false;
  }

  if(argc_err) {
    std::cerr << "Incorrect number of arguments for option " 
	      << opt << std::endl;
    return false;
  }

  return true;
}


void print_help()
{
  std::cout 
    << "Usage:\n"
    << "\ttrace <sources> -dR <dR> -dZ <dZ> -dR0 <dR0> -dZ0 <dZ0> /\n"
    << "\t   -p <pts> -t <trans> -s <steps> -a <angle> -phi0 <phi0> \n"
    << "\t   -pplot <pplot> -qplot <qplot>\n\n"
    << " <angle>   The toroidal angle of the plane in degrees (default = 0)\n"
    << " <dR>      R-spacing of seed points (default = major radius/(2*pts))\n"
    << " <dR0>     R-distance of first seed point from axis (default = dR)\n"
    << " <dZ>      Z-spacing of seed points (default = 0)\n"
    << " <dZ0>     Z-distance of first seed point from axis (default = dZ)\n"
    << " <phi0>    Toroidal angle of all seed points in degrees\n"
    << " <pts>     Number of seed points (default = 11)\n"
    << " <steps>   Integration steps per toroidal transit (default = 100)\n"
    << " <pplot>   If 1, Generate poincare plot data (default = 1)\n"
    << " <qplot>   If 1, Generate q-profile data (default = 1)\n"
    << " <trans>   Toroidal transits per seed point (default = 100)\n"
    << " <sources> May be one or more of the following:\n"
    << "\n  -m3dc1 <c1_file> <ts> <factor> <shift>\n"
    << "   * Loads field information from M3D-C1 data file\n"
    << "   <c1_file> filename of M3D-C1 hdf5 file (default = C1.h5)\n" 
    << "   <ts>      timeslice of M3D-C1 hdf5 file (default = -1)\n"
    << "   <factor>  factor by which to multiply field (default = 1)\n"
    << "   <shift>   toroidal phase shift of field, in degrees (default = 0)\n"
    << "\n  -geqdsk <gfile>\n"
    << "   * Loads field information from EFIT g-file\n"
    << "   <gfile>   filename of EFIT g-file\n"
    << "\n  -diiid-i <curr> <n> <phase>\n"
    << "   * Calculates field from DIII-D I-coil set\n"
    << "   <curr>    current (in Amps)\n"
    << "   <n>       toroidal modenumber\n"
    << "   <phase>   phase difference between top and bottom coil sets (deg)\n"
    << "\nExample:\n"
    << "\ttrace -m3dc1 C1.h5 -m3dc1 C1.h5 1 1e4\n\n"
    << " Loads field from time slices -1 and 1 from M3D-C1 file 'C1.h5'.\n"
    << " The field from time slice 1 is multiplied by 1e4.\n"
    << std::endl;
}

void print_parameters()
{
  std::cerr << "dR0 = " << dR0 << '\n'
	    << "dR = " << dR << '\n'
	    << "dZ0 = " << dZ0 << '\n'
	    << "dZ = " << dZ << '\n'
	    << "Number of surfaces = " << surfaces << '\n'
	    << "Toroidal transits per seed = " << transits << '\n'
	    << "Integrator steps per transit = " << steps_per_transit 
	    << std::endl;
}
