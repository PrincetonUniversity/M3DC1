//////////////////////////////////////////////////////////////////
// make_polar
//
//
// Creates a POLAR file for a generically shaped domain.
// Notes: Elongation is taken to be ~(r/a)
//        Triangularity is taken to be ~(r/a)^2
//
// Author: Nathaniel Ferraro
// Date:   3/26/09
//////////////////////////////////////////////////////////////////

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// defaults
const double R0 = 2.;
const double a0 = 1.;
const double e0 = 1.;
const double t0 = 0.;
const double s0 = 0.;
const double w10 = 0.1;
const double z00 = 0.;
const int nr0 = 20;
const int np0 = 50;

const char filename[] = "POLAR";

double R, a, e, t, s, z0;
double s1, f1, g1, w1;
int nr, np;

void print_usage();
void print_parameters();
int read_inputs(int, char**);
int read_input_file(const char*);
double calc_r(double);
void calc_g1();

int main(int argc, char* argv[])
{
  int i, j, k;
  double r, theta, sn, co, dx, dz, ddx, ddz, eps;
  double x, z, v, norm1, norm2, curv;

  FILE* f;


  // defaults
  R = R0;
  a = a0;
  e = e0;
  t = t0;
  s = s0;
  nr = nr0;
  np = np0;
  s1 = 0;
  f1 = 0;
  w1 = w10;
  z0 = z00;

  // read command line parameters
  if(read_inputs(argc,argv) != 0) {
    print_usage();
    return 0;
  }

  // output simulation parameters
  print_parameters();

  calc_g1();

  // open file
  f = fopen(filename, "w");

  if(f==NULL) {
    printf("Error: could not open '%s' for output\n", filename);
    return 1;
  }

  printf("Writing POLAR file...\n");

  // write output
  fprintf(f, "%d %d\n",nr,np);

  k = 0;
  for(i=0; i<nr; i++) {
    for(j=0; j<np; j++) {

      theta = 2.*M_PI*j/np;
      r = calc_r((double)i/(double)(nr-1.));
      sn = sin(theta);
      co = cos(theta);
      eps = r/a;
      
      x = r*(co - eps*eps*t*sn*sn) + s*(1.-eps*eps) + R;
      z = r*e*sn + z0;
      
      if(i==nr-1) {
	dx = a*(-sn - 2.*t*sn*co);
	dz = a*e*co;
	
	ddx =  a*(-co + 2.*t*(sn*sn - co*co));
	ddz = -a*e*sn;

	v = sqrt(dx*dx + dz*dz);
	
	curv = (dx*ddz - dz*ddx) / (v*v*v);

	norm1 =  dz/v;
	norm2 = -dx/v;
	
	fprintf(f, "%12.8f %12.8f %12.8f %12.8f %12.8f\n",
		x, z, norm1, norm2, curv);
      } else {
	fprintf(f, "%12.8f %12.8f\n", x, z);
      }
      k++;
    }
  }

  // close file
  fclose(f);

  printf("Done\n");

  return 0;
}

// f1 is the fraction of the maximum allowable g1 for packing
void calc_g1()
{
  int i;

  g1 = f1*w1/(1. - f1*w1*(atan((1.-s1)/w1) + atan(s1/w1)));
}

double calc_r(double f)
{
  double r;

  r = f;
  
  if(f1 != 0.) {
    r += g1*((f-1.)*atan(s1/w1) + f*atan((1.-s1)/w1) - atan((f-s1)/w1));
  }

  return a*sqrt(r);
}


void print_usage()
{
  printf("Usage: \n");
  printf("make_polar -R <R> -a <a> -s <s> -e <e> -t <t> -nr <nr> -np <np>\n");
  printf("           -s1 <s1> -f1 <f1> [input_file]\n");
  printf(" <R>: major radius (default=%g)\n", R0);
  printf(" <a>: minor radius (default=%g)\n", a0);
  printf(" <s>: Shafranov shift (default=%g)\n",s0);
  printf(" <e>: elongation (default=%g)\n", e0);
  printf(" <t>: triangularity (default=%g)\n", t0);
  printf(" <z0>: Z-coordinate of midplane (default=%g)\n", z00);
  printf(" <nr>: number of radial grid points (default=%d)\n", nr0);
  printf(" <np>: number of poloidal grid points (default=%d)\n", np0);
  printf(" <s1>: normalized position of packing region (0=axis, 1=boundary)\n");
  printf(" <w1>: normalized width of packing region\n");
  printf(" <f1>: packing amplitude (0=no packing, 1=infinite packing)\n");
}

void print_parameters()
{
  printf("Major radius:\t%g\n", R);
  printf("Minor radius:\t%g\n", a);
  printf("Shafranov shift:\t%g\n", s);
  printf("Elongation:\t%g\n", e);
  printf("Triangularity:\t%g\n", t);
  printf("Z-coordinate of midplane:\t%g\n", z0);
  printf("Radial grid points:\t%d\n", nr);
  printf("Poloidal grid points:\t%d\n", np);
  printf("Packing surface:\t%g\n", s1);
  printf("Packing amplitude:\t%g\n", f1);
  printf("Packing region width:\t%g\n", w1);
}

int read_inputs(int argc, char* argv[])
{
  int i;
  char* arg;

  for(i=1; i<argc; i++) {
    arg = argv[i];
	 
    if(0==strcmp(arg,"-h")){
      return 1; break;
    } else if(0==strcmp(arg, "-R")) {
      if(++i==argc) return 1;
      R = atof(argv[i]);
    } else if(0==strcmp(arg, "-a")) {
      if(++i==argc) return 1;
      a = atof(argv[i]);
    } else if(0==strcmp(arg, "-e")) {
      if(++i==argc) return 1;
      e = atof(argv[i]);
    } else if(0==strcmp(arg, "-s")) {
      if(++i==argc) return 1;
      s = atof(argv[i]);
    } else if(0==strcmp(arg, "-t")) {
      if(++i==argc) return 1;
      t = atof(argv[i]);
    } else if(0==strcmp(arg, "-z0")) {
      if(++i==argc) return 1;
      z0 = atof(argv[i]);
    } else if(0==strcmp(arg, "-s1")) {
      if(++i==argc) return 1;
      s1 = atof(argv[i]);
    } else if(0==strcmp(arg, "-f1")) {
      if(++i==argc) return 1;
      f1 = atof(argv[i]);
    } else if(0==strcmp(arg, "-w1")) {
      if(++i==argc) return 1;
      w1 = atof(argv[i]);
    } else if(0==strcmp(arg, "-np")) {
      if(++i==argc) return 1;
      np = atoi(argv[i]);
    } else if(0==strcmp(arg, "-nr")) {
      if(++i==argc) return 1;
      nr = atoi(argv[i]);
    } else {
      return read_input_file(arg);
    }
  }
  
  return 0;
}

int read_input_file(const char* infile)
{
  FILE* f;
  ssize_t nread;
  size_t len;
  char* line = NULL;
  char* tok = NULL;
  
  f = fopen(infile, "r");

  if(f==NULL) {
    printf("Error reading %s for input\n", infile);
    return 1;
  }

  while(nread = getline(&line, &len, f) != -1) {
    tok = strtok(line, " =\n!");

    if(0==strcmp(tok, "R")) {
      tok = strtok(NULL, " =\n!");
      R = atof(tok);
    } else if(0==strcmp(tok, "a")) {
      tok = strtok(NULL, " =\n!");
      a = atof(tok);
    } else if(0==strcmp(tok, "e")) {
      tok = strtok(NULL, " =\n!");
      e = atof(tok);
    } else if(0==strcmp(tok, "s")) {
      tok = strtok(NULL, " =\n!");
      s = atof(tok);
    } else if(0==strcmp(tok, "t")) {
      tok = strtok(NULL, " =\n!");
      t = atof(tok);
    } else if(0==strcmp(tok, "z0")) {
      tok = strtok(NULL, " =\n!");
      z0 = atof(tok);
    } else if(0==strcmp(tok, "s1")) {
      tok = strtok(NULL, " =\n!");
      s1 = atof(tok);
    } else if(0==strcmp(tok, "f1")) {
      tok = strtok(NULL, " =\n!");
      f1 = atof(tok);
    } else if(0==strcmp(tok, "w1")) {
      tok = strtok(NULL, " =\n!");
      w1 = atof(tok);
    } else if(0==strcmp(tok, "np")) {
      tok = strtok(NULL, " =\n!");
      np = atoi(tok);
    } else if(0==strcmp(tok, "nr")) {
      tok = strtok(NULL, " =\n!");
      nr = atoi(tok);
    } else {
      printf("Unrecognized option in input file: %s\n", tok);
    }
  }

  if(line) free(line);

  fclose(f);

  return 0;
}
