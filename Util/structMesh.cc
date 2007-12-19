#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>

using std::cout;
using std::ofstream;
using std::ostream;
using std::istream;
using std::endl;
using std::string;

/* 
 creates a "structured mesh" that can be used with SCOREC software
 it is structured in the fact that it has the same topological and
 geometric structure as a true structured mesh but does not take
 advantage of any of this structure
 
 note that the bounding box has minimum x and z coordinate 0 

 to be able to use this with M3D-C1, this file must be processed by
 <SCOREC SOFTWARE PATH>/mctk/Examples/PPPL/PPPL/test/DISCRETE/main struct.sms
 which will output a struct-dmg.sms and struct.dmg file which can be used to
 run M3D-C1.
*/

int main(int argc, char * argv[]) {
  string fileName("struct.sms");
  if(argc != 5) {
    cout << "Need to add the number of vertices in x and z directions and the bounding box size\n";
    cout << "./<executable name> <number of vertices in x direction> <number of vertices in z direction> <maximum x> <maximum z>\n";
    return 1;
  }
  int m = atoi(argv[1]); //note m here is x-direction
  int n = atoi(argv[2]); //note n here is y-direction
  cout << "There are " << m << " vertices in the x-dir and "
       << n << " vertices in the y-dir\n";

  FILE* outFile = fopen(fileName.data(), "w");      
  fprintf(outFile,"sms 2\n");
  fprintf(outFile, "%d %d %d %d %d\n", 
	  0, (m-1)*(n-1)*2,(m-1)*n+(n-1)*m+(m-1)*(n-1),m*n,m*n);  
  int i, j;
  int gentid, gentitytype, numEdges;
  double xDomainLength = atof(argv[3]);
  double yDomainLength = atof(argv[4]);
  cout << "The bounding box is x:0 - "<< xDomainLength <<
    ", y:0 - " << yDomainLength << endl;
  double dx = xDomainLength/(m-1.);
  double dy = yDomainLength/(n-1.);
  for(j=0;j<n;j++)  
    for(i=0;i<m;i++) 
      {
	if(i==0) {
	  if(j==0) {
	    gentid = 2;
	    gentitytype = 0;
	    numEdges = 3;
	  }
	  else if(j==(n-1)) {
	    gentid = 3;
	    gentitytype = 0;
	    numEdges = 2;
	  }
	  else {
	    gentid = 0;
	    gentitytype = 1;
	    numEdges = 4;
	  }
	}
	else if(i==(m-1)) {
	  if(j==0) {
	    gentid = 1;
	    gentitytype = 0;
	    numEdges = 2;
	  }
	  else if(j==(n-1)) {
	    gentid = 0;
	    gentitytype = 0;
	    numEdges = 3;
	  }
	  else {
	    gentid = 2;
	    gentitytype = 1;
	    numEdges = 4;
	  }
	}
	else if(j==0) {
	  gentid = 1;
	  gentitytype = 1;
	  numEdges = 4;
	}
	else if(j==(n-1)) {
	  gentid = 3;
	  gentitytype = 1;
	  numEdges = 4;
	}
	else {
	  gentid = 0;
	  gentitytype = 2;
	  numEdges = 6;
	}
	    
	fprintf(outFile, "%d %d %d\n%.15f %.15f 0 ",
	      gentid+1, gentitytype, numEdges, 
	      i*dx, j*dy);
	switch(gentitytype)
	  {
          case 1: //ofs << " 0";
	    fprintf(outFile, "0");
	    break;
          case 2: //ofs << " 0" << " 0" << " 0";
	    fprintf(outFile, "0 0 0");
	    break;
          default: break;
	  }
	fprintf(outFile, "\n");      
      } //done i,j vertex loop
  int numfaces;
  for(j=0;j<n;j++)  
    for(i=0;i<(m-1);i++) {
      if(j==0) {
	gentid = 1;gentitytype = 1;numfaces=1;
      }
      else if(j==(n-1)) {
	gentid = 3;gentitytype = 1;numfaces=1;
      }
      else {
	gentid = 0;gentitytype = 2;numfaces=2;
      }
      fprintf(outFile, "%d %d %d %d %d 0\n",  gentid+1, gentitytype, i+m*j+1, i+1+m*j+1, numfaces);
    }
  for(j=0;j<(n-1);j++)  
    for(i=0;i<m;i++) {
      if(i==0) {
	gentid = 0;gentitytype = 1;numfaces=1;
      }
      else if(i==(m-1)) {
	gentid = 2;gentitytype = 1;numfaces=1;
      }
      else {
	gentid = 0;gentitytype = 2;numfaces=2;
      }
      fprintf(outFile, "%d %d %d %d %d 0\n",  gentid+1, gentitytype, i+m*j+1, i+m+m*j+1, numfaces);
    }
  for(j=0;j<(n-1);j++)  
    for(i=0;i<(m-1);i++) 
      fprintf(outFile, "%d %d %d %d 2 0\n",  gentid+1, gentitytype, i+m*j+1, i+m+m*j+2);
  //done edges
  for(j=0;j<(n-1);j++)
    for(i=0;i<(m-1);i++) {
      fprintf(outFile, "1 2 3 %d %d %d 0\n", i+j*(m-1)+1, (m-1)*n+i+1+j*m+1,
	      -((m-1)*n+(n-1)*m+i+j*(m-1)+1));
      fprintf(outFile, "1 2 3 %d %d %d 0\n", -((m-1)*n+i+1+j*m), (m-1)*n+(n-1)*m+i+j*(m-1)+1,
	      -(i+j*(m-1)+1+m-1));
    }
  fclose(outFile);
  return 0;
}
