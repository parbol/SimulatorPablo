#include "Simulatorttbar.C"
#include <string>
#include <stdlib.h>

//To compile: g++ -o CallSimulator mt2_bisect.cpp CallSimulator.C `root-config --cflags --libs`


int main(int argc, char *argv[]) {

  if (argc != 11) {
      std::cout << "Usage: ./CallSimulatorttbar outputFileName.root topmass topwidth mintopmass maxtopmass wmass wmasswidth minwmass maxwmass numberOfEvents" << std::endl;
      return 1;
  }

  std::string outputname(argv[1]);
  float topmass = atof(argv[2]);
  float topwidth = atof(argv[3]);
  float mintop = atof(argv[4]);
  float maxtop = atof(argv[5]);
  float wmass = atof(argv[6]);
  float wwidth = atof(argv[7]);
  float minw = atof(argv[8]);
  float maxw = atof(argv[9]);
  int nevents = atoi(argv[10]);


  Simulatorttbar *a = new Simulatorttbar(nevents, topmass, topwidth, mintop, maxtop, wmass, wwidth, minw, maxw, outputname);

  delete a;

  return 0;


}

