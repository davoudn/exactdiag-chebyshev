#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
//#include <math.h>
#include <string>
#include <bitset>
#include <complex>
#include <stdio.h>
//#include "ioclass.h"
//#include "green.h"
#include "exactdiag.h"
//#include <exactdiag.h>
//#include <stdexcept.h>
//#include <complex>
#define         F77_FUNC(name, NAME)   name ## _
# define PRECISION 2
#define IM std::complex<double>(0.0,1.0)
 
using namespace std;




int main(int argc , char* argv[] )
{

InputInformation<double> inputInformation;

if(argc<11){
  cout<<" There is no input information, or the information is incomplete!"<<endl;
  cout<<" The correct entery has the following format: ";
  cout<<" ./state.out N N_up N_down h h_prime U V W EL Pluse_width sim_type"<<endl;
  cout<<" /* ********************************************************************************************************************"<<endl;
  cout<<" Where N: total number of sites(integer), N_up: number of spin up particles <= N,  N_up: number of spin downticles <= N,"<<endl;
  cout<<" h, h_prime : alternate hoppings (float), U: local interaction(float), V:nearest neighbour interaction(float),"<<endl;
  cout<<" W:next nearest neighbour interaction(float), EL: electric field strength(float), Pulse_width: the duration of simulation(float)"<<endl;
  cout<<" sim_type: the type of the simulation(integer), now, sim_type=1 corresponds to propagation, sim_type=2 corresponds to propagation with"<<endl;
  cout<<" randomized hopping parameters by a ranrom portion of h_prime."<<endl; 
  cout<<" /* ********************************************************************************************************************"<<endl;
}
else{
  inputInformation.nSites  = atoi(*(argv+1));
  inputInformation.ne_Up   = atoi(*(argv+2));
  inputInformation.ne_Down = atoi(*(argv+3));
  inputInformation.h       = atof(*(argv+4));
  inputInformation.h_prime = atof(*(argv+5));
  inputInformation.U       = atof(*(argv+6));
  inputInformation.V       = atof(*(argv+7));
  inputInformation.W = atof(*(argv+8));
  inputInformation.EL      = atof(*(argv+9));
  inputInformation.Pulse_width = atof(*(argv+10));

//  cout<<"ne_Up="<<inputInformation.ne_Up<<"ne_Down="<<inputInformation.ne_Down<<"  nSite="<<inputInformation.nSites;
//  cout<<"  h="<<inputInformation.h<<" h_prime="<<inputInformation.h_prime<<" U="<<inputInformation.U<<" Electric field="<<inputInformation.EL<<endl;

  if(atoi( *(argv+11) ) == 1){
    Cluster<double, complex<double> > cluster(inputInformation);
    cluster.TimeEvolution();
  }

  if(atoi( *(argv+11) ) == 2){
    inputInformation.flag_randomize=true;  
    Cluster<double, complex<double> > cluster(inputInformation);
    cluster.TimeEvolution();
  }
}
return 0;
}

