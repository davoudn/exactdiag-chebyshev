#include "math.h"
#include "matrix.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <bitset>
#include <string>
#include <complex>

using namespace std;

// input block class
template <typename T>
class InputInformation
{
public: 
InputInformation();
~InputInformation();

void ReadInformation( string inputfilename );
int  ne_Up,ne_Down,nSites,omegaMesh;
T fBeta,h,h_prime,U,E,EL,V,Pulse_width,W;
int  nStates;
bool flag_randomize;

private:
fstream inputfilestream;
};

template <typename T>
InputInformation<T>::InputInformation():flag_randomize(false){};

template <typename T>
InputInformation<T>::~InputInformation(){};

template <typename T>
void InputInformation<T>::ReadInformation(string inputfilename)
{

inputfilestream.open(inputfilename.data(), fstream::in);
inputfilestream>>ne_Up>>ne_Down>>nSites;
inputfilestream>>h>>h_prime>>U>>V>>EL;
inputfilestream.close();
flag_randomize=false;
nStates = Tarkib(nSites,ne_Up)*Tarkib(nSites,ne_Down);
return;
}





