#include "lanczos.h"
#include <limits>
#define ZERO numeric_limits<int>::min()
#define SPIN_DOWN 0
#define SPIN_UP 1

template <typename T>
T delta(T a1, T a2)
{
T temp = T(0);

if(a1 == a2)
  temp = T(1);

return temp;
}

inline int sign(int i)
{

if(i>=0) return 1;
else return -1;

}

int pow(int m, int n)
{
int temp=1;
for(int i=0;i<n;i++)
   temp*=m;
return temp;
}

inline int bitset2int(bitset<32> mybit)
{
int temp=0;

for(int i=0;i<32;i++) 
  {
   if(mybit.test(i)) temp+=pow(2,i);
  }
return temp;
}

inline bitset<32> int2bitset(int n)
{
bitset<32> tempbit(0);
int s=0;

for(int i=0;i<32;i++){
    s = pow(2,i);
    if( s&n ) tempbit.set(i);
  }

return tempbit;
}


inline int ElectronCount(int site1 ,int site2, int state )
{
bitset<32> tempstate(state);
int temp=0;
int site=0;

if(site1>site2) {site=site1;site1=site2;site2=site1;}

for(int i=site1+1;i<site2;i++)
   if(tempstate.test(i)) temp++;

return temp;
}


int c_(int site, int spin, int state)
{

  int state_ = sign(state)*state;
  bitset<32> mybit(state_);
  int temp=0,sign_=1;

  if( state == ZERO) return ZERO;

    if(mybit.test(site)) 
     mybit.set(site,0);
    else
     return ZERO; 
  sign_ = -pow(-1,ElectronCount(0,site,state_))*sign(state);

return bitset2int(mybit)*sign_;
}
 
int c_dag(int site, int spin, int state)
{
 int state_ = sign(state)*state;
 bitset<32> mybit(state_);
 int temp=0,sign_=1;

if(state ==ZERO) return ZERO;

   if(mybit.test(site)) 
    return ZERO;
   else
    {
      mybit.set(site,1);
      temp = bitset2int(mybit);
    }

sign_ = -pow(-1,ElectronCount(0,site,state_))*sign(state);

return temp*sign_;
}

inline int ElectronCount(int site, int Nsites, int spin, int state )
{
bitset<32> tempstate(state);
int temp=0;


for(int i=0;i<site;i++) 
  {
    if(spin == SPIN_DOWN)
     {if(tempstate.test(i)) temp++;}
    
    if(spin == SPIN_UP)
     {if(tempstate.test(Nsites+i)) temp++;} 
  }
return temp;
}

template <typename T>
void Resize_refill(matrix<T>& H, int size, T t)
{
vector<T>* temp = new vector<T>(size,T(0));
H.v.clear();
for(int i=0;i<size;i++) H.v.push_back(*temp);
delete temp;
return;
}

template <typename T>
void Resize_refill(CEll<T>& H, int size, T t)
{
vector< vector<T> > temp(size,vector<T>(0, T(0)) );
H.data.clear();
H.idx.clear();
H.data.swap( temp );
vector< vector<int> > temp1(size,vector<int>(0, 0) );
H.idx.swap( temp1 );
return;
}

/*
for(int i=0;i<N;i++) v1[i]=v2[i];
return;
}

template <typename T>
void memcpy(vector< vector<T> >& v1, T *v2,int cols,int rows)
{
vector<T> temp(cols,T(0));
for(int i=0;i<rows;i++) 
   {
     for(int j=0;j<cols;j++)
        temp[i] = v2[i*rows+j];
     v1.push_back(temp); 
   }

return;
}

template <typename T>
void memcpy(vector<T>& v1, T *v2,int cols)
{

for(int i=0;i<cols;i++) 
   v1.push_back(v2[i]);

return;
}

*/
/*
template <typename T>
void init(T *v, T val, int N)
{
for(int i=0;i<N;i++) v[i]=val; 
return;
}
*/
