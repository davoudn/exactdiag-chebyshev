#include "green.h"
#include "aux.h"
#include "math.h"
#include <bitset>
#include <map>
#include <string>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* Time */


using namespace std;

/*            Cluster class
**************************************
*/

template <typename T>
void Clear(vector< vector <T> >& v)
{
vector< vector<T> > temp(0,vector<T>(0,T(0)));
v.swap(temp);
return;
}

template <typename T>
void Clear(vector <T>& v)
{
vector<T> temp(0,T(0));
v.swap(temp);
return;
}

int mod(int i, int j)
{

return i - (i/j)*j;
}

/*                               cluster class                                   */
template<typename T , typename C>
class Cluster
{
public:
 Cluster(InputInformation<T>& inputInformation);
 ~Cluster();
 C    Gimp(int m,int m_, C om);
 void Gimp(wMesh<T,C> om, Gfunct<C>& G_imp);
 void SPDiagonalize(T Time_in);
 void Lapack_diag();
 void TimeEvolution();
 string monitorfilename,groundstatefilename;

private:
 int nStates, nSites, nTotSites, dim_up, dim_down, Maxiter;
 int ne_Up,ne_Down,Tarkib_,MaxTimesteps;
 T   t, V, U, W, h, h_prime, Time, Time_step, EL, pi, E0, Charge_gap, Pulse_width;
 matrix<C> h_rand, hOpings,J_ij; 
 vector<C> U_rand;
 vector< map<int,int> > Spinconfig;
  //lanczos vectors
 vector<C> Groundvecstate_0, initstate;
 CEll<C> SPHamiltonian, J, K, T_up, T_down, J_up, J_down;
 matrix<C> Hamiltonian;

 private: /* Private functions */

 inline int IsHope(int idx,int idx0, int state){
         bitset<32> temp(state);
         int temp1=0;

         if(temp.test(idx)==true && temp.test(idx0)==false){
           temp[idx]=0; temp[idx0]=1;
           int n = ElectronCount(idx,idx0,state);
           temp1 = -bitset2int(temp)*pow(T(-1),n);
         }

      return temp1;
      }

 inline int Mindex(int idx_up,int idx_down){
      return dim_up*idx_down+idx_up;
      }

 inline int IndexExtract(int state,int spin){
      int temp = -1;
      map<int,int>::iterator it;
      it=Spinconfig[spin].find(state);
      if( it !=Spinconfig[spin].end())
        temp = it->second;
      return temp;
 }

 void Setup_states();
 void MElements_best();
 void InitHopings_J(); 
 void MElements_nonsparse();
 vector<C> vec_c_(int site, int spin, vector<C>& vec_state);
 vector<C> vec_c_dag(int site, int spin, vector<C>& vec_state);

  // fractional grean functions
 vector<C> a_lesser, a_greater, a_0;
 vector<C> b_lesser, b_greater, b_0;
 
 void Current_melements(vector<C> &vec_left_, vector<C> &vec_in_, T* current, T* kinetic_energy);
 bool flag_randomize,lapack_evolution;
};

template<typename T,typename C>
Cluster<T,C>::Cluster(InputInformation<T>&  inputInformation):hOpings(inputInformation.nSites,inputInformation.nSites,C(0)),J_ij(inputInformation.nSites,inputInformation.nSites,C(0)),h_rand(inputInformation.nSites,inputInformation.nSites,C(0)),U_rand(inputInformation.nSites,C(0)), Spinconfig(2,map<int,int>())
{
char temp_str[70];
 
nStates = inputInformation.nStates;
nSites  = inputInformation.nSites;
ne_Up   = inputInformation.ne_Up;
ne_Down = inputInformation.ne_Down;
h = inputInformation.h;
h_prime = inputInformation.h_prime;
U = inputInformation.U;
V = inputInformation.V;
W = inputInformation.W;
EL = inputInformation.EL;
nTotSites = nSites;
Maxiter = 5000;
Time_step=0.001/EL;
MaxTimesteps = 10000;
pi = acos(-1);
Time = 0.0;
lapack_evolution=false;
Pulse_width= inputInformation.Pulse_width;
flag_randomize = inputInformation.flag_randomize;

if(flag_randomize == true){
  srand (0);
  T temp=T(0);
  for(int i=0;i<nSites-1;i++){
     temp = (T) (rand()%100)*pow(100.0,-1.0) ;
     cout << "pp"<<temp<<endl;
     h_rand[i][i+1] =  (h+h_prime*temp);
     h_rand[i+1][i] =  (h+h_prime*temp);
     h_rand[i][i] =  h_prime*temp;
     U_rand[i]= U+h_prime*temp;
  }
  h_rand[nSites-1][0] =  (h+h_prime*temp);
  h_rand[0][nSites-1] =  (h+h_prime*temp);
}

if(flag_randomize)
  sprintf( temp_str, "monitor-%d-%d-%d-%1.2f-%1.2f-%1.2f-%1.2f-%1.4f-%1.4f.out",inputInformation.nSites,inputInformation.ne_Up,inputInformation.ne_Down,inputInformation.h,inputInformation.h_prime,inputInformation.U,inputInformation.V,inputInformation.W,inputInformation.EL);
else
  sprintf( temp_str, "fin_rand_monitor-%d-%d-%d-%1.2f-%1.2f-%1.2f-%1.2f-%1.4f-%1.4f.out",inputInformation.nSites,inputInformation.ne_Up,inputInformation.ne_Down,inputInformation.h,inputInformation.h_prime,inputInformation.U,inputInformation.V,inputInformation.W,inputInformation.EL);

monitorfilename=temp_str;

cout<<" >>> --------------------- Parameters of Simulation --------------------<<<"<<endl;
cout<<" ==> Step number= 0"<<" ,Time="<<Time<<" ,Timestep="<<Time_step<<endl;
cout<<" ==> N="<<nSites<<", N_up="<<ne_Up<<", N_down="<<ne_Down<<endl;
cout<<" ==> U="<<U<<", V="<<V<<", W="<<W<<" ,h="<<h<<" ,h_prime="<<h_prime<<endl;
cout<<" ==> EL="<<EL<<", Pulse_width="<<Pulse_width<<endl;
cout<<" >>> -------------------------------------------------------------------<<<"<<endl;
cout<<endl;

};

template <typename T , typename C>
void Cluster<T,C>::InitHopings_J()
{

if( flag_randomize == true )
{
  for(int i=0;i<nSites-1;i++){
     hOpings[i][i+1] =  h_rand[i][i+1]*exp(-IM*pi*C(2)*(Time*EL));
     hOpings[i+1][i] =  h_rand[i+1][i]*exp( IM*pi*C(2)*(Time*EL));

     J_ij[i][i+1] =  IM*h_rand[i][i+1]*exp(-IM*pi*C(2)*(Time*EL));
     J_ij[i+1][i] =  -IM*h_rand[i+1][i]*exp( IM*pi*C(2)*(Time*EL));


   hOpings[nSites-1][0] =  h_rand[nSites-1][0]*exp(-IM*pi*C(2)*(Time*EL));
   hOpings[0][nSites-1] =  h_rand[0][nSites-1]*exp(IM*pi*C(2)*(Time*EL));

   J_ij[nSites-1][0] =  IM*h_rand[nSites-1][0]*exp(-IM*pi*C(2)*(Time*EL));
   J_ij[0][nSites-1] =  -IM*h_rand[0][nSites-1]*exp(IM*pi*C(2)*(Time*EL));
  }
}
else{
  for(int i=0;i<nSites-1;i++){
      if( mod(i+1,2) == 1  ){
          hOpings[i][i+1] =  h*exp(-IM*pi*C(2)*(Time*EL));
          hOpings[i+1][i] =  h*exp( IM*pi*C(2)*(Time*EL));

          J_ij[i][i+1] =  IM*h*exp(-IM*pi*C(2)*(Time*EL));
          J_ij[i+1][i] =  -IM*h*exp( IM*pi*C(2)*(Time*EL));

        }
      else{
          hOpings[i][i+1] =  h_prime*exp(-IM*pi*C(2)*(Time*EL));
          hOpings[i+1][i] =  h_prime*exp( IM*pi*C(2)*(Time*EL));

          J_ij[i][i+1] =  IM*h_prime*exp(-IM*pi*C(2)*(Time*EL));
          J_ij[i+1][i] =  -IM*h_prime*exp( IM*pi*C(2)*(Time*EL));

        }
  }
//if(flag_periodic == true){
   hOpings[nSites-1][0] =  h_prime*exp(-IM*pi*C(2)*(Time*EL));
   hOpings[0][nSites-1] =  h_prime*exp(IM*pi*C(2)*(Time*EL));

   J_ij[nSites-1][0] =  IM*h_prime*exp(-IM*pi*C(2)*(Time*EL));
   J_ij[0][nSites-1] =  -IM*h_prime*exp(IM*pi*C(2)*(Time*EL));
}

  vector<T> Er(hOpings.cols(),T(0));
  vector<C> z(hOpings.cols()*hOpings.cols(),C(0));
  T temp=T(0);

zheevr(hOpings , Er, z, hOpings.cols());

for(int i=0;i<ne_Up;i++) temp+= Er[i];
for(int i=0;i<ne_Down;i++) temp+= Er[i];
cout.precision(15);
cout<<endl;
cout<<" ==> Hopings ground state= "<<temp<<endl;
return;
}

template<typename T,typename C>
void Cluster<T,C>::Setup_states()
{
 bitset<32> temp_spinup(0),temp_spindown(0);
 int configcount = 0;
 dim_down=0;dim_up=0;
 Spinconfig[0].clear();
 Spinconfig[1].clear();

 Tarkib_=pow(T(2),nTotSites);
  
   for(int i=0;i<pow(T(2),nTotSites);i++){
      temp_spinup = int2bitset(i);
      if(temp_spinup.count() == ne_Up){ 
        Spinconfig[SPIN_UP].insert(pair<int,int>(i,dim_up));   
        dim_up++; 
      }
    }

   for(int i=0;i<pow(T(2),nTotSites);i++){
      temp_spindown = int2bitset(i);
      if(temp_spindown.count() == ne_Down){
        Spinconfig[SPIN_DOWN].insert(pair<int,int>(i,dim_down));
        dim_down++;
      }
    }
  
 dim_up = Spinconfig[SPIN_UP].size();
 dim_down = Spinconfig[SPIN_DOWN].size();
 nStates = dim_up*dim_down;
 cout<<" ==> Setting up the Hilbert space: ";
 cout<<"nStates="<<nStates<<", dim_up="<<dim_up<<", dim_down="<<dim_down<<endl;
 cout<<endl;
 Maxiter = 100;
 return;
}

template<typename T,typename C>
void Cluster<T,C>::MElements_best()
{

int temp1=0,temp2=0,temp_down=0,temp_up=0,temp_idx=0,temp_=0;
C tempinteraction=C(0), temp_h=C(0),temp_J=C(0);
bitset<32> tempbit_up,tempbit_down,tempbit_overlap;
// forming spin up and down hoping blocks
Resize_refill(T_up ,dim_up, C(0));
Resize_refill(T_down ,dim_down, C(0));
Resize_refill(SPHamiltonian ,dim_up*dim_down, C(0));

Resize_refill(J_up ,dim_up, C(0));
Resize_refill(J_down ,dim_down, C(0));

for(map<int,int>::iterator idx_down=Spinconfig[SPIN_DOWN].begin();idx_down!=Spinconfig[SPIN_DOWN].end();++idx_down)
   for(int idx=0;idx<nTotSites;idx++)
      for(int idx0=idx+1;idx0<nTotSites;idx0++)
         if(hOpings[idx][idx0]!=0.0){ 
             temp_down = IsHope(idx,idx0,idx_down->first);
             if(temp_down!=0){
                temp_idx = IndexExtract(abs(temp_down),SPIN_DOWN);
                temp_h = (C)sign(temp_down)*hOpings[idx][idx0];  
                temp_J=(C)sign(temp_down)*J_ij[idx][idx0];
                T_down.accumulate(idx_down->second,temp_idx,temp_h);
                J_down.accumulate(idx_down->second,temp_idx,temp_J);
             }
         }            
              
for(map<int,int>::iterator idx_up=Spinconfig[SPIN_UP].begin();idx_up!=Spinconfig[SPIN_UP].end();++idx_up)
   for(int idx=0;idx<nTotSites;idx++)
      for(int idx0=idx+1;idx0<nTotSites;idx0++)
         if(hOpings[idx][idx0]!=0.0){
             temp_up = IsHope(idx,idx0,idx_up->first);
             if(temp_up!=0){ 
                 temp_idx = IndexExtract(abs(temp_up),SPIN_UP);
                 temp_h = (C)sign(temp_up)*hOpings[idx][idx0];  
                 temp_J = (C)sign(temp_up)*J_ij[idx][idx0];
                 T_up.accumulate(idx_up->second,temp_idx,temp_h);
                 J_up.accumulate(idx_up->second,temp_idx,temp_J);
             }
         }
// diagonal terms 
for(map<int,int>::iterator idx_down=Spinconfig[SPIN_DOWN].begin();idx_down!=Spinconfig[SPIN_DOWN].end();++idx_down)
   for(map<int,int>::iterator idx_up=Spinconfig[SPIN_UP].begin();idx_up!=Spinconfig[SPIN_UP].end();++idx_up){
       // hoping termsrp
      temp1 = Mindex(idx_up->second,idx_down->second);               
      for(int i=0;i<T_down.idx[idx_down->second].size();i++){     
         temp2 = Mindex(idx_up->second,T_down.idx[idx_down->second][i]);           
         SPHamiltonian.accumulate(temp1,temp2,T_down.data[idx_down->second][i]);
         SPHamiltonian.accumulate(temp2,temp1,conj(T_down.data[idx_down->second][i]));
      }

      for(int i=0;i<T_up.idx[idx_up->second].size();i++){
         temp2 = Mindex(T_up.idx[idx_up->second][i], idx_down->second);
         SPHamiltonian.accumulate(temp1,temp2,T_up.data[idx_up->second][i]);
         SPHamiltonian.accumulate(temp2,temp1,conj(T_up.data[idx_up->second][i]));
      }

      tempbit_down = int2bitset(idx_down->first);
      for(int i=0;i<nTotSites;i++) 
         tempinteraction+=(C)tempbit_down[i]*hOpings[i][i];

      tempbit_up = int2bitset(idx_up->first);
      for(int i=0;i<nTotSites;i++) 
         tempinteraction+=(C)tempbit_up[i]*hOpings[i][i];

      tempbit_overlap = tempbit_down&tempbit_up;
// nearest neighbor columb interaction
        
      for(int i=1;i<nSites-1;i++){
         temp_+= tempbit_down[i]*(tempbit_down[i+1]+tempbit_down[i-1]);
         temp_+= tempbit_up[i]*(tempbit_up[i+1]+tempbit_up[i-1]);

         temp_+= tempbit_down[i]*(tempbit_up[i+1]+tempbit_up[i-1]);
         temp_+= tempbit_up[i]*(tempbit_down[i+1]+tempbit_down[i-1]);
      }
//
      temp_+= tempbit_down[0]*(tempbit_down[1]+tempbit_down[nSites-1]);
      temp_+= tempbit_up[0]*(tempbit_up[1]+tempbit_up[nSites-1]);
      temp_+= tempbit_down[0]*(tempbit_up[1]+tempbit_up[nSites-1]);
      temp_+= tempbit_up[0]*(tempbit_down[1]+tempbit_down[nSites-1]);
//
      temp_+= tempbit_down[nSites-1]*(tempbit_down[0]+tempbit_down[nSites-2]);
      temp_+= tempbit_up[nSites-1]*(tempbit_up[0]+tempbit_up[nSites-2]);
      temp_+= tempbit_down[nSites-1]*(tempbit_up[0]+tempbit_up[nSites-2]);
      temp_+= tempbit_up[nSites-1]*(tempbit_down[0]+tempbit_down[nSites-2]);       
      tempinteraction += (temp_/2)*V;
      temp_=0;
// next nearest neighbour interactopn

      for(int i=2;i<nSites-2;i++){
         temp_+= tempbit_down[i]*(tempbit_down[i+2]+tempbit_down[i-2]);
         temp_+= tempbit_up[i]*(tempbit_up[i+2]+tempbit_up[i-2]);

         temp_+= tempbit_down[i]*(tempbit_up[i+2]+tempbit_up[i-2]);
         temp_+= tempbit_up[i]*(tempbit_down[i+2]+tempbit_down[i-2]);
      }
//
      temp_+= tempbit_down[0]*(tempbit_down[2]+tempbit_down[nSites-2]);
      temp_+= tempbit_up[0]*(tempbit_up[2]+tempbit_up[nSites-2]);
      temp_+= tempbit_down[0]*(tempbit_up[2]+tempbit_up[nSites-2]);
      temp_+= tempbit_up[0]*(tempbit_down[2]+tempbit_down[nSites-2]);

//
      temp_+= tempbit_down[1]*(tempbit_down[3]+tempbit_down[nSites-1]);
      temp_+= tempbit_up[1]*(tempbit_up[3]+tempbit_up[nSites-1]);
      temp_+= tempbit_down[1]*(tempbit_up[3]+tempbit_up[nSites-1]);
      temp_+= tempbit_up[1]*(tempbit_down[3]+tempbit_down[nSites-1]);

//
      temp_+= tempbit_down[nSites-1]*(tempbit_down[1]+tempbit_down[nSites-3]);
      temp_+= tempbit_up[nSites-1]*(tempbit_up[1]+tempbit_up[nSites-3]);
      temp_+= tempbit_down[nSites-1]*(tempbit_up[1]+tempbit_up[nSites-3]);
      temp_+= tempbit_up[nSites-1]*(tempbit_down[1]+tempbit_down[nSites-3]);

// 
      temp_+= tempbit_down[nSites-2]*(tempbit_down[0]+tempbit_down[nSites-4]);
      temp_+= tempbit_up[nSites-2]*(tempbit_up[0]+tempbit_up[nSites-4]);
      temp_+= tempbit_down[nSites-2]*(tempbit_up[0]+tempbit_up[nSites-4]);
      temp_+= tempbit_up[nSites-2]*(tempbit_down[0]+tempbit_down[nSites-4]);

      tempinteraction += (temp_/2)*W;
      temp_=0.0;

      tempinteraction += U*(T)tempbit_overlap.count();
      SPHamiltonian.accumulate(temp1,temp1,tempinteraction);
      tempinteraction = C(0); 
   }
return;
}



template<typename T,typename C>
void Cluster<T,C>::MElements_nonsparse()
{

int temp1=0,temp2=0,temp_down=0,temp_up=0,temp_idx=0,temp_=0;
C tempinteraction=C(0), temp_h=C(0),temp_J=C(0);
bitset<32> tempbit_up,tempbit_down,tempbit_overlap;

// forming spin up and down hoping blocks
Resize_refill(T_up ,dim_up, C(0));
Resize_refill(T_down ,dim_down, C(0));
Resize_refill(Hamiltonian ,dim_up*dim_down, C(0));

Resize_refill(J_up ,dim_up, C(0));
Resize_refill(J_down ,dim_down, C(0));

for(map<int,int>::iterator idx_down=Spinconfig[SPIN_DOWN].begin();idx_down!=Spinconfig[SPIN_DOWN].end();++idx_down)
   for(int idx=0;idx<nTotSites;idx++)
      for(int idx0=idx+1;idx0<nTotSites;idx0++)
         if(hOpings[idx][idx0]!=0.0){ 
             temp_down = IsHope(idx,idx0,idx_down->first);
             if(temp_down!=0){
                temp_idx = IndexExtract(abs(temp_down),SPIN_DOWN);
                temp_h = (C)sign(temp_down)*hOpings[idx][idx0];  
                temp_J=(C)sign(temp_down)*J_ij[idx][idx0];
                T_down.accumulate(idx_down->second,temp_idx,temp_h);
                J_down.accumulate(idx_down->second,temp_idx,temp_J);
             }
         }            
              
for(map<int,int>::iterator idx_up=Spinconfig[SPIN_UP].begin();idx_up!=Spinconfig[SPIN_UP].end();++idx_up)
   for(int idx=0;idx<nTotSites;idx++)
      for(int idx0=idx+1;idx0<nTotSites;idx0++)
         if(hOpings[idx][idx0]!=0.0){
             temp_up = IsHope(idx,idx0,idx_up->first);
             if(temp_up!=0){ 
                 temp_idx = IndexExtract(abs(temp_up),SPIN_UP);
                 temp_h = (C)sign(temp_up)*hOpings[idx][idx0];  
                 temp_J = (C)sign(temp_up)*J_ij[idx][idx0];
                 T_up.accumulate(idx_up->second,temp_idx,temp_h);
                 J_up.accumulate(idx_up->second,temp_idx,temp_J);
             }
         }
           
// forming the hamitonian in full hilbert space
// second interactions ,interactions:
for(map<int,int>::iterator idx_down=Spinconfig[SPIN_DOWN].begin();idx_down!=Spinconfig[SPIN_DOWN].end();++idx_down)
   for(map<int,int>::iterator idx_up=Spinconfig[SPIN_UP].begin();idx_up!=Spinconfig[SPIN_UP].end();++idx_up){
       // hoping termsrp
      temp1 = Mindex(idx_up->second,idx_down->second);               
      for(int i=0;i<T_down.idx[idx_down->second].size();i++){     
         temp2 = Mindex(idx_up->second,T_down.idx[idx_down->second][i]);            
         Hamiltonian[temp1][temp2]+=T_down.data[idx_down->second][i];
         Hamiltonian[temp2][temp1]+=conj(T_down.data[idx_down->second][i]);
      }

      for(int i=0;i<T_up.idx[idx_up->second].size();i++){
         temp2 = Mindex(T_up.idx[idx_up->second][i], idx_down->second);
         Hamiltonian[temp1][temp2]+=T_up.data[idx_up->second][i];
         Hamiltonian[temp2][temp1]+=conj(T_up.data[idx_up->second][i]);
      }

// diagonal terms
      tempbit_down = int2bitset(idx_down->first);
      for(int i=0;i<nTotSites;i++) 
         tempinteraction+=(C)tempbit_down[i]*hOpings[i][i];

      tempbit_up = int2bitset(idx_up->first);
      for(int i=0;i<nTotSites;i++) 
         tempinteraction+=(C)tempbit_up[i]*hOpings[i][i];

      tempbit_overlap = tempbit_down&tempbit_up;
// nearest neighbor columb interaction
        
      for(int i=1;i<nSites-1;i++){
         temp_+= tempbit_down[i]*(tempbit_down[i+1]+tempbit_down[i-1]);
         temp_+= tempbit_up[i]*(tempbit_up[i+1]+tempbit_up[i-1]);

         temp_+= tempbit_down[i]*(tempbit_up[i+1]+tempbit_up[i-1]);
         temp_+= tempbit_up[i]*(tempbit_down[i+1]+tempbit_down[i-1]);
      }
//
      temp_+= tempbit_down[0]*(tempbit_down[1]+tempbit_down[nSites-1]);
      temp_+= tempbit_up[0]*(tempbit_up[1]+tempbit_up[nSites-1]);
      temp_+= tempbit_down[0]*(tempbit_up[1]+tempbit_up[nSites-1]);
      temp_+= tempbit_up[0]*(tempbit_down[1]+tempbit_down[nSites-1]);
//
      temp_+= tempbit_down[nSites-1]*(tempbit_down[0]+tempbit_down[nSites-2]);
      temp_+= tempbit_up[nSites-1]*(tempbit_up[0]+tempbit_up[nSites-2]);
      temp_+= tempbit_down[nSites-1]*(tempbit_up[0]+tempbit_up[nSites-2]);
      temp_+= tempbit_up[nSites-1]*(tempbit_down[0]+tempbit_down[nSites-2]);       
      tempinteraction += (temp_/2)*V;
      temp_=0;
// next nearest neighbour interactopn

      for(int i=2;i<nSites-2;i++){
         temp_+= tempbit_down[i]*(tempbit_down[i+2]+tempbit_down[i-2]);
         temp_+= tempbit_up[i]*(tempbit_up[i+2]+tempbit_up[i-2]);

         temp_+= tempbit_down[i]*(tempbit_up[i+2]+tempbit_up[i-2]);
         temp_+= tempbit_up[i]*(tempbit_down[i+2]+tempbit_down[i-2]);
      }
//
      temp_+= tempbit_down[0]*(tempbit_down[2]+tempbit_down[nSites-2]);
      temp_+= tempbit_up[0]*(tempbit_up[2]+tempbit_up[nSites-2]);
      temp_+= tempbit_down[0]*(tempbit_up[2]+tempbit_up[nSites-2]);
      temp_+= tempbit_up[0]*(tempbit_down[2]+tempbit_down[nSites-2]);

//
      temp_+= tempbit_down[1]*(tempbit_down[3]+tempbit_down[nSites-1]);
      temp_+= tempbit_up[1]*(tempbit_up[3]+tempbit_up[nSites-1]);
      temp_+= tempbit_down[1]*(tempbit_up[3]+tempbit_up[nSites-1]);
      temp_+= tempbit_up[1]*(tempbit_down[3]+tempbit_down[nSites-1]);

//
      temp_+= tempbit_down[nSites-1]*(tempbit_down[1]+tempbit_down[nSites-3]);
      temp_+= tempbit_up[nSites-1]*(tempbit_up[1]+tempbit_up[nSites-3]);
      temp_+= tempbit_down[nSites-1]*(tempbit_up[1]+tempbit_up[nSites-3]);
      temp_+= tempbit_up[nSites-1]*(tempbit_down[1]+tempbit_down[nSites-3]);

// 
      temp_+= tempbit_down[nSites-2]*(tempbit_down[0]+tempbit_down[nSites-4]);
      temp_+= tempbit_up[nSites-2]*(tempbit_up[0]+tempbit_up[nSites-4]);
      temp_+= tempbit_down[nSites-2]*(tempbit_up[0]+tempbit_up[nSites-4]);
      temp_+= tempbit_up[nSites-2]*(tempbit_down[0]+tempbit_down[nSites-4]);

      tempinteraction += (temp_/2)*W;
      temp_=0.0;

      tempinteraction += U*(T)tempbit_overlap.count();
      Hamiltonian[temp1][temp1]+=tempinteraction;
      tempinteraction = C(0); 
   }

return;
}


template <typename T,typename C>
void Cluster<T,C>::SPDiagonalize(T Time_in)
{
Time=Time_in;
fstream is;
T E_old=T(0),E_min=T(0),E_max=T(0);
vector<C> *temp;
int dim=0;
T E_=T(0);
Charge_gap=T(0);

//calculating  ground state of normal system(ne_Tot)
   Setup_states();
   Resize_refill(SPHamiltonian,nStates,C(0));
   InitHopings_J();
   MElements_best();

   Groundvecstate_0.clear();
   temp = new vector<C>(nStates,C(1,0));
   Groundvecstate_0.swap( *temp );
   delete temp;

   dim = SPHamiltonian.idx.size();
   Lanczos1(SPHamiltonian, dim, Maxiter,a_0, b_0, Groundvecstate_0,&E_min,&E_max,1);
   E0 = E_min;
   Charge_gap -=2.0*E0;
   cout<<" calculation of ground states finished"<<endl;

//calculating the a and b coeficiens for ne_To, the initial vector is calculated as c_dag_{i\sigma}|Groundstate>
   initstate = vec_c_(0, SPIN_UP, Groundvecstate_0);
   Resize_refill(SPHamiltonian, nStates, C(0));
   MElements_best();
   dim = SPHamiltonian.idx.size();

   Lanczos1(SPHamiltonian, dim, Maxiter ,a_lesser, b_lesser, initstate,&E_min,&E_max,0);
   Charge_gap +=E_min;

   initstate.clear();
   cout<<" calculation of system with one electron less finished"<<endl;
//coming back to original configuration state 
   ne_Up+=1;
   Setup_states();

//calculating the a and b coeficiens for ne_Tot-1, the initial vector is calculated as c_{i\sigma}|Groundstate>
   initstate = vec_c_dag(0,SPIN_UP,Groundvecstate_0);
   Resize_refill(SPHamiltonian,nStates,C(0));
   MElements_best();
   dim = SPHamiltonian.idx.size(); 
   Lanczos1(SPHamiltonian, dim, Maxiter ,a_greater, b_greater, initstate, &E_min,&E_max,0);
   Charge_gap +=E_min;   

   initstate.clear();
   cout<<" calculation of system with one electron more finished"<<endl;

// again comming back to the original configuration state 
   ne_Up-=1;
   Setup_states();
   if(Time_in==T(0))
     is.open(groundstatefilename.data(),fstream::out);
   else
     is.open(groundstatefilename.data(), std::fstream::out | std::fstream::app);
   is.precision(6);
   is<<fixed;
   cout<<"filename to write:"<<groundstatefilename.data()<<endl;

   is<<Time<<" "<<Charge_gap<<endl;
   is.close();
   cout<<"charge-gap="<<Charge_gap<<endl;

return;
}


template <typename T,typename C>
void Cluster<T,C>::TimeEvolution()
{
   T p = T(0), E_min=T(0), E_max=T(0), current=T(0), kinetic_energy=T(0), Time_wave=T(0), energy=T(0);
   int counter=0,temp_site=0;
   fstream is,is_wave,is_prob;
   int res=0;
   
   Time=T(0);
   cout<<">>>----------------------- Starting of time evolution --------------------------<<<"<<endl;
   Setup_states();
   InitHopings_J();
   Resize_refill(SPHamiltonian,nStates,C(0));
   MElements_best();

   Groundvecstate_0.clear();
   vector<C>*temp = new vector<C>(nStates,C(1));
   Groundvecstate_0.swap(*temp);
   delete temp;

   int dim = SPHamiltonian.idx.size();
   Lanczos1(SPHamiltonian, dim, 100, a_0, b_0, Groundvecstate_0, &E_min, &E_max,1);

   E_min -= 1.0;
   E_max += 1.0;
   E0 = E_min;  
  
   vector<C> Groundvecstate_cheb_(Groundvecstate_0);
   vector<C> psi_old(Groundvecstate_0);
   vector<C> psi_lanczos(Groundvecstate_0);

   C temp1 = C(0);
   static T A = (E_max-E_min)/2.0, probability=T(0);
   static T B = (E_max+E_min)/2.0;
   
   if (lapack_evolution == true){
      monitorfilename+=".lapack";
   }   

   is.open(monitorfilename.data(), fstream::out);
   is.precision(10);
   is<<fixed;

 
   while( Time < Pulse_width ){
      cout<<endl;
      counter++;
      Time += Time_step;
      Time_wave += Time_step;
      cout<<" ************************* Step number:"<<counter<<" ***************************"<<endl;
      cout<<" >>> --------------------- Parameters of Simulation --------------------<<<"<<endl;
      cout<<" ==> Time="<<Time<<" ,Timestep="<<Time_step<<endl;
      cout<<" ==> N="<<nSites<<", N_up="<<ne_Up<<", N_down="<<ne_Down<<endl;
      cout<<" ==> U="<<U<<", V="<<V<<", W="<<W<<" ,h="<<h<<" ,h_prime="<<h_prime<<endl;
      cout<<" ==> EL="<<EL<<", Pulse_width="<<Pulse_width<<endl;
      cout<<" >>> -------------------------------------------------------------------<<<"<<endl;

      if(lapack_evolution == false){
        InitHopings_J();
        MElements_best();
        Hscale(SPHamiltonian,E_min,E_max);
        Chebychev_step(SPHamiltonian, Groundvecstate_cheb_, Time_step,1.0e-16,E_min,E_max);
      }
      else{
        InitHopings_J();
        MElements_nonsparse();
        Lapack_step(Hamiltonian, Groundvecstate_cheb_, Time_step, 400);
      }

      energy = abs(Groundvecstate_cheb_*(SPHamiltonian*Groundvecstate_cheb_));
      psi_old= Groundvecstate_cheb_;

      energy = real(Groundvecstate_cheb_*(SPHamiltonian*Groundvecstate_cheb_));
      Resize_refill(SPHamiltonian,nStates,C(0));
      Current_melements(Groundvecstate_cheb_,Groundvecstate_cheb_, &current, &kinetic_energy);
       
      probability=abs(Groundvecstate_0*Groundvecstate_cheb_);

// ******* calculating order parameters and storing them ********  
     cout<<" ==>>"<<Time<<" "<<current<<" "<<kinetic_energy<<" "<<energy<<" "<<probability<<endl; 
     cout<<endl;
     is<<Time<<" "<<-current/((double)nSites)<<" "<<kinetic_energy/((double)nSites)<<" "<<energy/((double)nSites)<<" "<<probability<<endl;
   }

return;
}

template <typename T,typename C>
void Cluster<T,C>::Current_melements(vector<C> &vec_left_, vector<C> &vec_in_, T* current, T* kinetic_energy)
{

int temp1=0 ,temp2=0 ,temp_down=0 ,temp_up=0;
T tempinteraction=T(0);
bitset<32> tempbit_up, tempbit_down, tempbit_overlap;

Resize_refill(J, dim_up*dim_down, C(0));

for(int idx_down=0;idx_down<Spinconfig[SPIN_UP].size();idx_down++)
   for(int idx_up=0;idx_up<Spinconfig[SPIN_DOWN].size();idx_up++){
      temp1 = Mindex(idx_up,idx_down);              
      for(int i=0;i<T_down.idx[idx_down].size();i++){     
         temp2 = Mindex(idx_up,T_down.idx[idx_down][i]);            
         J.accumulate(temp1,temp2, J_down.data[idx_down][i]);
         J.accumulate(temp2,temp1, conj(J_down.data[idx_down][i]));
      }

      for(int i=0;i<T_up.idx[idx_up].size();i++){
         temp2 = Mindex(T_up.idx[idx_up][i], idx_down);
         J.accumulate(temp1, temp2, J_up.data[idx_up][i]);
         J.accumulate(temp2, temp1, conj(J_up.data[idx_up][i]));
      }
   }
*current = real(vec_left_*(J*vec_in_));
Resize_refill(J, dim_up*dim_down, C(0));

Resize_refill(K, dim_up*dim_down, C(0));

for(int idx_down=0;idx_down<Spinconfig[SPIN_DOWN].size();idx_down++)
   for(int idx_up=0;idx_up<Spinconfig[SPIN_UP].size();idx_up++){
      temp1 = Mindex(idx_up,idx_down);
      for(int i=0;i<T_down.idx[idx_down].size();i++){     
         temp2 = Mindex(idx_up,T_down.idx[idx_down][i]);            
         K.accumulate(temp1,temp2, T_down.data[idx_down][i]);
         K.accumulate(temp2,temp1, conj(T_down.data[idx_down][i]));
      }

      for(int i=0;i<T_up.idx[idx_up].size();i++){
         temp2 = Mindex(T_up.idx[idx_up][i], idx_down);
         K.accumulate(temp1,temp2, T_up.data[idx_up][i]);
         K.accumulate(temp2,temp1, conj(T_up.data[idx_up][i]));
      }
   }

*kinetic_energy = real(vec_in_*(K*vec_in_));
Resize_refill(K, dim_up*dim_down, C(0)); 
return;
}


template <typename T, typename C>
C Cluster<T,C>::Gimp(int m,int m_, C omega)
{
 return Frac_Green(a_lesser,b_lesser,E0,-1.0,omega-0.025*IM) + Frac_Green(a_greater,b_greater,E0,1.0,omega-0.025*IM);
}

template<typename T,typename C>
Cluster<T,C>::~Cluster(){};


template<typename T, typename C>
vector<C> Cluster<T,C>::vec_c_dag(int site, int spin, vector<C>& vec_state)
{
  vector< map<int,int> > old_(Spinconfig);
  int old_dimup = dim_up, old_dimdown = dim_down, tempidx=0;

  if(spin == SPIN_UP){

     ne_Up++;
     if( ne_Up> nSites) {cout<<" N_up could not be larger than nSites!!"<<endl;return vec_state;}
     Setup_states();
     vector<C> vectemp(nStates,C(0));
     int temp=0,sgn=0;

     for(map<int,int>::iterator idx_down=old_[SPIN_DOWN].begin();idx_down!=old_[SPIN_DOWN].end();++idx_down)
     for(map<int,int>::iterator idx_up=old_[SPIN_UP].begin();idx_up!=old_[SPIN_UP].end();++idx_up){
        temp = c_dag(site,0,idx_up->first);
        tempidx = IndexExtract(abs(temp),spin); 

        if( tempidx != -1 ) 
          vectemp[Mindex(tempidx,idx_down->second)] += (C)sign(temp)*vec_state[idx_up->second+old_dimup*idx_down->second];         
     }
     return vectemp;   
   } else {
  
     ne_Down++;
     if( ne_Down > nSites) {cout<<" N_down could not be larger than nSites!!"<<endl;return vec_state;}
     Setup_states();
     vector<C> vectemp(nStates,C(0));
     int temp=0,sgn=0;

     for(map<int,int>::iterator idx_down=old_[SPIN_DOWN].begin();idx_down!=old_[SPIN_DOWN].end();++idx_down)
     for(map<int,int>::iterator idx_up=old_[SPIN_UP].begin();idx_up!=old_[SPIN_UP].end();++idx_up){
        temp = c_dag(site,0,idx_down->first);
        tempidx = IndexExtract(abs(temp),spin);      
         
        if( tempidx != -1 )  
          vectemp[Mindex(idx_up->second,tempidx)] += (C)sign(temp)*vec_state[idx_up->second+old_dimup*idx_down->second];   
     }
     return vectemp;
   }

}

template<typename T, typename C>
vector<C> Cluster<T,C>::vec_c_(int site, int spin, vector<C>& vec_state)
{

  vector< map<int,int> > old_(Spinconfig); 
  int old_dimup = dim_up, old_dimdown = dim_down, tempidx=0;

  if(spin == SPIN_UP){

     ne_Up--;
     if( ne_Up<0) {cout<<" N_up could not be less than 0!!"<<endl;return vec_state;}
     Setup_states();
     vector<C> vectemp(nStates,C(0));
     int temp=0,sgn=0;

     for(map<int,int>::iterator idx_down=old_[SPIN_DOWN].begin();idx_down!=old_[SPIN_DOWN].end();++idx_down)
       for(map<int,int>::iterator idx_up=old_[SPIN_UP].begin();idx_up!=old_[SPIN_UP].end();++idx_up){
         temp = c_(site,0,idx_up->first);
         tempidx = IndexExtract(abs(temp),spin);

         if( tempidx != -1 )
           vectemp[Mindex(tempidx,idx_down->second)] += (C)sign(temp)*vec_state[idx_up->second+old_dimup*idx_down->second];
       }
     return vectemp; 
   } else {

     ne_Down--;
     if( ne_Down < 0) {cout<<" N_down could not be less than 0!!"<<endl;return vec_state;}
     Setup_states();
     vector<C> vectemp(nStates,C(0));
     int temp=0,sgn=0;

     for(map<int,int>::iterator idx_down=old_[SPIN_DOWN].begin();idx_down!=old_[SPIN_DOWN].end();++idx_down)
       for(map<int,int>::iterator idx_up=old_[SPIN_UP].begin();idx_up!=old_[SPIN_UP].end();++idx_up){
         temp = c_dag(site,0,idx_down->first);
         tempidx = IndexExtract(abs(temp),spin);
         
         if( tempidx != -1 ) 
           vectemp[Mindex(idx_up->second,tempidx)] += (C)sign(temp)*vec_state[idx_up->second+old_dimup*idx_down->second];
       }
   return vectemp;
   }

}


 
template<typename T,typename C> 
void Cluster<T,C>::Lapack_diag()
{

vector<C> Evecs(nStates*nStates,C(0)), temp_left_(nStates,C(0));
vector<T> Er(nStates,T(0));

Setup_states();
MElements_nonsparse();
zheevr(Hamiltonian , Er, Evecs, nStates);

int i=0;

temp_left_.assign(&Evecs[i*nStates],&Evecs[i*nStates+nStates]);
cout<<"********** Full diag Eigenvalues ***************"<<endl;
cout<<temp_left_<<endl;
cout<<"********** End of Full Diag ****************"<<endl;

i=2;
temp_left_.assign(&Evecs[i*nStates],&Evecs[i*nStates+nStates]);

cout<<temp_left_<<endl;

return ;
}

