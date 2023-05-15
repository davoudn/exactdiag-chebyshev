#include <vector>
using namespace std;

template <typename T>
class Ell
{
public:
Ell(int rw,int nz);
~Ell();
void accumulate(int row, int col, T dta);

int cols;
int rows;
vector<int> idx;
vector<T> data;

};

template <typename T>
Ell<T>::Ell(int rw, int nz  ):idx(rw*nz,-1),data(rw*nz,T(0))
{
 rows=rw;
 cols=nz;
};


template <typename T>
Ell<T>::~Ell(){};

template <typename T>
void Ell<T>::accumulate(int row, int col, T dta)
{
int count=0;
 for(int i=0;i<cols;i++) 
   {
     if(idx[row*cols+i] != -1) count++;
     if(idx[row*cols+i] == col ) {data[row*cols+i]+= dta;count=-1;break;}
   }

 if(count != -1){
  data[row*cols+count]+= dta;
  idx[row*cols+count] = col;
 }
return;
}
// *********************************************************************
template <typename T>
class CEll
{
public:
CEll(int dim);
CEll();
~CEll();


int cols;
int rows;
vector< vector<int> > idx;
vector< vector<T> > data;
void accumulate(int row, int col, T dta);
void resize(int dim);
int size();
void CEllRow_CEllCol_major();
int tot_nze();
CEll& operator+= (const T& a);
CEll& operator-= (const T& a);
CEll& operator*= (const T& a);
CEll& operator= (const CEll<T>& A);

};

template <typename T>
CEll<T>::CEll(int dim):idx(dim,vector<int>(0,0)),data(dim,vector<T>(0,T(0))){rows=dim;};

template <typename T>
CEll<T>::CEll(){};

template <typename T>
CEll<T>::~CEll(){};

template <typename T>
void CEll<T>::resize(int dim)
{
idx.resize(dim, vector<int>(0,0));
data.resize(dim, vector<T>(0,T(0)));
return;
}

template <typename T>
int CEll<T>::tot_nze()
{
int temp=0;

for(int i=0;i<this->idx.size();i++)
   for(int j=0;j<this->idx[i].size();j++)
      temp++;
return temp;
}

template <typename T>
void CEll<T>::accumulate(int row, int col, T dta)
{
int count = 0;
for(int i=0;i<idx[row].size();i++) if(idx[row][i] == col ) {data[row][i]+= dta;count=-1;break;}

 if(count != -1){
  data[row].push_back(dta);
  idx[row].push_back(col);
 }
return;
}

template<typename T>
int CEll<T>::size()
{

return this->idx.size();
}

template <typename T>
void CEll<T>::CEllRow_CEllCol_major()
{

CEll<T> *temp= new CEll(this->size()); 

for(int i=0;i<this->size();i++)
   for(int j=0;j<this->idx[i].size();j++)
      temp->accumulate(this->idx[i][j],i,this->data[i][j]);

this->resize(temp->size());

for(int i=0;i<temp->size();i++)
   for(int j=0;j<temp->idx[i].size();j++)
      this->accumulate(temp->idx[i][j],i,temp->data[i][j]);

delete temp;

return;
}

template<typename T>
CEll<T>& CEll<T>::operator*= (const T& a)
{
  for(int i=0;i<data.size();i++)
     for(int j=0;j<data[i].size();j++)  data[i][j]*= a;
  return *this;
}

template<typename T>
CEll<T>& CEll<T>::operator= (const CEll<T>& A)
{
  this->resize(A.data.size());  
  for(int i=0;i<A.data.size();i++)
     for(int j=0;j<A.data[i].size();j++){this->accumulate(i,A.idx[i][j],A.data[i][j]);}
  return *this;
}

template<typename T>
CEll<T>& CEll<T>::operator-= (const T& a)
{
  for(int i=0;i<data.size();i++)
     for(int j=0;j<data[i].size();j++)
        if(idx[i][j] == i)  data[i][j]-= a;

  return *this;
}

template<typename T>
CEll<T> operator* (const CEll<T>& A , const CEll<T>& B)
{ 
  T tsum = T(0);  
  CEll<T> ctemp(A.data.size());
 
  for(int i=0;i<A.idx.size();i++){
     for(int j=0;j<A.idx[i].size();j++)
        for(int k=0;k<B.idx[ A.idx[i][j] ].size();k++)
           ctemp.accumulate(i,B.idx[ A.idx[i][j] ][k], A.data[i][j]*B.data[ A.idx[i][j] ][k]);
  }
return ctemp;
}
// ******************************************************************* //

// ******************************************************************* //
template<typename T>
class CSR
{

public:
CSR(int nnz_,int N_);
CSR(){};
~CSR();

void clear();
int N,nnz;
vector<T> sa; //factor 2 for complex 
vector<int> isa;
vector<int> jsa;
  /*!!!!!!!!!!!!!!!!! Others */
  // !!!!!!!!!! form CSR arrays isa,jsa,sa 
};
template<typename T>
CSR<T>::CSR(int nnz_,int N_)
{
  N = N_;
  nnz=nnz_;
};

template<typename T>
CSR<T>::~CSR(){};

template<typename T>
void CSR<T>::clear()
{
this->sa.clear();
this->jsa.clear();
this->isa.clear();
return;
}

// ******************************************************************* //

// ******************************************************************* //
template<typename T>
vector<T> operator* (const CEll<T>& A , const vector<T>& V)
{
  T tsum = T(0);    
  vector<T> vtemp;
  for(int i=0;i<A.idx.size();i++){
     for(int j=0;j<A.data[i].size();j++)  
        tsum+= A.data[i][j]*V[ A.idx[i][j] ];
     vtemp.push_back(tsum);
     tsum = T(0);
  }
return vtemp;
}

template<typename T>
void MV_Multi(CEll<T>& A , T *V, T *W)
{
 T tsum = T(0);
 for(int i=0;i<A.idx.size();i++){
    for(int j=0;j<A.data[i].size();j++)  
       tsum+= A.data[i][j]*V[ A.idx[i][j] ];
    W[i] = tsum;
    tsum = T(0);
 }
        return ;
}

template<typename T>
T VV_Multi(T *V, T *W, int N)
{
 T temp=T(0);
     for(int i=0;i<N;i++)
        temp+= conj(V[i])*W[i];
return temp;
}

template<typename T>
vector<T> SV_Multi(T a, T *V, int N)
{
 vector<T> temp(N,T(0));
     for(int i=0;i<N;i++)
        temp[i]+= a*V[i];
return temp;
}


/*  new sparse matrix */

template <typename T>
class Sparse_hamiltonian
{
public:
Sparse_hamiltonian(int dim);
Sparse_hamiltonian();
~Sparse_hamiltonian();

int cols;
int rows;
CEll<T> T_up,T_down;
vector<int> Spinupconfig,Spindownconfig;
T U;
vector<T> data(int idx);
vector<T> index(int idx);

void resize(int dim);
};


template <typename T>
Sparse_hamiltonian<T>::Sparse_hamiltonian(){};

template <typename T>
Sparse_hamiltonian<T>::~Sparse_hamiltonian(){};

template <typename T>
vector<T> Sparse_hamiltonian<T>::data(int idx)
{
return;
}

template <typename T>
vector<T> Sparse_hamiltonian<T>::index(int idx)
{
return;
}

template<typename T>
 vector<T> operator* (const Sparse_hamiltonian<T>& A , const vector<T>& V)
{
        T tsum = T(0);
       
        vector<T> vtemp,temp_data,temp_idx;

        for(int i=0;i<A.size();i++)
        {
         temp_data = A.data(i);
         temp_idx = A.index(i);
         for(int j=0;j<temp_data.size();j++)  tsum+= temp_data[j]*V[ temp_idx[j] ];
                  vtemp.push_back(tsum);
                  tsum = T(0);
                }
        return vtemp;
}

