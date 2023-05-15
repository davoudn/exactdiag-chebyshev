#include <iostream>
#include "spmatrix.h"
#include "ioclass.h"
#include <math.h>



//#include "complex.h"
using namespace std;
template <typename T>
int Lanczos(matrix<T>& H, int dimension, int Maxiter, vector<T>& a, vector<T>& b, vector<T>& v0);

template <typename T>
int Lanczos(matrix<T>& H, int dimension, int Maxiter, vector<T>& a, vector<T>& b, vector<T>& v0, T* E_min, T* E_max);

template <typename T>
int Lanczos(CEll<T>& H, int dimension, int Maxiter, vector<T>& a, vector<T>& b, vector<T>& v0, T* E_min, T* E_max);


template <typename T>
void ab2matrix(matrix<T>& c, vector<T> a, vector<T> b); 
// lanczos                   

template <typename T>
void ab2matrix(matrix<T>& c, vector<T> a, vector<T> b)
{
  c[a.size()][a.size()] = a[a.size()];

c[0][0] = a[0];
c[a.size()-1][a.size()-1] = a[a.size()-1];
 for(int i=1;i<a.size()-1;i++)
 {
    c[i][i] = a[i];
    c[i][i+1] = b[i]; 
    c[i+1][i] = b[i];
 }
return;
}

template <typename T>
matrix<T> ab2matrix(vector<T>& a,vector<T>& b)
{
matrix<T> temp(a.size(),a.size(),T(0));

 for(int i=0;i<a.size();i++)
 {
    temp[i][i] = a[i];
    if(i!= a.size()-1)
    {
      temp[i][i+1] =  b[i+1];
      temp[i+1][i] =  b[i+1];
    }
 }
return temp;
}
                     
template <typename T>
T norm(vector<T>& v)
{
T temp=T(0);
for(int i=0;i<v.size();i++) temp+=pow(abs(v[i]),2.0);

return  pow(temp,0.5);
}

template <typename T>
void normalize(vector<T>& v)
{
T temp=T(0);
for(int i=0;i<v.size();i++) temp+=pow(abs(v[i]),2.0);
temp = pow(pow(temp,0.5),-1.0);
v=temp*v;

return;
}

template <typename T>
void fill(vector<T>& v, T t)
{

for(int i=0;i<v.size();i++) v[i]=t;

return;
}

template <typename T>
int Lanczos(matrix<T>& H, int dimension, int Maxiter, vector<T>& a, vector<T>& b, vector<T>& v0, T* E_min, T* E_max)
{
 vector<T> v1(dimension, T(0));
 vector<T> w(v0);
 vector<T> z(Maxiter*Maxiter,T(0));
 vector<T> e(Maxiter,T(0));

 fstream is;
 is.open("convergence", fstream::out);
 
 vector<T> a1,b1;

 int i=0,i_=0;
 static double eps = 0.000000001;
 T e_new = T(0); 
 T e_old = T(0);

 a.clear();b.clear();
 b.push_back(T(0));

 normalize(v0);

/* first lanczos proceduer for the calculating treediagonal matrix coeficents a[i],b[i] */
 v1 = H*v0;
 a.push_back(v1*v0);
 v1 = v1 - a[0]*v0;
 b.push_back(norm(v1));
 
 cout<<"Iteration   Ground-State energy     Relative error"<<endl;
 cout.precision(13);
 cout<<fixed<<endl; 
 for(i=1;i<Maxiter;i++)
 {  
    v0 = H*v1 - b[i]*v0;
    a.push_back(v1*v0); 
    v0 = v0 - a[i]*v1;
    b.push_back(norm(v0));     
    normalize(v0);
    v1.swap(v0);
/* a block for convergence check   */
       
    e.resize(a.size(), T(0));
    z.resize(a.size()*a.size(),T(0));
 
    a1=a;b1=b;
    dstemr(a1, b1, a1.size(), e, z, a1.size());
      
    e_new = e[0]; 
    is<<i<<" "<<abs(e_old-e_new);      
    for(int s=0;s<e.size();s++) 
       is<<" "<<e[s];
    is<<endl;
        
    cout<<"Ground-State energy     Relative error"<<endl;
    cout<<i<<"       "<<e[0]<<"                  "<<abs(e_old-e_new)<<"b[i]="<<b[i]<<" a[i]="<<a[i]<<endl;
    if(abs(e_old-e_new)<eps) 
    {
       cout<<"converged"<<endl;
       break;
    }
    e_old = e_new;
 }     
 i_ = i;
 is.close();
// calculation of ground state of the tridiagonal matrix 
 a1 = a;
 b1 = b;
 dstemr(a1, b1, a.size(), e, z, a.size());//a.size());
/* calculation of ground state of the system it is down by restarting the lanczos algourithm the difference is that
   now the the ground state of the tridiagonal system exists thus the contribution of each lanczos vector is known
   and the only thing that has to be down is to multiplay these contribution to the lanczos vectors and sum up at each iteration
   it needs only one vector more that should be allocated 
*/
 v0 = w;
 w = vector<T>(dimension,T(0));
 normalize(v0);
 w = w + z[0]*v0; // w is the ground state of the system that should be computed
 v1 = H*v0;
 v1 = v1 - a[0]*v0;
 normalize(v1);

 w = w + z[1]*v1;

 for(i=1;i<i_;i++)
 {
    v0 = H*v1 - b[i]*v0;
    v0 = v0 - a[i]*v1;
    normalize(v0);
    v1.swap(v0);
    if(i<Maxiter-1) 
      w = w + z[i+1]*v1;     
 }
 normalize(w);
 v0 = w;
 *E_min = e[0];
 *E_max = e[0];

  cout.precision(20);
  cout<<"ITERATIONS=="<<i_<<endl;
  cout<<"Lanczos ground state=="<<e[0]<<endl;
  cout<<"abs(w*(H*w)-e[0])="<<abs(w*(H*w)-e[0])<<endl;

return 0;
}
// *******************************************************************************


template <typename T>
int Lanczos(CEll<T>& H, int dimension, int Maxiter, vector<T>& a, vector<T>& b, vector<T>& v0, T* E_min, T *E_max)
{
 vector<T> v1(dimension, T(0));
 vector<T> w(v0);
 vector<T> z(Maxiter*Maxiter,T(0));
 vector<T> e(Maxiter,T(0));

 fstream is;
 is.open("convergence", fstream::out);
 
 vector<T> a1,b1;

 int i=0,i_=0;
 static double eps = 1.0e-8;
 T e_new = T(0); 
 T e_old = T(0);

 a.clear();b.clear();
 b.push_back(T(0));

 normalize(v0);
// first lanczos proceduer for the calculating treediagonal matrix coeficents a[i],b[i] 
 v1 = H*v0;
 a.push_back(v1*v0);
 v1 = v1 - a[0]*v0;
 b.push_back(norm(v1));
 normalize(v1);
 
// cout<<"Iteration   Ground-State energy     Relative error"<<endl;
 cout.precision(13);
 cout<<fixed; 
 for(i=1;i<Maxiter;i++)
 {  
    v0 = H*v1 - b[i]*v0;
    a.push_back(v1*v0); 
    v0 = v0 - a[i]*v1;
    b.push_back(norm(v0));     
    normalize(v0);
    v1.swap(v0);
// a block for convergence check   
    e.resize(a.size(), T(0));
    z.resize(a.size()*a.size(),T(0));
 
    a1=a;b1=b;
    dstemr(a1, b1, a1.size(), e, z, a1.size());
      
    e_new = e[0]; 
    is<<i<<" "<<abs(e_old-e_new);      
    for(int s=0;s<e.size();s++) 
       is<<" "<<e[s];
    is<<endl;

    cout<<"Ground-State energy     Relative error"<<endl;
    cout<<i<<"       "<<e[0]<<"                  "<<abs(e_old-e_new)<<"b[i]="<<b[i]<<" a[i]="<<a[i]<<endl;
      
    if(abs(e_old-e_new)<eps) 
      break;
    e_old = e_new;
 }     
 i_ = i;
 is.close();
 a1 = a;
 b1 = b;
 dstemr(a1, b1, a.size(), e, z, a.size());//a.size());
// calculation of ground state of the system it is down by restarting the lanczos algourithm the difference is that
//   now the the ground state of the tridiagonal system exists thus the contribution of each lanczos vector is known
//   and the only thing that has to be down is to multiplay these contribution to the lanczos vectors and sum up at each //iteration
//   it needs only one vector more that should be allocated 

 v0 = w;
 w = vector<T>(dimension,T(0));
 normalize(v0);
 w = w + z[0]*v0; // w is the ground state of the system that should be computed
 v1 = H*v0;
 v1 = v1 - a[0]*v0;
 normalize(v1);
 w = w + z[1]*v1;

 for(i=1;i<i_;i++)
 {
    v0 = H*v1 - b[i]*v0;
    v0 = v0 - a[i]*v1;
    normalize(v0);
    v1.swap(v0);
    if(i<Maxiter-1) 
      w = w + z[i+1]*v1;     
 }
 normalize(w);
 v0 = w;
 *E_min = e[0];
 *E_max = e[e.size()];

return 0;
}

 
template <typename T, typename C>
int Lanczos1(CEll<C>& H, int dimension, int Maxiter, vector<C>& a, vector<C>& b, vector<C>& v0, T* E_min,T* E_max, int calculate_groundstate)
{
 vector<C> v1(dimension, C(0));
 vector<C> w(v0);
 vector<C> z(1,C(0));
 vector<T> e(Maxiter,T(0));
 matrix<C> abmatrix(1,1,C(0));
 fstream is;
 vector<C> a1,b1;
 T e_new = T(0); 
 T e_old = T(0);

 is.open("convergence", fstream::out);
 int i=0,i_=0;
 static double eps = 1.0e-13;
 a.clear();b.clear();
 b.push_back(C(0));

 normalize(v0);
 
//first lanczos proceduer for the calculating treediagonal matrix coeficents a[i],b[i] 
 v1 = H*v0;
 a.push_back(v1*v0);
 v1 = v1 - a[0]*v0;
 b.push_back(norm(v1));
 normalize(v1);
 
 for(i=1;i<Maxiter;i++){  
       v0 = H*v1 - b[i]*v0;
       a.push_back(v1*v0); 
       v0 = v0 - a[i]*v1;
       b.push_back(norm(v0));     
       normalize(v0);
       v1.swap(v0);
//a block for convergence check   
       z.resize(a.size()*a.size(), C(0));
       e.resize(a.size(),T(0));

       abmatrix = ab2matrix(a,b);
       zheevr(abmatrix, e, z,a.size());
       a1=a;b1=b;
       e_new = e[0]; 
//     is<<i<<" "<<abs(e_old-e_new);      
//     for(int s=0;s<e.size();s++) 
//        is<<" "<<e[s];
//     is<<endl;
       if(abs(e_old-e_new)<eps)  
         break;
       e_old = e_new;
 }     
 i_ = i;
 is.close();

// calculation of ground state of the tridiagonal matrix 
 a1 = a;
 b1 = b;
 abmatrix = ab2matrix(a,b);
 zheevr(abmatrix, e, z,a.size());

 if(calculate_groundstate ==1)
 {
    v0 = w;
    w = vector<C>(dimension,C(0));
    normalize(v0);
    w = w + z[0]*v0; // w is the ground state of the system that should be computed
    v1 = H*v0;
    v1 = v1 - a[0]*v0;
    normalize(v1);
    w = w + z[1]*v1;

    for(i=1;i<i_;i++)
    {
       v0 = H*v1 - b[i]*v0;
       v0 = v0 - a[i]*v1;
       normalize(v0);
       v1.swap(v0);
       if(i<Maxiter-1) 
         w = w + z[i+1]*v1;     
    }
    normalize(w);
    v0 = w;
    cout<<endl;
    cout<<" >>> Some usefull information about Lanczos run: <<<"<<endl;
    cout<<" ==> abs(<psi_0|H|psi_0>-e[0])="<<abs(w*(H*w)-e[0])<<endl;
 } 
 *E_min = e[0];
 *E_max = e[e.size()-1];
 cout.precision(15);
 cout<<" ==> Emin(e[0])="<<*E_min<<", Emax="<<*E_max<<endl;
 cout<<" ==> First exited="<<e[1]<<", "<<"Gap="<<e[1]-e[0]<<endl;
 cout<<" ==> Iterations="<<i_<<endl;
    cout<<" >>> ------------------------------------------- <<<"<<endl;

return 0;
}


template <typename T, typename C>
int Lanczos_fast(CEll<C>& H, int dimension, int Maxiter, vector<C>& a, vector<C>& b, vector<C>& v0, T* E_min,T* E_max, T eps, int calculate_groundstate)
{
 vector<C> v1(dimension, C(0));
 vector<C> w(v0);
 vector<C> z(1,C(0));
 vector<T> e(Maxiter,T(0));
 matrix<C> abmatrix(1,1,C(0));
 vector< vector<C> > basis(0,vector<C>(v0.size(),C(0)));
 fstream is;
 vector<C> a1,b1;
 T e_new = T(0);
 T e_old = T(0);

 is.open("convergence", fstream::out);
 int i=0,i_=0;
 a.clear();b.clear();
 b.push_back(C(0));

 v0 = vector<C>(v0.size(),C(1));
 normalize(v0);
 basis.push_back(v0);
//first lanczos proceduer for the calculating treediagonal matrix coeficents a[i],b[i] 
 v1 = H*v0;
 a.push_back(v1*v0);
 v1 = v1 - a[0]*v0;
 b.push_back(norm(v1));
 normalize(v1);
 basis.push_back(v1);

 for(i=1;i<Maxiter;i++){
       v0 = H*v1 - b[i]*v0;
       a.push_back(v1*v0);
       v0 = v0 - a[i]*v1;
       b.push_back(norm(v0));
       normalize(v0);
       v1.swap(v0);
       basis.push_back(v1);
//a block for convergence check   
       z.resize(a.size()*a.size(), C(0));
       e.resize(a.size(),T(0));

       abmatrix = ab2matrix(a,b);
       zheevr(abmatrix, e, z,a.size());
       a1=a;b1=b;
       e_new = e[0];
       is<<i<<" "<<abs(e_old-e_new);
       for(int s=0;s<e.size();s++)
          is<<" "<<e[s];
       is<<endl;
       if(abs(e_old-e_new)<eps)
         break;
       e_old = e_new;
 }
 i_ = i;
 is.close();

// calculation of ground state of the tridiagonal matrix 
 a1 = a;
 b1 = b;
 abmatrix = ab2matrix(a,b);
 zheevr(abmatrix, e, z,a.size());

 if(calculate_groundstate ==1){
    w = vector<C>(dimension,C(0));

    for(i=0;i<a.size();i++)
       if(i<Maxiter-1)
         w = w + z[i]*basis[i];
    
    normalize(w);
    v0 = w;
    cout<<"abs(w*(H*w)-e[0])="<<abs(w*(H*w)-e[0])<<endl;
 }
 *E_min = e[0];
 *E_max = e[e.size()-1];

 cout<<"ITERATIONS="<<i_<<endl;
 cout<<"Emin"<<*E_min<<" Emax="<<*E_max<<endl;

 return 0;
}



template<typename T,typename C>
void Lapack_step(matrix<C>&  H, vector<C>& Groundvecstate_0_, T time_step,int state_num)
{
  cout<<endl;
  cout<<"***** Starting Lapack evolution *****"<<endl;
  int N = H.cols();

//zheevr(H, temp_Es ,temp_Evecs, state_num);
/*   ****************************************************************** */

int n = H.cols();
char jobz,range,uplo;
int il,iu,m,ldz,lda,lwork,lrwork,info,liwork;
T vl,abstol,vu;
T *w,*rwork;
C *a,*z,*work;
int *isuppz,*iwork;

jobz = 'V';
range= 'I';
uplo=  'U';
a=ctof(H);
lda= n;
il=1;iu=state_num;
abstol=0.0;
info=0;
ldz=n;
z=new C[n*n];
lwork=2*n;
//lrwork = 2*n;
lrwork = 24*n;
liwork=10*n;
work = new C[lwork];
iwork = new int[liwork];
w =new T[n];
rwork = new T[lrwork];
isuppz = new int[2*n];
T *temp_Es = new T[n];
C *temp_Evecs = new C[n*n];

zheevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m,
   &temp_Es[0],&temp_Evecs[0],&ldz,isuppz,work,&lwork, rwork, &lrwork, iwork,&liwork,&info);

delete z;
delete work;
delete iwork;
delete w;
delete isuppz;
delete a;


/*   ****************************************************************** */
// time evolution part
  C temp1=C(0);
  vector<C> temp(N,C(0));
  for(int i=0;i<state_num;i++)
  {
     temp1 = VV_Multi(&temp_Evecs[i*N],&Groundvecstate_0_[0],N);
     temp  = temp + SV_Multi(exp(-IM*time_step*temp_Es[i])*temp1,&temp_Evecs[i*N],N);
  }
Groundvecstate_0_ = temp;
cout<<"Ground state energy of new hamiltoni"<<temp_Es[0]<<endl;
cout<<"Norm =="<<norm(Groundvecstate_0_)<<endl;
cout<<"***** End of Lapack evolution *****"<<endl;
cout<<endl;

delete temp_Es;
delete temp_Evecs;


return ;
}

/*   kernel polynomial 
************************************************
*/
template<typename T,typename C>
void Chebychev_step(CEll<C>& H,  vector<C>& Groundvecstate_0_, T time_step,T tol, T E_min_in, T E_max_in)
{
 cout<<" >>>---- Strating Chebychev evolution ----<<<"<<endl;
 T temp_bessel=T(0);
 C temp1=C(0);
 int Maxiter =100, dim=H.idx.size(),temp_bessel_count=0;
 vector<C> vinit(dim,C(1));
 vector<C> a_0,b_0;
 vector<T> alpha;
 
 T A = (E_max_in-E_min_in)/2.0;
 T B = (E_max_in+E_min_in)/2.0;

 cout<<" ==> A= "<<A<<" B="<<B<<endl;
    
 while( 0<1 ){
      temp_bessel = jn(temp_bessel_count, A*time_step); 
      if(abs(temp_bessel)>tol){
           alpha.push_back(temp_bessel);
           temp_bessel_count++;
      }
      else 
        break;
 }
     cout<<" ==> Number of bessel evaluations="<<temp_bessel_count<<endl;
 if(temp_bessel_count != 0)
   Groundvecstate_0_ = chebishef_seri_eval(H, Groundvecstate_0_, alpha);

 temp1 = exp(-IM*B*time_step);
 Groundvecstate_0_=temp1*Groundvecstate_0_;

 cout<<" ==> Norm="<<norm(Groundvecstate_0_)<<endl;
 cout<<" >>>-----  End of Chebychev evolution  -----<<<"<<endl;
return;
}


template<typename T,typename C>
vector<C> Chebychev_propagation(CEll<C>& H,  vector<C>& Groundvecstate_0_, T time_in, T time_step, T tol, T E_min_in, T E_max_in, vector< vector<C> >& work_space)
{
 cout<<"***** Strating Chebychev evolution *****"<<endl;
 T temp_bessel=T(0);
 vector<C> a_0,b_0;
 int Maxiter =100, dim=H.idx.size(),temp_bessel_count=0,diff=0;
 vector<T> alpha;
 vector<C> temp_v(Groundvecstate_0_);
 C temp1=C(0);
 
 T A = (E_max_in-E_min_in)/2.0;
 T B = (E_max_in+E_min_in)/2.0;
    
 cout<<"A= "<<A<<" B="<<B<<endl;

 while( 0<1 ){
      temp_bessel = jn(temp_bessel_count, A*time_step);
      if(temp_bessel>tol){
        alpha.push_back(temp_bessel);
        temp_bessel_count++;
      }
      else
        break;
 }
 cout<<" Number of bessel evaluations="<<temp_bessel_count<<endl;
 
 for(T temp_time=time_step;temp_time<=time_in;temp_time+=time_step){
    if(time_step != time_in)
      work_space.clear();
    if(temp_bessel_count != 0){
      diff = temp_bessel_count-work_space.size();
      chebishef_vecs_eval(H,temp_v, work_space,diff);
    }
    temp_v.clear();temp_v.resize(H.size());
    for(int i=0; i<work_space.size(); i++)
       temp_v  = temp_v +  (C(2)*pow(-IM,(C)i)*alpha[i])*work_space[i];
        
    temp1 = exp(-IM*B*time_step);
    temp_v=temp1*temp_v;
 }
 cout<<"Norm =="<<norm(Groundvecstate_0_)<<endl;
 cout<<"***** End of Chebychev evolution *****"<<endl;
return temp_v;
}

template<typename T>
vector<T> Tn(CEll<T>& H , vector<T>& v0 , vector<T>& v1) 
 {   
   vector<T> temp( v0.size() , T(0));
   temp = (H*v1);
   temp = T(2)*temp - v0;
   return temp;
 }


template<typename T>
void chebishef_vec_filling(matrix<T>& H ,vector<T>& vi ,vector<T>& vj ,T* chebishef_coef, int Evol_maxiter)
{
     T m ,m0 , m1, mtemp = T(0);
     int sz = vi.size();
     int n = Evol_maxiter;
     vector<T> vi0(sz , T(0));
     vector<T> vi1(sz , T(0));
     vector<T> vtemp(sz , T(0));
 
     m0 = vj*vi;
     chebishef_coef[0]=m0;
     vi0 = vi;
     vi1 = H*vi0; m1 = vj*vi1;
     chebishef_coef[1] = m1;

     for(int i=2  ;i<n ;i++){
        vtemp = Tn(H ,vi0 , vi1 );
        vi0=vi1;
        vi1=vtemp;
        mtemp =vj*vtemp;
        chebishef_coef[i] = mtemp;
     }
return;
}

template<typename T>
void chebishef_coef_filling(matrix<T>& H ,vector<T>& vi ,vector<T>& vj ,T* chebishef_coef, int Evol_maxiter)
{
     T m ,m0 , m1, mtemp = T(0);
     int sz = vi.size();
     int n = Evol_maxiter;
     vector<T> vi0(sz , T(0));
     vector<T> vi1(sz , T(0));
     vector<T> vtemp(sz , T(0));

     m0 = vj*vi;
     chebishef_coef[0]=m0;
     vi0 = vi;
     vi1 = H*vi0; m1 = vj*vi1;
     chebishef_coef[1] = m1;

     for(int i=2  ;i<n ;i++){
        vtemp = Tn(H ,vi0 , vi1 );
        vi0=vi1;
        vi1=vtemp;
        mtemp =vj*vtemp;
        chebishef_coef[i] = mtemp;
     }   
return;
}

template<typename T,typename C>
vector<C> chebishef_seri_eval(CEll<C>& H, vector<C>& temp_v0, vector<T> alpha)
{
        int n = H.idx.size();

        vector<C> temp_v1(n, C(0));
        vector<C> temp_v(n, C(0));

        temp_v = ((C)alpha[0])*temp_v0;
        temp_v1 = H*temp_v0;
        temp_v = temp_v + C(2)*((-IM)*alpha[1])*temp_v1;
        
        for(int i=2; i<alpha.size(); i++){
           temp_v0 = Tn(H, temp_v0, temp_v1);
           temp_v  = temp_v +  (C(2)*pow(-IM,(C)i)*alpha[i])*temp_v0;
           temp_v1.swap(temp_v0);
        }

return temp_v;
}

template<typename C>
void chebishef_vecs_eval(CEll<C>& H,vector<C>& temp_v0, vector< vector<C> >& work_space,int diff)
{
 int n = H.idx.size();
 vector<C> temp_v1(n, C(0));  
 bool flag = true;

 if(work_space.size() == 0){
   temp_v1 = H*temp_v0;
   work_space.push_back(temp_v0);
   work_space.push_back(temp_v1);
   flag = false;
 }

 if(flag)
   for(int i=0; i<diff; i++){
      temp_v0 = Tn(H, work_space[work_space.size()-2], work_space[work_space.size()-1]);
      work_space.push_back(temp_v0); 
   }

return;
}


template<typename T>
T Aw(T w ,vector<T>& chebishef_coef)
{
 int sz = chebishef_coef.size();
 T temp,s = T(0);
 temp = chebishef_coef[0] ;
 if(abs(w)>1.0) 
   return 0.0;
 else
 {
   for(int i=1;i<sz;i++) 
   {
      s++;
      temp+=2.0*chebishef_coef[i]*ch_kernel(3.0, s ,(T) sz)*cos(s*acos(w));
   }
   temp *= pow(sqrt(1.0-w*w),-1.0) ;
 }
 return temp;
}

template<typename T, typename C>
void Hscale(CEll<C>& H,T e_min , T e_max )
{
  T a = (e_max - e_min)/2.0;
  T b = (e_min + e_max)/2.0;

H -= b; 
H *= (1.0/a);
  return ;
}


template<typename T>
T kpmfdelta(T x, T N)
{
  T v = 1.0;
  double landa = 2.0;
  for(double i=1;i<N; i++)
     v += 2*cos(i*acos(x))*ch_kernel(landa,i,N)*cos(i*acos(0));       
  return pow(3.14*sqrt(1.0-x*x),-1.0)*v;
}

template<typename T>
T ch_kernel(T landa ,T i ,T N)
{
 T temp= T(0);
 temp =  sinh(landa*(1.0-(i/N)))/sinh(landa);
 return temp;
}

