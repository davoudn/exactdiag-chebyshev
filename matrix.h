#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <complex>
#include <math.h>
#include <string>
#define 	F77_FUNC(name, NAME)   name ## _
# define PRECISION 2 

using namespace std;

template <typename T> class matrix;
template <typename T>
class matrix
{
public:
	matrix(int rows=1 , int cols =1 , T t = T());
	matrix(const matrix<T>& A);
	~matrix();
	//matrix& operator+= (const matrix& A);
	matrix& operator+= (const matrix& A);
	matrix& operator-= (const matrix& A);
	matrix& operator= (const matrix& A);
	matrix& operator*= (const matrix& A);
	matrix  operator-() const;
        matrix  operator+() const;
        vector< vector<T> > v;
    vector<T>& operator[] (int i);
    const vector<T>& operator[] (int i) const;
    int rows() const;
    int cols() const;
    void reserve(int height,int lenth) ;
   	
//private:
 //   vector< vector<T> > v;

};


/* External Lapack call from fortran       */
extern "C" {    
  void zgetri_(int *n, complex<double> *a, int *lda, int *ipiv,
      complex<double> *work, int *lwork, int *info);
}

extern "C" void dgeev_(char *jobvl, char *jobvr, int *n, double *a,
    int *lda, double *wr, double *wi, double *vl,
    int *ldvl, double *vr, int *ldvr,
    double *work, int *lwork, int *info);


extern "C" void dsyevr_(char *jobz, char *range, char *uplo, int *n,
    double *a, int *lda, double *vl,double *vu,int *il,int *iu ,double *abstol, int *m, 
    double *w,double *z, int *ldz,int *isuppz,double *work,int *lwork ,int *iwork,int *liwork,int *info);


extern "C" void sstevx_(char *jobz, char *range, int *n,double *d,double *e, double *vl,double *vu,int *il,int *iu ,double *abstol, int *m,
    double *w,double *z, int *ldz ,double *work ,int *iwork ,int *ifail ,int *info);

extern "C" void sstev_(char *jobz, int *n,double *d,double *e ,double *z, int *ldz ,double *work ,int *info);

extern "C" void dstemr_(char *jobz, char *range, int *n, double *d, double *e, double *vl, double *vu, int *il, int *iu, int *m, double *w, double *z, int *ldz, int *nzc, int *isuppz, bool *tryrac, double *work, int *lwork, int *iwork, int *liwork, int *info);


extern "C" void zheevr_( char* jobz, char* range, char* uplo, int *n,
      complex<double> *a, int* lda, double *vl, double *vu, int *il, int *iu, double *abstol, int *m,
      double *w, complex<double> *z, int *ldz, int *isuppz, complex<double> *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork, int *info);




/* non class Arithmetic operator definitons 
*************************************
*/
template<typename T>
matrix<T> operator+ (const matrix<T>& A , const matrix<T>& B);

template<typename T>
matrix<T> operator- (const matrix<T>& A , const matrix<T>& B);

template<typename T>
matrix<T> operator* (const matrix<T>& A ,const matrix<T>& B);

template<typename T>
vector<T> operator* (const matrix<T>& A ,const vector<T>& V );

template<typename T>
matrix<T> operator* (const T& a ,const matrix<T>& A );

template<typename T,typename C>
vector<C> operator* (const T& a ,const vector<C>& A );

template<typename T>
vector<T> operator- (const vector<T>& A ,const vector<T>& B );

template<typename T>
vector<T> operator+ (const vector<T>& A ,const vector<T>& B );

template<typename T>
T operator* (const vector<T>& A ,const vector<T>& B );


template<typename T>
matrix<T> operator- (const matrix<T>& A , const T& b);


/*Standard input output operator definitions
**************************************
*/
template <typename T>
ostream& operator<< (ostream& , matrix<T>& A);

template <typename T>
ostream& operator<< (ostream& , vector<T>& V);

template <typename T>
istream& operator>> (istream& , matrix<T>& A);


/*    non class function definitions 
***************************************
*/


/*  lapack matrix actions   
****************************************************
*/

template<typename T>
void dgeev(matrix<T>& H, vector<T>& Er, vector<T>& Ei, matrix<T>& Evecs);

template<typename T>
matrix<T> Inv(const matrix<T>& A);

template<typename T>
void dsyevr(matrix<T>& H,vector<T>& Er,matrix<T>& Eves,int Enu);

template<typename T>
void dstemr(vector<T>& diag, vector<T>& subdiag,int dim, vector<T>& E ,matrix<T>& Evecs, int Enu);

template<typename T,typename C>
void zheevr(matrix<C>& H,vector<T>& Er,matrix<C>& Evecs,int Enu);

template <typename T>
T* ctof(const matrix<T>& in);

template <typename T>
void ftoc(T *in, matrix<T>& out);

/* lapack for  c only not c++ */

template<typename T>
void dsyevr1v(T **H,T *Er,int n,T **Evecs,int Enu);

template<typename T>
void dsyevr1(T **H,T *Er,int n,int Enu);

template <typename T>
T* ctof1(T **in,int n);

template <typename T>
void ftoc1(T *in, T **out,int n);

template <typename T> 
T** AllocateDA1(  int nRows, int nCols);

template <typename T>
void FreeDA(T **dArray,int rows,int cols);





/*  other functions
**************************************
*/
template<typename T>
void v_sort(vector<T>& v,matrix<T>& M);

template<typename T>
void cv_sort(T* x,int n);

void i_time(double clk);

int factoriel(int n);
int Tarkib(int N,int n);



/*    class function body definitions 
***************************************
*/
template<typename T>
matrix<T>::matrix(int row ,int col , T t ): v(row , vector<T>(col , t)){};

template<typename T>
matrix<T>::matrix(const matrix<T>& A ): v(A.v){};

template<typename T>
matrix<T>::~matrix(){};

template<typename T>
void matrix<T>::reserve(int height,int lenth ) 
{
 this->v.reserve(height);
 for(int i=0;i<height;i++) this->v[i].reserve(lenth);
}


template<typename T>
matrix<T>& matrix<T>::operator= (const matrix<T>& A)
{
	if(&v != &A.v)
	{
		vector < vector<T> > u;
		v.swap(u) ;
		v=A.v;
	}
	return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator*= (const matrix<T>& A )
{
	int rows = A.rows();
	int cols =A.cols();
	for(int i=0;i<rows;i++)
	 {
		 for(int j=0;j<cols;j++)
		 {
			 for(int s=0;s<rows;s++) v[i][j] += v[i][s]*A.v[s][j];
		 }

	 }
	return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator+= (const matrix<T>& A)
{
	int rows = A.rows();
	int cols =A.cols();
	for(int i=0;i<rows;i++)
	 {
		 for(int j=0;j<cols;j++)
		 {
			 v[i][j] += A.v[i][j];
		 }
	 }
	return *this;
}

template<typename T>
matrix<T>& matrix<T>::operator-= (const matrix<T>& A)
{
	int rows = A.rows();
	int cols =A.cols();
	for(int i=0;i<rows;i++)
	 {
		 for(int j=0;j<cols;j++)
		 {
			 v[i][j] -= A.v[i][j];
		 }
	 }
	return *this;
}

template <typename T>
vector<T>& matrix<T>::operator[] (int i)
{
  return v[i];
}

template <typename T> 
const vector<T>& matrix<T>::operator[]  (int i) const
{
 return v[i];
}


template <typename T>
int matrix<T>::rows() const 
{ 	return v.size();}

template <typename T>
int matrix<T>::cols() const 
{	return v[0].size();}


/* non class member body operator/functions  definitions
******************************************************
*/


/*                  Arithmethic operators
******************************************************
*/
template <typename T>
 matrix<T> operator+ (const matrix<T>& A , const matrix<T>& B)
{

	matrix<T> temp(A);
	temp += B;
	return temp;
}

template< typename T>
 matrix<T> operator* (const matrix<T>& A , const matrix<T>& B)
{
	matrix<T> temp(A); 
	temp *= B;
	return temp;
}

template<typename T>
 vector<T> operator* (const matrix<T>& A , const vector<T>& V)
{
	T tsum = T(0);
	int col = A.cols();
	vector<T> vtemp;
	for(int i=0;i<A.rows();i++)
    	{
         for(int j=0;j<col;j++)  tsum+= A[i][j]*V[j];
		  vtemp.push_back(tsum);
		  tsum = T(0);
		}
	return vtemp;
}


template<typename T, typename C>
 matrix<C> operator* (const T& a , const matrix<C>& A)
 { 
	 int row = A.rows();
	 int col = A.cols();
	 matrix<C> temp(A);

	 for(int i=0;i< row ; i++)
		  for(int j=0;j<col ;j++) temp[i][j] = a * A[i][j];

	 return temp;
 }

template<typename T>
vector<T> operator- (const vector<T>& A ,const vector<T>& B )
{
	int sz = A.size();
	vector<T> temp(sz,A[0]); 
	for(int i=0 ;i<sz;i++) temp[i] = A[i]-B[i];
	return temp;
}


template<typename T>
vector<T> operator+ (const vector<T>& A ,const vector<T>& B )
{
        int sz = A.size();
        vector<T> temp(sz,A[0]);
        for(int i=0 ;i<sz;i++) temp[i] = A[i]+B[i];
        return temp;
}



template<typename T>
T operator* (const vector<T>& A ,const vector<T>& B )
{
	int sz = A.size();
        T temp = T(0);
	for(int i=0;i< sz; i++) temp +=  conj(A[i])*B[i];
        return temp;      
}

template<typename T>
matrix<T> operator- (const matrix<T>& A , const T& b)
{
        matrix<T> temp(A);
        int sz = A[0].size();
        for(int i=0 ; i<sz ; i++) temp[i][i] -= b;
        return temp; 
}

template<typename T>
vector<T> operator* (const T& a ,const vector<T>& A )
{
 int sz = A.size();
 vector<T> temp(sz , T(0));
 for(int i=0;i<sz ;i++) temp[i] = a * A[i];
 return temp;
}

template<typename T>
vector<T> conj(vector<T>& A)
{
vector<T> temp(A);

for(int i=0;i<temp.size();i++)
   temp[i]=conj(temp[i]);

return temp;
}

template<typename T>
vector<T> site(int i, int n, T t)
{
 vector<T> temp(n,T(0));
 temp[i] = t;
 return temp;
}


/*                  STD input output operators
******************************************************
*/
 template<typename T>
 ostream& operator<< (ostream& output , matrix<T>& A)
 {
	 int row = A.rows();
	 int col = A.cols();
         output.precision(5);
	 for(int i=0;i<row;i++)
	  {
	    for(int j=0;j<col;j++)
		output<< fixed <<A[i][j]<<" ";
                output<<endl;		
	  }
}

 template<typename T>
 ostream& operator<< (ostream& output , vector<T>& V)
 {
	 for(int i=0;i<V.size();i++)
	 {
             output<<V[i]<<" ";
//           output<<endl;
      
         }
         output<<endl;
	 return output;
 }

template<typename T>
 istream& operator>> (istream& input , matrix<T>& A)
 {
	 for(int i=0;i<A.rows();i++)
	  {
	    for(int j=0;j<A.cols();j++)
				input >> A[i][j];
      }
	 return input;
 }


/*  Lapack extern function usage 
*************************************
*/


template <typename T>
T* ctof(const matrix<T>& in)
{
  int rows=in.rows();
  int cols=in.cols();

  T *out;
  int i, j;

  out = new T[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}
 

template <typename T>
void ftoc(T *in, matrix<T>& out)
{
  int i, j;
  int rows=out.rows();
  int cols=out.cols();
  for (i=0; i<rows; i++){ 
    for (j=0; j<cols; j++){ 
      out[i][j] = in[i+j*cols];
    }
  }
}



template<typename T>
void dgeev(matrix<T>& H, vector<T>& Er, vector<T>& Ei, matrix<T>& Evecs)
{
  int n=H.cols();
  char jobvl, jobvr;
  int lda,  ldvl, ldvr, lwork, info;
  T *a, *vl, *vr, *work;
  T Er_temp[n];
  T Ei_temp[n];
  jobvl = 'N';    
  jobvr = 'V';
  lda = n;
  a = ctof(H);
  ldvl = n;
  vl = new T[n*n];
  ldvr = n;
  vr = new T[n*n];
  work = new T[4*n];
  lwork = 4*n;

  dgeev_(&jobvl, &jobvr, &n, a, &lda, &Er[0], &Ei[0], vl, &ldvl, vr, &ldvr, work, &lwork, &info);

  ftoc(vr, Evecs);

  delete a;
  delete vl;
  delete vr;
  delete work;
}

template<typename T>
void dstemr(vector<T>& diag, vector<T>& subdiag,int dim ,vector<T>& E, vector<T>& Evecs, int Enu)
{

int n = dim;
char jobz,range;
bool tryrac;
int il,info,iu,ldz,nzc,liwork,lwork,m;
T vl,vu;
int *isuppz,*iwork;
T *work;//,*w,*z;

jobz = 'V';
range='A';
info=0;
ldz=n;
nzc=n;
tryrac = true ;
lwork = 18*n;
liwork = 10*n;

work = new T[lwork];
iwork = new int[liwork];
isuppz = new int[2*n];

//

dstemr_(&jobz, &range, &n, &diag[0], &subdiag[1], &vl, &vu, &il, &iu, &m, &E[0], &Evecs[0], &ldz, &nzc, isuppz, &tryrac, work, &lwork, iwork, &liwork, &info);

//

delete work;
delete iwork;
//delete w;
delete isuppz;

return;
}

template<typename T>
void dsyevr(matrix<T>& H,vector<T>& Er,matrix<T>& Evecs,int Enu)
{
int n = H.cols();
char jobz,range,uplo;
int il,iu,m,ldz,lda,lwork,liwork,info;
T vl,abstol,vu;
T *w,*z,*work,*a;
int *isuppz,*iwork;

jobz = 'V';
range= 'I';
uplo=  'U';
a=ctof(H);
lda= n;
il=1;iu=Enu;
abstol=0.0;
info=0;
ldz=n;
z=new T[n*n];
lwork=28*n;
liwork=12*n;
work = new T[lwork];
iwork = new int[liwork];
w =new T[n];
isuppz = new int[2*n];


dsyevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m,
   &Er[0],z,&ldz,isuppz,work,&lwork ,iwork,&liwork,&info);
ftoc(z,Evecs);


delete z;
delete work;
delete iwork;
delete w;
delete isuppz;

}
//

template<typename T,typename C>
void zheevr(matrix<C>& H,vector<T>& Er,vector<C>& Evecs,int Enu)
{

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
il=1;iu=Enu;
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


zheevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m,
   &Er[0],&Evecs[0],&ldz,isuppz,work,&lwork, rwork, &lrwork, iwork,&liwork,&info);

delete a;
delete z;
delete work;
delete iwork;
delete w;
delete isuppz;

}



//
template<typename T>
void dsyevr(matrix<T>& H,vector<T>& Er,int Enu)
{
int n = H.cols();
char jobz,range,uplo;
int il,iu,m,ldz,lda,lwork,liwork,info;
T vl,abstol,vu;
T *w,*z,*work,*a;
int *isuppz,*iwork;

jobz = 'N';
range= 'A';
uplo=  'U';
a=ctof(H);
lda= n;
il=1;iu=Enu;
abstol=0.0;
info=0;
ldz=n;
z=new T[n*n];
lwork=28*n;
liwork=12*n;
work = new T[lwork];
iwork = new int[liwork];
w =new T[n];
isuppz = new int[2*n];


dsyevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m,
   &Er[0],z,&ldz,isuppz,work,&lwork ,iwork,&liwork,&info);

delete z;
delete work;
delete iwork;
delete w;
delete isuppz;

}




template <typename T>
matrix<T> Inv(const matrix<T>& A)
{
  matrix<T> B(A);
  int N=A.rows();

  if(N==2){
    T delta=A[0][0]*A[1][1]-A[0][1]*A[1][0];
    B[0][0]=1.0/delta*A[1][1];
    B[0][1]=-1.0/delta*A[0][1];
    B[1][0]=-1.0/delta*A[1][0];
    B[1][1]=1.0/delta*A[0][0];
  }
  else{

    T *matrix;
    matrix = ctof(A);

    int *ipiv = new int[N];
    int info;
    complex<double> query_work[2];
    complex<double> *work = 0;
   
    int lwork = -1;

    zgetrf_(&N, &N, matrix, &N, ipiv, &info);  
    zgetri_(&N, matrix, &N, ipiv, &query_work[0],&lwork, &info);
    lwork=10;
    
    work = new complex<double>[lwork];
    zgetri_(&N, matrix, &N, ipiv, work, &lwork, &info);

    ftoc(matrix, B);

    delete matrix;
    delete work;
    delete ipiv;
  }
  return B;
}

/*  ctof1 ftoc1 for c only not c++ matrix class or vectors */

template <typename T>
T* ctof1(T **in,int n)
{
  int rows=n;
  int cols=n;

  T *out;
  int i, j;

  out = new T[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}

template <typename T>
void ftoc1(T *in, T **out,int n)
{
  int i, j;
  int rows=n;
  int cols=n;
  for (i=0; i<rows; i++){ 
    for (j=0; j<cols; j++){ 
      out[i][j] = in[i+j*cols];
    }
  }
}

template<typename T>
void dsyevr1v(T **H,T *Er,int n,T **Evecs,int Enu)
{
char jobz,range,uplo;
int il,iu,m,ldz,lda,lwork,liwork,info;
T vl,abstol,vu;
T *w,*z,*work,*a;
int *isuppz,*iwork;

jobz = 'V';
range= 'I';
uplo=  'U';
a=ctof1(H,n);
lda= n;
il=1;iu=Enu;
abstol=0.0;
info=0;
ldz=n;
z=new T[n*n];
lwork=28*n;
liwork=12*n;
work = new T[lwork];
iwork = new int[liwork];
w =new T[n];
isuppz = new int[2*n];


dsyevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m, 
   Er,z,&ldz,isuppz,work,&lwork ,iwork,&liwork,&info);

ftoc1(z,Evecs,n);

delete z;
delete work;
delete iwork;
delete w;
delete isuppz;
return;
}




template<typename T>
void dsyevr1(T **H,T *Er,int n,int Enu)
{
char jobz,range,uplo;
int il,iu,m,ldz,lda,lwork,liwork,info;
T vl,abstol,vu;
T *work,*a,*z;
int *isuppz,*iwork;

jobz = 'N';
range= 'I';
uplo=  'U';
a=ctof1(H,n);
lda= n;
il=1;iu=Enu;
abstol=0.0;
info=0;
ldz=n;
lwork=28*n;
liwork=12*n;
work = new T[lwork];
iwork = new int[liwork];

isuppz = new int[2*n];


dsyevr_(&jobz,&range,&uplo, &n, a,&lda, &vl,&vu,&il,&iu ,&abstol, &m, 
   Er,z,&ldz,isuppz,work,&lwork ,iwork,&liwork,&info);

delete work;
delete iwork;
delete isuppz;
return;
}

void i_time(double clk)
{
  int seconds,minuts,hours;
  seconds = (int) clk/CLOCKS_PER_SEC;

  minuts  = (int) seconds/60;
  seconds = (int) seconds%60;
  hours   = (int) minuts/60;
  minuts  = (int) minuts%60;

  printf("hours:%d  minuts:%d seconds:%d \n",hours,minuts,seconds);
}

int factoriel(int n)
{
int temp = 1;
for(int i=0;i<n;i++) temp*=(n-i); 
return temp;
}

int Tarkib(int N, int n)
{
return factoriel(N)/(factoriel(n)*factoriel(N-n));
}



 /*                      MAIN PROG                           */

