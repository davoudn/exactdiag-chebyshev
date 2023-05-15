#include <vector>
#include <complex>

using namespace std;
#define IM std::complex<double>(0.0,1.0)
#define RE std::complex<double>(1.0,0.0)
/*                                                   */



/*                green function                     */
template <typename T>
class Gfunct
{
public:
Gfunct(int orbit,int egridno);
~Gfunct();

int size();
std::vector<T>& operator() (int m,int m_);
Gfunct<T>& operator= (const Gfunct<T>& A);

vector< vector<T> > a;
int orbitalno;
};


template <typename T>
Gfunct<T>::Gfunct(int orbit,int egridno):a(orbit*orbit, vector<T>(egridno, T(0))) {orbitalno = orbit; };

template <typename T>
Gfunct<T>::~Gfunct(){};

template <typename T>
int Gfunct<T>::size()
{
 
return a[0].size();

}


template <typename T>
vector<T>& Gfunct<T>::operator() (int m,int m_)
{
  return a[m+orbitalno*m_];
}

template <typename T>
Gfunct<T>& Gfunct<T>::operator= (const Gfunct<T>& A)
{
 for(int j=0;j<a.size();j++)  a[j] = A.a[j];
 return *this;
}

template <typename T>
Gfunct<T> operator- (const Gfunct<T>& A,const Gfunct<T>& B)
{
Gfunct<T> temp(A.orbitalno,A.a[0].size());

for(int m=0;m<A.orbitalno;m++)
 for(int m_=0;m_<A.orbitalno;m_++)
    for(int i=0;i<A.a[0].size();i++)
       temp(m,m_)[i] = A.a[m+A.orbitalno*m_][i]-B.a[m+A.orbitalno*m_][i];
    
return temp;
}


/*    Fractional green function                      */

template <typename T,typename C>
class Fractional_GreenFunction{
public:
 Fractional_GreenFunction();
~Fractional_GreenFunction();
vector<T> a;
vector<T> b;
T E0;
double gl;
C operator() (C omega);
};


template <typename T , typename C >
Fractional_GreenFunction<T,C>::Fractional_GreenFunction(){};


template <typename T , typename C >
Fractional_GreenFunction<T,C>::~Fractional_GreenFunction(){};

template<typename T , typename C >
C Fractional_GreenFunction<T,C>::operator() (C z)
{
int n = b.size();
C temp0,temp1;
E0 = T(0);

{
temp0 =  pow((z*gl + (E0   - a[n-1])), -1.0)*pow(b[n-1],2.0);

 for(int i=n-1;i>1;i--)
   {
    temp1 = z*gl + (E0   - a[i-1]) - temp0;
    temp0 = pow(temp1, -1.0)*pow(b[i-1],2.0);
   }

temp0 = gl*pow(z*gl + (E0 - a[0]) - temp0,-1.0);
}

return temp0;
}


/*                omega mesh class                   */
template <typename T, typename C>
class wMesh
{
public:
wMesh(T beta, int sz);
wMesh(T width, int meshsize, T eta);
~wMesh();
int size();

vector<C> a;
C operator[](int i);
};

template <typename T, typename C>
wMesh<T,C>::wMesh(T beta, int sz):a(sz,C(0))
{
  for(int i=0;i<sz;i++) a[i] = acos(-1.0)*(2.0*(T)i +1.0)*pow(beta,-1)*IM;
};

template <typename T, typename C>
wMesh<T,C>::wMesh(T width, int meshsize, T eta):a(meshsize,C(0))
{
  for(int i=0;i<meshsize;i++) a[i] = width*pow((T)meshsize,-1.0)*((T)(i-0.5*meshsize)*RE-eta*IM);
};



template <typename T, typename C>
wMesh<T,C>::~wMesh(){};

template <typename T, typename C>
C wMesh<T,C>::operator[](int i)
{
return a[i];
}

template <typename T, typename C>
int wMesh<T,C>::size()
{
return a.size();
}


template<typename T , typename C >
C Frac_Green(vector<C>& a,vector<C>& b, T E0, T gl,C z) 
{
int n = b.size();
C temp0,temp1;

{
temp0 =  pow((z +gl*(E0   - a[n-1])), -1.0)*pow(b[n-1],2.0);

 for(int i=n-1;i>1;i--)
   {
    temp1 = z + gl*(E0   - a[i-1]) - temp0;
    temp0 = pow(temp1, -1.0)*pow(b[i-1],2.0);
   }

temp0 = pow(z + gl*(E0 - a[0]) - temp0,-1.0);
}

return temp0;
}









