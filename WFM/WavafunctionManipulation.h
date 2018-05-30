#include <iostream>
#include <vector>
#include "stdio.h"
#include "math.h"
#include "assert.h"
#include "stdlib.h"
#include <iterator>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <list>
#include "string.h"
#include <float.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sort_vector.h>



using namespace std;

#define PI  3.141592653589793
#define PI_Square_Root 1.772453850905515881919

#define Malloc  malloc
#define Matrix(M,x,y)   *(M.data+(y)*M.col_size+(x))
#define MatrixP(M,x,y)   (M.data+(y)*M.col_size+(x))
#define Sigmoid(x) (1.0/(1.0+exp(-x)))
#define D_Sigmoid(x) (x*(1.0-x))
#define  Square(x)  ((x)*(x))
#define  CoefficientAttanch(X,Y,Z,V){ \
                    if((X+Y+Z)==3)\
                    { if(X==1&&Y==1&&Z==1)\
                       { V=V*6.0;}\
                      if((X==2&&Y==1&&Z==0)||(X==1&&Y==2&&Z==0)||(X==2&&Y==0&&Z==1)||(X==0&&Y==2&&Z==1)||(X==1&&Y==0&&Z==2)||(X==0&&Y==1&&Z==2))\
                       { V=V*3.0;}\
                    }\
                    if((X+Y+Z)==2)\
                    {  if((X==1&&Y==1&&Z==0)||(X==1&&Y==0&&Z==1)||(X==0&&Y==1&&Z==1))\
                       { V=V*2.0;   }\
                    }\
                    if(((X+Y+Z)==4))\
                    {   if( (X==3&&Y==1)||(Y==3&&X==1)||(X==3&&Z==1)||(Y==3&&Z==1)||(Z==3&&X==1)||(Z==3&&Y==1)){ V=V*4.0; }\
                        if( (X==2&&Y==2)||(Z==2&&X==2)||(Y==2&&Z==2)){ V=V*6.0; }\
                        if( (X==2&&Y==1&&Z==1)||(Z==1&&X==1&&Y==2)||(Z==2&&X==1&&Y==1) ){ V=V*12.0; }\
                     }\
                  }

typedef  struct Matrix
{
	int col_size;
	int row_size;
	double*data;
}Matrix;


#define Matrix_Init(M) {M.data =(double*)Malloc(M.col_size*M.row_size*sizeof(double));}
#define Matrix_set_zero(M){ memset(M.data,0,M.col_size*M.row_size*sizeof(double));}//
#define Matrix_Copy(M,N){ \
						if((M.col_size == N.col_size)&&(M.row_size==N.row_size))\
							{		memcpy(M.data,N.data,M.row_size*M.col_size*sizeof(double));\
								}\
							}
#define Matrix_Free(M){free(M.data);}


typedef struct Coordinate{
 double X;
 double Y;
 double Z;
 unsigned int Index;
}Coordinate;

typedef struct STFtype{
int Center;
unsigned int functype;
std::vector<double>exps;
std::vector<double>contracts;

}STFtype;

  typedef struct gsl_Type
  {
    gsl_matrix * m;
    gsl_matrix * evec;
    gsl_eigen_symmv_workspace * w;
    gsl_vector *eval;
    gsl_permutation *perm;
    Matrix rotationM;
  }gsl_Type;



typedef struct GTFtype{

 int Center;
 unsigned int functype;
 double exps;

}GTFtype;
inline bool operator<(const GTFtype &a,const GTFtype &b){ if(a.Center!=b.Center){ return a.Center<b.Center;}else{ if(a.exps!=b.exps){ return a.exps<b.exps;}else{ return a.functype <b.functype; } } }
inline bool operator>(const GTFtype &a,const GTFtype &b){ if(a.Center!=b.Center){ return a.Center>b.Center;}else{ if(a.exps!=b.exps){ return a.exps>b.exps;}else{ return a.functype >b.functype; } } }
 typedef struct Operator
 {
   int OperatorType;
   double OperatorParameter;

 }Operator;

#define STYPE 0
#define PTYPE 1
#define DTYPE 2
#define FTYPE 3
#define S_T 0
#define X_T 1
#define Y_T 2
#define Z_T 3
#define XX_T 4
#define YY_T 5
#define ZZ_T 6
#define XY_T 7
#define XZ_T 8
#define YZ_T 9
#define XXX_T 10
#define YYY_T 11
#define ZZZ_T 12
#define XXY_T 13
#define XXZ_T 14
#define YYZ_T 15
#define XYY_T 16
#define XZZ_T 17
#define YZZ_T 18
#define XYZ_T 19
#define ZZZZ_T 20
#define YZZZ_T 21
#define  YYZZ_T 22
#define  YYYZ_T 23
#define  YYYY_T 24
#define  XZZZ_T 25
#define  XYZZ_T 26
#define  XYYZ_T 27
#define  XYYY_T 28
#define  XXZZ_T 29
#define  XXYZ_T 30
#define  XXYY_T 31
#define  XXXZ_T 32
#define  XXXY_T 33
#define  XXXX_T 34

#define  decipherType(ix,iy,iz,E){\
                        if(ix==0&&iy==0&&iz==0){E=S_T;}\
                        if(ix==1&&iy==0&&iz==0){E=X_T;}\
                        if(ix==0&&iy==1&&iz==0){E=Y_T;}\
                        if(ix==0&&iy==0&&iz==1){E=Z_T;}\
                        if(ix==2&&iy==0&&iz==0){E=XX_T;}\
                        if(ix==0&&iy==2&&iz==0){E=YY_T;}\
                        if(ix==0&&iy==0&&iz==2){E=ZZ_T;}\
                        if(ix==1&&iy==0&&iz==1){E=XZ_T;}\
                        if(ix==1&&iy==1&&iz==0){E=XY_T;}\
                        if(ix==0&&iy==1&&iz==1){E=YZ_T;}\
                        if(ix==3&&iy==0&&iz==0){E=XXX_T;}\
                        if(ix==0&&iy==3&&iz==0){E=YYY_T;}\
                        if(ix==0&&iy==0&&iz==3){E=ZZZ_T;}\
                        if(ix==2&&iy==1&&iz==0){E=XXY_T;}\
                        if(ix==2&&iy==0&&iz==1){E=XXZ_T;}\
                        if(ix==0&&iy==2&&iz==1){E=YYZ_T;}\
                        if(ix==1&&iy==2&&iz==0){E=XYY_T;}\
                        if(ix==1&&iy==0&&iz==2){E=XZZ_T;}\
                        if(ix==0&&iy==1&&iz==2){E=YZZ_T;}\
                        if(ix==1&&iy==1&&iz==1){E=XYZ_T;}\
                        if(ix==0&&iy==0&&iz==4){E=ZZZZ_T;}\
                        if(ix==0&&iy==1&&iz==3){E=YZZZ_T;}\
                        if(ix==0&&iy==2&&iz==2){E=YYZZ_T;}\
                        if(ix==0&&iy==3&&iz==1){E=YYYZ_T;}\
                        if(ix==0&&iy==4&&iz==0){E=YYYY_T;}\
                        if(ix==1&&iy==0&&iz==3){E=XZZZ_T;}\
                        if(ix==1&&iy==1&&iz==2){E=XYZZ_T;}\
                        if(ix==1&&iy==2&&iz==1){E=XYYZ_T;}\
                        if(ix==1&&iy==3&&iz==0){E=XYYY_T;}\
                        if(ix==2&&iy==0&&iz==2){E=XXZZ_T;}\
                        if(ix==2&&iy==1&&iz==1){E=XXYZ_T;}\
                        if(ix==2&&iy==2&&iz==0){E=XXYY_T;}\
                        if(ix==3&&iy==1&&iz==0){E=XXXY_T;}\
                        if(ix==4&&iy==0&&iz==0){E=XXXX_T;}\
                     }

#define  X_phi_ 1
#define  Y_phi_ (1<<6)
#define  Z_phi_ (1<<11)

#define   X_P  1
#define   Y_P  2
#define   Z_P  3
#define   XX_P 4
#define   YY_P 5
#define   ZZ_P 6
#define   XY_P 7
#define   XZ_P 8
#define   YZ_P 9


static int  typeTx[56]={0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,0,1,1,0,1,0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5};
static int  typeTy[56]={0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3,0,1,2,0,1,0};
static int  typeTz[56]={0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,1,0,2,2,1,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0,5,4,3,2,1,0,4,3,2,1,0,3,2,1,0,2,1,0,1,0,0};
static int  typeSize[8]={1,3,6,10,15,21,28,36};
static int  typeSizeBase[8]={0,1,4,10,20,35,56,84};


int ft (int i)
{
  int value=1;

    if(i==0)
        return 1;
    for(;i>1;i--)
    {
       value=value*i;
    }

  return value;
}

template <typename T>
std::vector<size_t> ordered(std::vector<T> const& values,bool ascending) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));

    if(ascending)
    {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] < values[b]; }
         );
    }else
    {
        std::sort(
            begin(indices), end(indices),
            [&](size_t a, size_t b) { return values[a] > values[b]; }
         );
    };
    return indices;
}

extern double overlapIntegral(unsigned int itype1,double exp1,unsigned int itype2,double exp2,Coordinate A ,Coordinate B);
extern int  RotateOperator(std::vector<GTFtype>&b2,std::vector< Coordinate >&AtomN,std::vector< double >&MO_Old,std::vector< double >&MO_New,Matrix rotationM);
extern 	void InitInvariantContext(gsl_Type * Context);
extern void FreeInvariantContext(gsl_Type * Context);
extern  void  InvariantDirection(std::vector<GTFtype>bs,std::vector<double> MO_i,std::vector< Coordinate > AtomN_i,std::vector<double> MO_o,std::vector< Coordinate > AtomN_o,gsl_Type * Context);
	

 






