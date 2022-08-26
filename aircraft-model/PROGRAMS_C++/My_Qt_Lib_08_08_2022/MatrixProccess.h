//---------------------------------------------------------------------------

#ifndef MatrixProccessH
#define MatrixProccessH
class TComp;
//---------------------------------------------------------------------------
void MtrxMultMatrx( double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez) ;
void MtrxMultMatrxTransp(double *parrA,int nRowsA, int nColsA, double * parrB,int nRowsB, double *parrRez) ;
void MtrxTranspMultMatrx(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez) ;
void MtrxSumMatrx(double *parrA, double * parrB,int nRows, int nCols, double *parrRez) ;
void MtrxMinusMatrx(double *parrA, double * parrB,int nRows, int nCols, double *parrRez);
void MatrxMultScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez);
void MatrxDivideScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez);
void MatrTransp( double *parrA, int nRows, int nCols, double *parrRez);
void MtrxMultMatrx_MultMatrx(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut)  ;
void MtrxMultMatrx_MultMatrxTransp(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut);
void MtrxMultMatrx_MultMatrxTransp(double *parrInp0,double *parrInp1,double *parrInp2,int nRows,int nCols, double *parrOut);

void calcF_D_FTransp_(double *parrFInp,int nRows,int nCols,double *parrDInp, double *parrOut);

void calcF_Mult_FTransp(double *parrA,int nRowsA, int nColsA, double *parrOut);

void OuterProduct(double *pVect0 , double *pVect1, double *pVectRez) ;

double ScalProduct(double *pVect0 , double *pVect1, const int len) ;

double Norm3( double *arrA) ;

double Norm3_A_Minus_B( double *arrA,double *arrB);

double OuterProduct_2(double *pVect0 , double *pVect1);

double MinDoubleArray(double *parrDoubleInp, const int lenArrayInp, int *pNumArgMin) ;

double MaxDoubleArray(double *parrDoubleInp, const int lenArrayInp, int *pNumArgMax) ;

int MaxIntArray(int *parrDoubleInp, const int lenArrayInp, int *pNumArgMax);

bool   InverseMtrx3(double *arrInp, double *arrOut);

double   calcDet3( double *arrInp);

bool   SolvFrobenius6(double *arrInpMtrxA00,double *arrInpMtrxA01
	   ,double *arrInpMtrxA10,double *arrInpMtrxA11, double *arrInpVectB,double *arrOutX);

 // вычисление YT * D * Y для матрицы 3х3
double   calcYT_D_Y(double *arrY, double *arrD );

int CalcProperVectors_And_Numbers_R3(double * arrKInp,double *arrV , double *arrLamb)  ;

void   swap(double *a0, double *a1)  ;

void   swap_vect(double *pa0, double *pa1, const int len) ;

double   NormVect3(double *p)  ;

double   NormVect2(double *p);

void   NormalizeVect3(double *p);

int  CalcProperVectors2(double * arrKInp,double *arrV , double *arrLamb)  ;

bool   SolvLinEq2(double *arrA, double *arrB,double *arrX);

bool   InverseMtrx2(double *arrA, double *arrOut);

double NormVect(double *arr, const int lenarr);

//bool   InverseFrobMtrx4(double *arrInp, double *arrOut);
//bool   InverseFrobMtrx5(double *arrInp, double *arrOut);
void buildAddMinor(double *parrMtrxInp, const int idim
	, const int iRowElem,const int iColElem, double *parrMinor) ;
double   calcDet4( double *arrInp);

bool   InverseMtrx4(double *arrInp, double *arrOut);

bool   InverseMtrx5(double *arrInp, double *a4rrOut);

double   calcDet5( double *arrInp);

bool   InverseMtrx( double *arrInp, const int nDimen, double *arrOut) ;

double   calcDet6( double *arrInp);

bool   InverseMtrx6(double *arrInp, double *arrOut);

double   calcDet7( double *arrInp);

bool   InverseMtrx7(double *arrInp, double *arrOut);

double   calcDet8( double *arrInp);

double   calcDet9( double *arrInp);

bool   InverseMtrx8(double *arrInp, double *arrOut);

bool   InverseMtrx9(double *arrInp, double *arrOut);

bool GaussMeth(double *parrA,int nDimA,  double * parrB, double *parrRez);

int fncFindNumVedushiaStroka(double *arrAUn, int  nDimA, int i0);

void Swap (int *iarrNums, int iNumCols, int j,int  i);

void Swap (double *parrA, int iNumCols, int j,int  i);

bool fncPodstanovka( double *arrAUn, int nDimA, int i, int j);

void formMatrxE(const int IdimE, double *arrOut);

void createOrthogBasis_dim3 (double *arrVectMissV, double *arrTargV, double *arrBasis);

double  calcYT_D_Y(double *arrY, double *arrD, const int IDim );

double  calcAngBetweenVect (double *arr1, double *arr2,const  int len);

double MinColDoubleMtrx(double *parrDoubleMtrx, const int nRows, const int nCols
	 , const int NumCol, int *pNumArgMin);

void calcF_D_FTransp(double *parrFInp,double *parrDInp,int nDim, double *parrOut);
void flipArray(double *parr, const int len);
void flipTwoDimArray(double *parr, const int lenrows, const int lencols);
bool NormalizeVect(double *p, const int LEnp);
void caclMtrxPolyn2( double a, double b, double c
                 ,double *arrA0, const int dimA, double *arrARez);

double Sp( double *arrA0, const int dimA);
bool solvCharactEq_5(double *arrA, TComp *cmparrRoots, const double VAlMin, const double VAlMax);
void createCrilovMtrx(double *arrA, const int dimA, double *arrx, double *arrCrilov);
bool solvCharactEq_4(double *arrA, TComp *cmparrRoots);

double   calcDet_A_minus_LambdaE_D3( double *arrInp, const double VAlLamb);

void fillE(double *arrE, const int dim);
bool solvCharactEq_4_New(double *arrA, TComp *cmparrRoots);
void calcCharactPolinomCoeffs_3dim(double *arrKInp,double * arrq);
void calcPolinomCoeffs_DopMinor(double *arrKInp,double * arrq);
void excludeRow(double *arrA,int qRows, int qCols, int numRow, double *arrARez);
void excludeCol(double *arrA,int qRows, int qCols, int numCol, double *arrARez);
void swapCols(double *arrA,int qRows, int qCols, int numCol0,int numCol1);
double NormSquareVect(double *arr, const int lenarr);

void changeCol(double *arrC, const int numRow,  const int numCol
               ,const int numChangedCol ,double *arrNewCol);

int getSegmentNum(double *arrX, const int len,const double x);

double NormSquareVect(const double *arr, const int lenarr);

double ScalProduct(const double *pVect0 ,const  double *pVect1, const int len);

bool InverseMtrx_GaussMeth(double *parrA,int nDimA, double *parrInvA);

void getMainMInor(double *arrInp, const int nDimen, const int nMinorDimention, double *arrT);

bool   isPositiveDefinite( double *arrInp, const int nDimen);


void MtrxMultMatrx( long double *parrA,int nRowsA, int nColsA, long double * parrB,int nColsB, long double *parrRez) ;
void MtrxMultMatrxTransp(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nRowsB, long double *parrRez) ;
void MtrxTranspMultMatrx(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nColsB, long double *parrRez) ;
void MtrxSumMatrx(long double *parrA, long double * parrB,int nRows, int nCols, long double *parrRez) ;
void MtrxMinusMatrx(long double *parrA, long double * parrB,int nRows, int nCols, long double *parrRez);
void MatrxMultScalar(long double *parrA, int nRows, int nCols, long double valScal,long double *parrRez);
void MatrxDivideScalar(long double *parrA, int nRows, int nCols, long double valScal,long double *parrRez);
void MatrTransp(long double *parrA, int nRows, int nCols, long double *parrRez);
 void OuterProduct(long double *pVect0 , long double *pVect1, long double *pVectRez) ;
long double ScalProduct(long double *pVect0 , long double *pVect1, const int len) ;
long double Norm3( long double *arrA) ;
long double   NormVect2(long double *p);
long double OuterProduct_2(long double *pVect0 , long double *pVect1);
long double MinDoubleArray(long double *parrDoubleInp, const int lenArrayInp, int *pNumArgMin) ;
long double MaxDoubleArray(long double *parrDoubleInp, const int lenArrayInp, int *pNumArgMax) ;
long double   calcDet3( long double *arrInp);
 bool   InverseMtrx3(long double *arrInp, long double *arrOut);
 long double   calcDet3( long double *arrInp);
 // вычисление YT * D * Y для матрицы 3х3
long double   calcYT_D_Y(long double *arrY, long double *arrD );
int CalcProperVectors_And_Numbers_R3(long double * arrKInp,long double *arrV , long double *arrLamb)  ;
void   swap(long double *a0, long double *a1)  ;
void   swap_vect(long double *pa0, long double *pa1, const int len) ;
long double   NormVect3(long double *p)  ;
void   NormalizeVect3(long double *p);
int  CalcProperVectors2(long double * arrKInp,long double *arrV , long double *arrLamb)  ;
bool   SolvLinEq2(long double *arrA, long double *arrB,long double *arrX);
long double NormVect(long double *arr, const int lenarr);

void buildAddMinor(long double *parrMtrxInp, const int idim
	, const int iRowElem,const int iColElem, long double *parrMinor) ;
long double   calcDet4( long double *arrInp);
bool   InverseMtrx4(long double *arrInp, long double *arrOut);
bool   InverseMtrx5(long double *arrInp, long double *arrOut);
long double   calcDet5( long double *arrInp);
bool   InverseMtrx(long double *arrInp, const int nDimen, long double *arrOut);

bool   InverseMtrx(long double *arrInp, const int nDimen, long double *arrOut) ;
long double    calcDet6( long double  *arrInp);
bool   InverseMtrx6(long double  *arrInp, long double  *arrOut);

long double    calcDet7( long double  *arrInp);
long double    calcDet8( long double  *arrInp)  ;
bool   InverseMtrx7(long double *arrInp, long double *arrOut);

bool   InverseMtrx8(long double *arrInp, long double *arrOut);
bool   InverseMtrx2(long double *arrA,long  double *arrOut);

bool GaussMeth(long double *parrA,int nDimA,  long double * parrB, long double *parrRez);
int fncFindNumVedushiaStroka(long double *arrAUn, int  nDimA, int i0);
void Swap (long double *parrA, int iNumCols, int j,int  i);
bool fncPodstanovka( long double *arrAUn, int nDimA, int i, int j);
void MtrxMultMatrx_MultMatrx(long double *parrInp0,long double *parrInp1,long double *parrInp2,int nDim,long  double *parrOut);
void MtrxMultMatrx_MultMatrxTransp(long double *parrInp0,long double *parrInp1,long double *parrInp2,int nDim,long  double *parrOut);



//---------------------------------------------------------------------------
void MtrxMultMatrx(float *parrA,int nRowsA, int nColsA, float * parrB,int nColsB, float *parrRez) ;
void MtrxMultMatrxTransp(float *parrA,int nRowsA, int nColsA, float * parrB,int nRowsB, float *parrRez) ;
void MtrxTranspMultMatrx(float *parrA,int nRowsA, int nColsA, float * parrB,int nColsB, float *parrRez) ;
void MtrxSumMatrx(float *parrA, float * parrB,int nRows, int nCols, float *parrRez) ;
void MtrxMinusMatrx(float *parrA, float * parrB,int nRows, int nCols, float *parrRez);
void MatrxMultScalar(float *parrA, int nRows, int nCols, float valScal,float *parrRez);
void MatrxDivideScalar(float *parrA, int nRows, int nCols, float valScal,float *parrRez);
void MatrTransp(float *parrA, int nRows, int nCols, float *parrRez);
void MtrxMultMatrx_MultMatrx(float *parrInp0,float *parrInp1,float *parrInp2,int nDim, float *parrOut)  ;
void MtrxMultMatrx_MultMatrxTransp(float *parrInp0,float *parrInp1,float *parrInp2,int nDim, float *parrOut);
void OuterProduct(float *pVect0 , float *pVect1, float *pVectRez) ;
float ScalProduct(float *pVect0 , float *pVect1, const int len) ;
float Norm3( float *arrA) ;
float OuterProduct_2(float *pVect0 , float *pVect1);
float MinfloatArray(float *parrfloatInp, const int lenArrayInp, int *pNumArgMin) ;
float MaxfloatArray(float *parrfloatInp, const int lenArrayInp, int *pNumArgMax) ;

bool   InverseMtrx3(float *arrInp, float *arrOut);
float   calcDet3( float *arrInp);

 // вычисление YT * D * Y для матрицы 3х3
float   calcYT_D_Y(float *arrY, float *arrD );
int CalcProperVectors_And_Numbers_R3(float * arrKInp,float *arrV , float *arrLamb)  ;
void   swap(float *a0, float *a1)  ;
void   swap_vect(float *pa0, float *pa1, const int len) ;
float   NormVect3(float *p)  ;
float   NormVect2(float *p);
void   NormalizeVect3(float *p);
int  CalcProperVectors2(float * arrKInp,float *arrV , float *arrLamb)  ;
bool   SolvLinEq2(float *arrA, float *arrB,float *arrX);
bool   InverseMtrx2(float *arrA, float *arrOut);
float NormVect(float *arr, const int lenarr);
//bool   InverseFrobMtrx4(float *arrInp, float *arrOut);
//bool   InverseFrobMtrx5(float *arrInp, float *arrOut);
void buildAddMinor(float *parrMtrxInp, const int idim
    , const int iRowElem,const int iColElem, float *parrMinor) ;
float   calcDet4( float *arrInp);
bool   InverseMtrx4(float *arrInp, float *arrOut);
bool   InverseMtrx5(float *arrInp, float *arrOut);
float   calcDet5( float *arrInp);
bool   InverseMtrx( float *arrInp, const int nDimen, float *arrOut) ;
float   calcDet6( float *arrInp);
bool   InverseMtrx6(float *arrInp, float *arrOut);
float   calcDet7( float *arrInp);
bool   InverseMtrx7(float *arrInp, float *arrOut);
float   calcDet8( float *arrInp);
bool   InverseMtrx8(float *arrInp, float *arrOut);
bool GaussMeth(float *parrA,int nDimA,  float * parrB, float *parrRez);
int fncFindNumVedushiaStroka(float *arrAUn, int  nDimA, int i0);
void Swap (int *iarrNums, int iNumCols, int j,int  i);
void Swap (float *parrA, int iNumCols, int j,int  i);
bool fncPodstanovka( float *arrAUn, int nDimA, int i, int j);
void formMatrxE(const int IdimE, float *arrOut);
void createOrthogBasis_dim3 (float *arrVectMissV, float *arrTargV, float *arrBasis);
float  calcYT_D_Y(float *arrY, float *arrD, const int IDim );
float  calcAngBetweenVect (float *arr1, float *arr2,const  int len);
float MinColfloatMtrx(float *parrfloatMtrx, const int nRows, const int nCols
     , const int NumCol, int *pNumArgMin);
void calcF_D_FTransp(float *parrFInp,float *parrDInp,int nDim, float *parrOut);
void flipArray(float *parr, const int len);
void flipTwoDimArray(float *parr, const int lenrows, const int lencols);
bool NormalizeVect(float *p, const int LEnp);
void caclMtrxPolyn2( float a, float b, float c
                 ,float *arrA0, const int dimA, float *arrARez);

float Sp( float *arrA0, const int dimA);

void formMatrxE(const int IdimE, float *arrOut);

//---------------------------------------------------------------
template <class Obj>
int SIGNUM__(Obj &obj)
{
    if (obj ==0)
    {
        return 0;
    }
    return (obj > 0)?1:-1;
}



#endif
