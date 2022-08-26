//---------------------------------------------------------------------------

#ifndef MatrixProccessH
#define MatrixProccessH
//---------------------------------------------------------------------------
void MtrxMultMatrx(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez) ;
void MtrxMultMatrxTransp(double *parrA,int nRowsA, int nColsA, double * parrB,int nRowsB, double *parrRez) ;
void MtrxTranspMultMatrx(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez) ;
void MtrxSumMatrx(double *parrA, double * parrB,int nRows, int nCols, double *parrRez) ;
void MtrxMinusMatrx(double *parrA, double * parrB,int nRows, int nCols, double *parrRez);
void MatrxMultScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez);
void MatrxDivideScalar(double *parrA, int nRows, int nCols, double valScal,double *parrRez);
void MatrTransp(double *parrA, int nRows, int nCols, double *parrRez);
void MtrxMultMatrx_MultMatrx(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut)  ;
void MtrxMultMatrx_MultMatrxTransp(double *parrInp0,double *parrInp1,double *parrInp2,int nDim, double *parrOut);
void OuterProduct(double *pVect0 , double *pVect1, double *pVectRez) ;
double ScalProduct(double *pVect0 , double *pVect1, const int len) ;
double Norm3( double *arrA) ;
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
bool   InverseMtrx5(double *arrInp, double *arrOut);
double   calcDet5( double *arrInp);
bool   InverseMtrx( double *arrInp, const int nDimen, double *arrOut) ;
double   calcDet6( double *arrInp);
bool   InverseMtrx6(double *arrInp, double *arrOut);
double   calcDet7( double *arrInp);
bool   InverseMtrx7(double *arrInp, double *arrOut);
double   calcDet8( double *arrInp);
bool   InverseMtrx8(double *arrInp, double *arrOut);
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
void   NormalizeVect(double *p, const int lenp);
double _sqrt_(const double a);




void MtrxMultMatrx(long double *parrA,int nRowsA, int nColsA, long double * parrB,int nColsB, long double *parrRez) ;
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






#endif
