//---------------------------------------------------------------------------

#ifndef CompH
#define CompH
class TComp
  {

	 public:
	 double m_Re; // действительная часть
	 double m_Im ; // мнимая часть
	 TComp();
	 // конструктор копирования
	 TComp  (const TComp &R) ;
	  // оператор присваивания
     TComp &operator=(const TComp  &R);
	 // парам констр
	 TComp (  const double x,const double y) ;

	TComp  operator +(const TComp &cmp0);
	TComp  operator *(const TComp &cmp0);
    TComp  operator /(const TComp &cmp0);
	const TComp & operator +=(const TComp &cmp0);
	const TComp & operator -=(const TComp &cmp0);
	const TComp & operator *=(const TComp &cmp0) ;
	TComp  operator -(const TComp &cmp0);

	 double modul();
	 double phase();
	 TComp Sqrt_();
	 TComp Ln();
	 TComp Sopr();
	 void root3( TComp *cmparrRoots) ;

	 static void HermiteSoprMatr(TComp*parrA, int nRows, int nCols, TComp*parrRez);
	 static bool   InverseCmpMtrx2( TComp*arrA, TComp*arrOut);
	 static TComp scalProd(TComp *cmparrA, TComp *cmparrB, const int LEn);
	 static bool angBetveenComps(TComp cmpInp1, TComp cmpInp2, double *pvalOut);
	 static bool projectVect_3Dim(TComp *cmparrBasis0, TComp *cmparrBasis1, TComp *cmparrInp
  ,TComp *pcmpAlf0, TComp *pcmpAlf1, TComp *pcmparrPerp);


  };
TComp exp_(const TComp cmp);
TComp exp_(const double fi);
void MtrxMultMatrx(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nColsB, TComp*parrRez);
void MtrxMultMatrxTransp(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nRowsB, TComp*parrRez) ;
void MtrxTranspMultMatrx(TComp*parrA,int nRowsA, int nColsA, TComp* parrB,int nColsB, TComp*parrRez) ;
void MtrxSumMatrx(TComp*parrA, TComp* parrB,int nRows, int nCols, TComp*parrRez) ;
void MtrxMinusMatrx(TComp*parrA, TComp* parrB,int nRows, int nCols, TComp*parrRez);
void MatrxMultScalar(TComp*parrA, int nRows, int nCols, TComp valScal,TComp*parrRez);
void MatrxDivideScalar(TComp*parrA, int nRows, int nCols, TComp valScal,TComp*parrRez);
void MatrTransp(TComp*parrA, int nRows, int nCols, TComp*parrRez);
void MatrRe(TComp*parrA, int nRows, int nCols, double *parrRez);
void MatrIm(TComp*parrA, int nRows, int nCols, double *parrRez);
void   swap(TComp *a0, TComp *a1);
 bool   InverseMtrx3(TComp *arrInp, TComp *arrOut);
 TComp   calcDet3( TComp *arrInp);
 bool   InverseMtrx4(TComp *arrInp, TComp *arrOut);
 TComp   calcDet4( TComp *arrInp);
 void buildAddMinor(TComp *parrMtrxInp, const int ndim
	, const int iRowElem,const int jColElem, TComp *parrMinor);
 void createHankel( int lenS, TComp *cmparrS, TComp *cmparrHankel );
 void HermiteSoprMatr(TComp*parrA, int nRows, int nCols, TComp*parrRez);
 bool   InverseMtrx2( TComp*arrA, TComp*arrOut);
 bool fncLinFrac(TComp cmpa, TComp cmpb, TComp cmpc, TComp cmpd, TComp cmpZ, TComp *pcmpRez);
 bool  findCircleParams (TComp cmp0,TComp  cmp1,TComp  cmp2
	, TComp  *pcmpCentre, double *pvalRadius);
bool   SolvLinEq2__(double *arrA, double *arrB,double *arrX);

 TComp transfTrigonForm(double valMod, double valArg);

 TComp scalProd(TComp *cmparr0,TComp *cmparr1, int lenarr);

 TComp  vectNorm(TComp *cmparr0, int lenarr) ;
 double vectNorm_(TComp *cmparr0, int lenarr);

 bool   InverseMtrx( TComp *arrInp, const int nDimen,  TComp *arrOut);



 class TCompLong//:public TObject
  {

	 public:
	 long double m_Re; // действительная часть
	 long double m_Im ; // мнимая часть
	 TCompLong();
	 // конструктор копирования
	 TCompLong  (const TCompLong &R) ;
	  // оператор присваивания
	 TCompLong operator=(TCompLong  R);
	 // парам констр
	 TCompLong (  const double x,const double y) ;

     double modul();

  };
  #endif
