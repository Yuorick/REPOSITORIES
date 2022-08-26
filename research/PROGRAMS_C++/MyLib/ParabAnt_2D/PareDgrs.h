//---------------------------------------------------------------------------

#ifndef PareDgrsH
#define PareDgrsH
#include "SincDgr.h"
class TSincDgr;
class Comp;
class TURPolyLine;
class TPareDgrs
{
public:
 // длина волны
 double	 mLambda ;
 // угол сдвига парциальных диаграмм  радианы
double  mAngSdvig ;

 // массив диаграмм в паре
 TSincDgr mSincDgr;



 TPareDgrs() ;


// конструктор копирования
 TPareDgrs(const TPareDgrs &R) ;
 TPareDgrs operator=(TPareDgrs  R2) ;
 // парам констр
__fastcall  TPareDgrs(const double Lambda,double RadAngSdvig, TSincDgr SincDgr);

void  createGraphs_For_BearingFuncs( wchar_t *wchFoldName1);

double   calcFirstCoefTalorBearFunc_For_PartDiagr();

double   calcSecondDerivBearFunc_For_PartDiagr(const double tet) ;

double   calcThirdCoefTalorBearFunc_For_PartDiagr();

int  estimateARSMethPartDiagr( TComp *cmparrPartSZv
		  ,  double *pvalEstPartDiaRadAngTarg,  TComp *pcmpPartDiaKTarg
		   , double *pvalPartDiaDisp, double *pvalPsi2);

void  createMtrxCorrForMeasures(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp,  double *arrMtrxCorr );

/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

TComp  fncDgrF0(const double VAlRadTetta );

double fncModDgr0(const double VAlRadTetta );

double fncArgDgr0(const double VAlRadTetta );

double fncDerivModDgr0(const double VAlRadTetta );

TComp  fncDgrF1(const double VAlRadTetta );

double fncModDgr1(const double VAlRadTetta );

double fncArgDgr1(const double VAlRadTetta );

double fncDerivModDgr1(const double VAlRadTetta );
//
TComp  fncDgrF0(const double VAlRadTetta , TURPolyLine &plnModulGraph, TURPolyLine &plnArgGraph);

double fncModDgr0(const double VAlRadTetta , TURPolyLine  &plnGraph);

double fncArgDgr0(const double VAlRadTetta , TURPolyLine  &plnGraph);

double fncDerivModDgr0(const double VAlRadTetta , TURPolyLine  &plnGraph);

 double fncBearingCoeff(TURPolyLine &plnGraph);
/*


TComp  fncDgrF1(const double VAlRadTetta , TURPolyLine  plnGraph);

double fncModDgr1(const double VAlRadTetta , TURPolyLine  plnGraph);

double fncArgDgr1(const double VAlRadTetta , TURPolyLine  plnGraph);

double fncDerivModDgr1(const double VAlRadTetta , TURPolyLine  plnGraph); */
//

double fncBearingCoeff() ;

bool createGraphs_GuarantSystError_2Targs( wchar_t *wchFoldName1, double valRadEps, double valRo
  , double valRadDiap);

static double max_(double a, double b);

void imitateMeasureArr(TComp cmpK, double valEps, TComp *cmparrMeas);

double estimateARSM(TComp *cmparrMeas);

bool createGraphs_GuarantSystError_2Targs( wchar_t *wchFoldName1
  , TURPolyLine &plnModulGraph0, TURPolyLine &plnModulArg0, TURPolyLine &plnModulGraph1, TURPolyLine &plnModulArg1
  ,double valRo1, double valRo2, double valRadEps,  double valRadDiap);

};

#endif
