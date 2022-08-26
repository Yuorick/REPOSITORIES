//---------------------------------------------------------------------------

#ifndef TripleH
#define TripleH
#include "SincDgr.h"
class TSincDgr;
class Comp;
class TURPolyLine;
class TTriple
{
public:
 // длина волны
 double	 mLambda ;
 // угол сдвига парциальных диаграмм  радианы
double  marrRadAngSdvig[2] ;

 // массив диаграмм в веере
 TSincDgr marrDgr[3];



 TTriple() ;


// конструктор копирования
 TTriple(const TTriple &R) ;
 TTriple operator=(TTriple  R2) ;
 // парам констр
__fastcall  TTriple(const double Lambda,double *arrRadAngSdvig, TSincDgr *arrDgr);

int estimate( TComp *cmparrS, double *valEstAngTarg
, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

int  calculateGenAngs(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

 int  findRootsFgr(wchar_t *wchFileName,double *arr_b,  double *arrRoots);

double  fncFGr(double *arr_b, double valmu);

double  findRootMethChord_For_Fgr(double *arr_b, double  valX0, double valX1);

bool calcVect_b(TComp *cmparrS, double *arr_b);

bool calcMtrxCorrGenAngs(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp , double *arrMtrxCorrGenAngs );

void createMtrxCorrForMeasures(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, double *arrMtrxCorr );

void createSystErrGraphs( wchar_t *Fold) ;

double fnc_dTet_po_dDel1(double valRadAngTargCur, double valRadAngAntpCur);

bool calcTrueVect_b(double valRadAngTarg, double valRadAngAntp, double *arr_b);

double fnc_dFGr_po_dTet(double valRadTet, double *arr_b);

double fnc_dTet_po_dDel2(double valRadAngTargCur, double valRadAngAntpCur)   ;

double fnc_dTet_po_dDel3(double valRadAngTargCur, double valRadAngAntpCur);

double fnc_dTet_po_da1(double valRadAngTargCur, double valRadAngAntpCur) ;

double fnc_dTet_po_da2(double valRadAngTargCur, double valRadAngAntpCur);

double fnc_dTet_po_da3(double valRadAngTargCur, double valRadAngAntpCur) ;

double fnc_dFGr_po_dLambda(double valRadTet, double *arr_b);

double fnc_dTet_po_dLambda(double valRadAngTargCur, double valRadAngAntpCur);

static void checkPlaneHypot(wchar_t *wchFoldReport
  ,TURPolyLine  plnModulGraph0,TURPolyLine  plnArgGraph0,TURPolyLine  plnModulGraph1
  ,TURPolyLine  plnArgGraph1,TURPolyLine  plnModulGraph2,TURPolyLine  plnArgGraph2);

};


 double max_(double a, double b )
 {
	 if (a > b)
	 {
	  return a;
	 }
	 return b;
 }
 double min_(double a, double b )
 {
	 if (a < b)
	 {
	  return a;
	 }
	 return b;
 }

 double max_(double a, double b, double c )
 {
	double d =  max_(a,  b )  ;

	 return max_(d, c ) ;
 }

 double min_(double a, double b, double c )
 {
	double d =  min_( a, b );

	 return min_(d, c ) ;
 }

 void calc_dB0_po_dS(TComp *cmparrS, double val_b0, double val_b1, double *arrdB_po_dS);

 void calc_dB1_po_dS(TComp *cmparrS, double val_b0, double val_b1, double *arrdB_po_dS);
#endif
