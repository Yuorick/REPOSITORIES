//---------------------------------------------------------------------------

#ifndef Veer4H
#define Veer4H
#include "SimpleSumDgrm.h"
#include "SimpleDifDgrm.h"
#include "Comp.h"
class TComp;
class TSimpleSumDgrm ;
class TSimpleDifDgrm ;
class TVeer4
{
public:

 TSimpleSumDgrm  marrSumDgr[2];
// TSimpleDifDgrm  marrDifDgr[2];
 // ??????? ???? ??????? ?????????? ?????????
 double mTang;
 // ????? ???????????? ?????????? ?????????
 double	 marrScnDif[2];

 __fastcall  TVeer4() ;
// ??????????? ???????????
__fastcall  TVeer4 (const TVeer4 &R2) ;
 // ????? ??????
 __fastcall TVeer4(const double Tang,  double *arrScnDif, const double KSum
	   , const double ScnSum, const double Tension );
// ???????? ????????????
TVeer4   operator=(TVeer4  R2) ;

//--------------------------------------------------------------------------
int  _fastcall solv4(TComp *pSDifZv,TComp *pSSumZv, TComp *pZTarg, TComp *pZAnt, double *palfTrg, double *palfAnp ) ;

int _fastcall calc_vectG_and_mtrxH (TComp *pSDifZv,TComp *pSSumZv
		,  double alfTrg, double alfAnp, double* arr_FGreek, double* arr_dFGreek ) ;

void _fastcall calc_CoeffOtrag(TComp *pSDifZv,double alfTrg,double  alfAnp,double* valK11, double* valK12
	 ,double* valK21,double* valK22,double* arr_gradK11, double* arr_gradK12, double* arr_gradK21
	 ,double* arr_gradK22, double*arr_HessK11, double* arr_HessK12, double* arr_HessK21,double* arr_HessK22);


static void _fastcall calcFncs_Fi1_and_Fi2( TComp  cmpX1,  TComp cmpX2,TComp& cmpFi1,TComp & cmpFi2);

static void _fastcall calcVectrsGrad_Fi1_and_Fi2( TComp   cmpX1,  TComp  cmpX2
		 ,TComp* cmpArrGradFi1, TComp* cmpArrGradFi2) ;

static  void _fastcall calcMtrxHess_Fi1_and_Fi2( TComp  cmpX1, TComp  cmpX2
	   ,TComp* cmpArrHessFi1,TComp* cmpArrHessFi2) ;
};
#endif
