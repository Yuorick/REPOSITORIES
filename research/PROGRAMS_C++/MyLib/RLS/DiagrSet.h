//---------------------------------------------------------------------------

#ifndef DiagrSetH
#define DiagrSetH

class  TComp;
class TDiagr;
class TDiagrSet
{
public:



 // массив диаграмм
 TDiagr *mpDiagr;
 // к-во диаграмм
 int mNumDiagr;

  ~TDiagrSet();
   TDiagrSet();
// конструктор копирования
 TDiagrSet(const TDiagrSet &R) ;
  // оператор присваивания
 TDiagrSet operator=(TDiagrSet  R2) ;
 // парам конструктор
 TDiagrSet(const int NumDiagr, double *arrAlfaDiagr);
  // парам конструктор
 TDiagrSet(const int NumDiagr, TDiagr *pDiagr);

 //-----------------------------------------------------------------------------------------------------------------
  int  _fastcall MMP_GroupRelax(TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt, double *palfTrg, double *palfAnp ) ;

  double  _fastcall calcFncObj_MMP(TComp *pcmpSZv,TComp cmpZTrg, TComp cmpZAnp, double alfTrg, double alfAnp );

  bool  _fastcall find_ZTarg_and_ZAnt( TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnp, double alfTrg, double alfAnp );

  int  _fastcall find_AlfTarg_and_AlfAnt( TComp *pcmpSZv, TComp ZTarg, TComp ZAnt, double *alfTrg, double *alfAnp );

  void  _fastcall calc_F_and_dF_po_dAlf ( TComp *pcmpSZv, TComp ZTarg, TComp ZAnt
, double alfTrg, double alfAnp, double* arr_FGreek,double* arr_dFGreek );

  int  _fastcall solvNewtonMeth_Razn(const double valSigNoise,TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr );

  int  _fastcall calc_vectG_and_mtrxH_NewtonMeth_Razn (TComp *pcmpSZv, double alfTrg, double alfAnp
		,double* arr_FGreek, double*  arr_dFGreek , TComp*  pcmpZTarg,TComp* pcmpZAnp  );

  void  _fastcall calc_vectF_from_Alfa( TComp *pcmpSZv, double alfTrg, double alfAnp
  ,double* arr_FGreek, TComp*  pcmpZTarg,TComp* pcmpZAnp  )  ;

  void  _fastcall imitateMeasures( TComp *pcmpS,TComp *pcmpSZv, double alfTrg, double alfAnp
  ,TComp cmpZTarg, TComp cmpZAnt,double valNoiseSkz );

  void  _fastcall createDiagrGraphs(wchar_t *wchFoldName1 );

  void  _fastcall calcMtrxCorrel(const double valSigNoise,const double alfTrg, const double alfAnp
  , const TComp cmpZTarg,const TComp cmpZAnt,double* arr_dFGreekInv,double* arrMtrxCorr) ;

  void  _fastcall calcMtrxCorrel_Mistake(const double valSigNoise,const double alfTrg, const double alfAnp
  , const TComp cmpZTarg,const TComp cmpZAnt,double* arr_dFGreekInv,double* arrMtrxCorr);

  int  _fastcall solvNewtonMeth_RaznPerebor(const double valSigNoise
  , const double valXMax,const double valStepX, TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr );

  int  _fastcall solveMinMMP_Perebor_Na_Setke(const double valSigNoise
  , const int iQUANT_STEPS, TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr, double *parrObj , double *pvalMinObj );

 static double  _fastcall dblArrMin( double *arr,const int LENArr, int *pinum) ;



};


#endif
