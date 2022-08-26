// ---------------------------------------------------------------------------

#ifndef HomingSituationH
#define HomingSituationH
#include "MainForm.h"
#include "Bomb.h"
#include "ShipTarg.h"
#include "BombTraj.h"
#include "Filt.h"
#include "CntlFuncPar.h"
#include "Glonass.h"
class  TBombTraj;
class  TShipTarg;
class  TCntlFuncPar;
class  TGlonass ;
class TURPointXY;
class TURPolygon;

///

class THomingSituation {
public:
	// корабль цель
	TShipTarg mShipTarg;
	// траектория АУБ
	TBombTraj mBombTraj;
	TCntlFuncPar mCntlFuncPar ;
	TFilt mFilt ;
	TGlonass mGlonass ;
	long double mTCur;


	// путь у папке с отчетом
	wchar_t *mpwcharrFoldReport;

	~THomingSituation();
	// конструктор по умолчанию
	THomingSituation();
	// конструктор копирования
	THomingSituation(const THomingSituation &R);

	// оператор присваивания
	THomingSituation operator = (THomingSituation R2);

	// парам конструктор 1
	THomingSituation(TShipTarg ShipTarg //
		, TBombTraj BombTraj //
		, TCntlFuncPar CntlFuncPar
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		);
   // парам конструктор 2
  THomingSituation(TShipTarg ShipTarg //
		, TBombTraj BombTraj //
		, TCntlFuncPar CntlFuncPar
		,TGlonass Glonass
		,long double TCur
		, wchar_t *pwcharrFoldReport // путь у папке с отчетом
		);

	long double fncAngVisir();
	long double fncDist();

	void   ImitateZamer(long double &valRZv,long double &valFiVisZv,long double &valTZv
	 , long double &valSig_BMO_R, long double &valSig_MMO_R,
	   long double &valSig_BMO_Fi, long double &valSig_MMO_Fi);

 
	void fncMoveBomb_TO_ZeroAlt_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz );

	long double fncCalcU(long double TCur);

	long double fncCalcDeltaAngVisir() ;

	bool fncSolveOCP_for_MaxDist_MethNewton(TCntlFuncPar &CntlFuncPar);

	bool  InitialApprox(long double* arrX);

	void calc_Fgr( long double* arrZ, long double* arrFgr);

	void calc_Fgr_and_H( long double* arrZ, long double* arrFgr, long double* mtrxH);

	void fncCalc_GradX1_and_GradX2_po_dTay__( long double valTEnd,int &quantTPerecl
	,long double* GradX1_po_dTay,long double* GradX2_po_dTay  );

	void fncCalc_GradX1_and_GradX2_po_dTay( long double valTEnd,int &quantTPerecl
	,long double* GradX1_po_dTay,long double* GradX2_po_dTay  )  ;


	void fncMove_to_EndTime_and_CollectMtrxs( long double* arrMtrxL, long double* arrMtrxB, int &quantTPerecl);

	void fncEilerStep_VS_and_MtrxL_and_MtrxB( const long double valStepInt
		 ,long double* mtrxL,long double* mtrxB ) ;

	bool IsSolutionTrue(long double* arrX, long double* arrXt);

  	THomingSituation fncFindOptimalControl_for_Dist_MethodPerebora_1Point(TCntlFuncPar &CntlFuncPar
	   ,const long double valShagT, long double& valDGorizOpt) ;

	static int fncSign (long double a);

	void fncMoveClas_TO_BlastPoint_AND_ShowGraphs( wchar_t *wcharrPath, long double &valDHoriz );
	TBombTraj fnc_Check_Possibility_Takeover_Perebor_2Point(const long double valShagT, long double& valDGorizOpt);

	bool fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_1Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn);

	bool fncMoveBomb_TO_TakeoverPoint_AND_ShowGraphs(wchar_t *wcharrPath
	  ,const long double valDeltaAimPoint,const long double valFunnelAng) ;

	bool fnc_FindPossibleAngle_For_Takeover_Perebor(const int iQuantPointsOverswitch
   , const long double valDeltaAimPoint ,const long double valFunnelAng,const long double valShagT
   , THomingSituation &HomingSituationReturn) ;

	void fncMoveClass_TO_FixedTime_AND_ShowGraphs( wchar_t *wcharrPath, const long double valFixedTime );

	bool fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_2Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn);

	 bool fnc_Check_Possibility_Takeover_For_FixedTraj_Perebor_3Point(const long double valDeltaAimPoint
	 ,const long double valFunnelAng,const long double valShagT, THomingSituation &HomingSituationReturn) ;

	  bool fnc_FindMaxDist_For_Takeover_Perebor(const int iQuantPointsOverswitch,const long double valShagT
   , THomingSituation &HomingSituationReturn, long double &valMaxDistTakeover) ;

	bool fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_1Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover);

	bool fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_2Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover);

	bool fnc_FindMaxDist_For_Takeover_For_FixedTetta0_Perebor_3Point(const long double valShagT
   , THomingSituation &HomingSituationReturn0, long double &valMaxDistTakeover);

	bool findTakeoverDist_AND_ShowGraphs(wchar_t *wcharrPath ,long double &valDistTakeover)  ;

	bool IsGeom(long double &valDistTakeover);

	bool THomingSituation:: IsHoming();

	static void fncCalcProbability_AND_ShowGraphs( wchar_t *pwchStat, const int NR, TBomb Bomb,long double valTBegin
, long double valTet0, long double valV0, long double valY0, long double valX0, long double valStepInt, long double valKNav
, TShipTarg ShipTarg ,long double *arrCfW, TGlonass Glonass, long double &valProbab ) ;

  bool fncMoveHomeHead_TO_BlastPoint_AND_ShowGraphs_StatMeth( wchar_t *wcharrPath
   , TBomb BombModel, long double &valProb,  long double &valSig_AngVis ) ;

  long double fncCalcU__( long double TCur);

  long double fncCalcU1(long double valEstFiVisTochka, long double TCur);
  void ImitateZamerGlonass(long double  *arrBombVSZv);

  long double calcPsi1();

  long double calcPsi2();

  THomingSituation findOptAng_for_MaxDist_AndShowGraphs(wchar_t *pwchMaxDistOutFile, const long double valAngStep, const int quantSteps) ;

  static void createGraph_OptMaxDist_fom_Height(long double mV0,long double valStepInt, wchar_t *mpwchMaxDistOutFile
	, long double valHMin,long double valHStep,int quantStep, long double &valDist, long double &valTPerecl, long double &valAng);

  static bool createGraph_MaxDist_For_Takeover_from_Height(long double valV0,long double valStepInt, TBomb Bomb
	 ,wchar_t *pwchMaxDistZahvOutFile, int QuantOverSwitchPnts	, long double valHMin
	 ,long double valHStep,int quantSteps, long double &valHeight,int &iTypeOfTraj
	 , long double &valMaxDistTakeover,long double &valTetta0, long double &valTPerecl
	 , long double &valTZahv, long double &valVXahv, long double &valYZahv,long double & valXZahv
	 , long double &valTettaZahv);

 bool findTakeoverDist_AND_CreatePictures(wchar_t *wcharrPath, long double &valDistTakeover) ;

// void drawPictures(wchar_t * pwchShapeUAB, wchar_t *pwcharrPath, const double valRastigenie);
void drawPictures( wchar_t *pwcharrPath, const double valRastigenie);

 void drawJetPictures( wchar_t *wcharrPath);

 static void createTransformedPlgShapeFile(wchar_t *pwchShapeFileInp, wchar_t *pwchShapeFileOut, const double  valAng
	, const TURPointXY pntCentre, const TURPointXY pntSdvig,const double valRastigenie) ;

// void drawEMB_Blast_Pictures(wchar_t *pwchShapeJet, wchar_t *wcharrPath)  ;
 void drawEMB_Blast_Pictures( wchar_t *wcharrPath)  ;

 //void drawVessel(wchar_t *pwchShapeFileInp, wchar_t *pwchShapeFileOut) ;

// static void drawGlonass(wchar_t *pwchShapeImages, wchar_t *pwcharrPathOut, double valRastigenie, double valXSdvig, double valYSdvig ) ;

 bool fncCreateGraph_Prob_from_FiDefeat( wchar_t *wcharrPath );
 void fncMoveBomb_TO_ZeroAlt_AND_ShowPictures( wchar_t *wcharrPath, long double &valDHoriz );


 void createPicturesConusZahv_and_ConusHoming(wchar_t *pwchOutFile);
  void createPicturesConusZahv(wchar_t *pwchOutFile);
  double     fnc_b( double valTet)  ;
  void createPicturesConusHoming(wchar_t *pwchOutFile);
  double    fnc_tetta( double valV) ;
  double     fnc_a( double valTet) ;
  void     fnc_f_and_df( double valV   ,double valTet0  ,double &valF,double &valDF);
 void createGraph_f(wchar_t *wcharrPath, double valV) ;
 void  createEllipsGraph(wchar_t *wchFileName);
 double THomingSituation::calcSigAngVisExtrapol();
 void createGraphs_P_from_RDefeat_FiDefeat(double *parrDist,
			double *parrSigAlfVisir,double *arrSigR, const int lenArrSigAlf , wchar_t *wcharrPath00) ;

 double calcSigRExtrapol();
};


#endif
