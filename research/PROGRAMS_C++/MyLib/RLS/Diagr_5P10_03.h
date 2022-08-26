//---------------------------------------------------------------------------

#ifndef Diagr_5P10_03H
#define Diagr_5P10_03H
class TComp;
double findPartDiagrWidth07_GeneralizedAng();
double findPartDiagrWidth07_Rad(const double VAlLamb,const double VAlApert);

void  createGraphsPartial_and_Sum_Diagrams(wchar_t *wchFoldName1, const double valLamb
	, const double  valAppert , const double  valAngSdvig );

double transformAngToGeneralizedAng (const double valLamb , const double  valAppert
, const double  valAng  ) ;

double transformGeneralizedAngToAng (const double valLamb , const double  valAppert
, const double  GeneralizedAng  ) ;

double findTet07_For_SumDiagr_5P10_03( const double  valAngSdvig);

double findFirstZero_For_SumDiagr_5P10_03( const double  valAngSdvig);

double findSecondZero_For_SumDiagr_5P10_03( const double  valGenAngSdvig) ;


double fncSumDiagr_5P10_03( const double  valAngSdvig, double tetGeneralized) ;

double fncDerivSumDiagr_5P10_03( const double  valAngSdvig, double tetGeneralized);

double findSumDiagrWidth07_GeneralizedAng(const double  valAngSdvig);

double findSumDiagrWidth07_Rad(const double  valAngSdvig, const double VAlLamb,const double VAlApert)  ;

double findCrossLevel_For_SumDiagr_5P10_03( const double  valAngSdvig) ;

void ImitateMeasureArrayPartialDiagrams_5P10_03( double valLamb,double valAppert , double valAngSdvig
		, double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp, const double VAlNoiseSkz
		,double mAmplFactSig, TComp *cmparrPartS, TComp *cmparrPartSZv);

void ImitateMeasureArraySumDiagrams_5P10_03(  TComp *cmparrPartS, TComp *cmparrPartSZv,  TComp *cmparrSumS, TComp *cmparrSumSZv);

int EstGenAngsThreeSumDiagr_5P10_03(wchar_t *wchFoldName1,  TComp *cmparrS, double valGenAngSdvig
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1) ;

double fncFGr_SumDiagr(double *arr_b, double valAngSdvig, double valmu);

int estimateMethThreeSumDiagr_5P10_03(wchar_t *wchFoldName1, double valNoiseSkz
  , double valmAmplFactSig, double valLambda,double valAppert , TComp *cmparrS
, int iNumRayTriple, double valGenAngSdvig
  , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 ) ;

bool calcMtrxCorrGenAngs_Meth3SumDiagr(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz,double valmAmplFactSig, double *arrMtrxCorrGenAngs )  ;

void createPtesentGraphs_for_SumDiagr( wchar_t *wchrPresntSumDiagrams, double valLambda
	   , double valAppert , double valAngSdvig, int iNumRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleX );

void createGraphSumDiagr_from_rad(wchar_t *FileName, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz, double  valSMeshenieVert, double valScaleY );

void createGraphFGreece_from_rad_For_3SumDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void createGraphsSKZ_from_AngDiffer_For_3SumDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp , double valNoiseSkz,double valmAmplFactSig);

void createGraphs_for_Compare_Diagrams(wchar_t *FileName, double valLambda
	,double  valAppert ,double  valAngSdvig ) ;

 int findRootsFgr_For_3SumDiagr_5P10_03(  double valGenAngSdvig, double *arr_b,  double *arrRoots);

double 	  findRootMethChord_For_Fgr_For_3SumDiagr_5P10_03(  double valGenAngSdvig
	 , double *arr_b, double  valX0, double valX1)  ;

void createGraphFGreece_For_3SumDiagr_from_GenAng(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 ) ;

bool calcMtrxCorrGenAngs_Meth3PartDiagr(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz, double valmAmplFactSig, double *arrMtrxCorrGenAngs );

int estimateMethThreePartDiagr_5P10_03(wchar_t *wchFoldName1, double valNoiseSkz
  , double valmAmplFactSig, double valLambda,double valAppert , TComp *cmparrS
, int iNumRayTriple, double valGenAngSdvig
  , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

int EstGenAngsThreePartDiagr_5P10_03(TComp *cmparrS, double valGenAngSdvig
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

bool calcVect_b(TComp *cmparrS, double *arr_b);

double fncFGr_PartDiagr(double *arr_b, double valGenAngSdvig, double valmu) ;

double 	  findRootMethChord_For_Fgr_For_3PartDiagr_5P10_03(  double valGenAngSdvig
	 , double *arr_b, double  valX0, double valX1);

 int findRootsFgr_For_3PartDiagr_5P10_03(  double valGenAngSdvig, double *arr_b,  double *arrRoots) ;

void createGraphs_for_Compare_PartDiagrams(wchar_t *Fold ,double  valGenAngSdvig  ) ;

void createGraphFGreece_For_3PartDiagr_from_GenAng(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valGenAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void createPtesentGraphs_for_PartDiagr( wchar_t *wchrPresntSumDiagrams, double valLambda
	   , double valAppert , double valAngSdvig, int iNumRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleY );

void createGraphPartDiagr_from_rad(wchar_t *FileName, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY );

void createGraphFGreece_from_rad_For_3PartDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void createGraphsSKZ_from_AngDiffer_For_3PartDiagr(wchar_t *Fold, double valLambda
	,double  valAppert ,double  valAngSdvig,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp , double valNoiseSkz, double valmAmplFactSig);

void createMtrxCorrForMeasuresForAnsambleOfPartDiagrs(double valGenTargEps,  double valGenAntpEps , double valGenAngSdvig
  ,TComp cmpKTarg ,TComp cmpKAntp , double valNoiseSkz, double valmAmplFactSig, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr );
























#endif
