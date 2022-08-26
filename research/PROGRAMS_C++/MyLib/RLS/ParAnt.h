//---------------------------------------------------------------------------

#ifndef ParAntH
#define ParAntH
#include <stdio.h>
#include <stdlib.h>
class TComp;
class TURPointXY;
class TURPolygon;
//ЧАСТОТА       A11       A12       A21       A22

const double constArrMtrxTransf[] = {
//ЧАСТОТА           A11             A12          A21           A22
9.4000000000, -0.3515524743, -0.3515524743, 0.5741652505, -0.5741652505,
9.4500000000, -0.3439380720, -0.3439380720, 0.5643911959, -0.5643911959,
9.5000000000, -0.3362811185, -0.3362811185, 0.5545952986, -0.5545952986,
9.5500000000, -0.3285821185, -0.3285821185, 0.5447773585, -0.5447773585,
9.6000000000, -0.3208415984, -0.3208415984, 0.5349349363, -0.5349349363,
9.6500000000, -0.3130601067, -0.3130601067, 0.5250631437, -0.5250631437,
9.7000000000, -0.3052382145, -0.3052382145, 0.5151641596, -0.5151641596,
9.7500000000, -0.2973765155, -0.2973765155, 0.5052305373, -0.5052305373
};


//ТАБЛИЦА ТОЧЕК СДВИГА НАЧАЛА КООРДИНАТ ПО ЧАСТОТАМ


const double constArrPointsSdvig[] = {
//ЧАСТОТА           X0              Y0
9.4000000000, -0.9222627816, -0.9222627816,
9.4500000000, -0.9537500811, -0.9537500811,
9.5000000000, -0.9868512460, -0.9868512460,
9.5500000000, -1.0216896229, -1.0216896229,
9.6000000000, -1.0584014120, -1.0584014120,
9.6500000000, -1.0971373843, -1.0971373843,
9.7000000000, -1.1380648828, -1.1380648828,
9.7500000000, -1.1813701616, -1.1813701616
} ;


//       ЭЛЛИПТИЧЕСКИЙ СЛУЧАЙ
const int LEnArrCoefPolinEllipticCase = 7;
const int QUantCasesEllipticType = 0 ;
//  ТАБЛИЦА КОЭФФИЦИЕНТОВ ПОЛИНОМОВ СТЕПЕНИ 6
 double constArrCoeffPolinEllipticCase[] = {
//ЧАСТОТА       Диап.арг.     Разреш,мрад      valKSlop            1           X^1           X^2          X^3            X^4           X^5            X^6
9.5500000000, 3.0572411341, -30.3632903479, -0.0776537071, -0.0003975492, 1.0088794325, -0.0199986270, -0.2992760815, 0.1789688375, -0.0441963169, 0.0041319748,
9.6000000000, 2.9914865847, -30.3417271104, -0.1381808151, -0.0003232309, 1.0070716172, -0.0124840149, -0.3132749606, 0.1893222200, -0.0475439604, 0.0045318288,
9.6500000000, 2.9462198231, -30.3203873262, -0.1788722381, -0.0002707524, 1.0056148782, -0.0061434212, -0.3242932610, 0.1974896329, -0.0502622756, 0.0048686182,
9.7000000000, 2.9087669553, -30.2992675398, -0.2114259167, -0.0002211934, 1.0042139567, 0.0000751816, -0.3347152559, 0.2052636071, -0.0529091384, 0.0052056755,
9.7500000000, 2.8755104918, -30.2783643666, -0.2391756828, -0.0001671033, 1.0027509241, 0.0066085616, -0.3454441790, 0.2133289980, -0.0557073459, 0.0055699232
};




//   ГИПЕРБОЛИЧЕСКИЙ СЛУЧАЙ
//ТАБЛИЦА ОБЛАСТЕЙ ОПРЕДЕЛЕНИЯ ДЛЯ ЛЕВОЙ ВЕТВИ

// Левая ветвь гиперболы. Формат таблицы:
// первый столбец - частота
// второй столбец - tmin -нижняя грантца области определения переменной полинома P(t)
// третий столбец - tmax -верхняя грантца области определения переменной полинома P(t)
// четвертый столбец - TettaMin -нижняя граница области значений полинома P(t)
// пятый столбец - TettaMax - верхняя граница области значений полинома P(t)


const double constArrAreaCoeffHypCaseLeft[] = {
//ЧАСТОТА     tmin       tmax        TettaMin    TettaMax
9.4000000000, 0.0000000000, -14.5062874471, 0.0000000000, 18.2118797579,
9.4500000000, 0.0000000000, -9.3569192538, 0.0000000000, 17.8270038917,
9.5000000000, 0.0000000000, -11.0036008516, 0.0000000000, 17.4406601541,
9.5500000000, 0.0000000000, -8.9013086547, 0.0000000000, 17.0638522006,
9.6000000000, 0.0000000000, -9.1840545996, 0.0000000000, 16.6855076253,
9.6500000000, 0.0000000000, -9.2970111828, 0.0000000000, 16.3110837192,
9.7000000000, 0.0000000000, -9.1374540566, 0.0000000000, 15.9405198534,
9.7500000000, 0.0000000000, -8.7723928106, 0.0000000000, 15.5737566426
};


//ТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ПОЛИНОМА НЕЧЕТНОЙ СТЕПЕНИ
//ДЛЯ ЛЕВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]

const int LEnArrCoefPolinHypCaseLeft = 7;
const int QUantCasesHyperbolicType = 8 ;
const double constArrCoeffPolinHypCaseLeft[] = {
//ЧАСТОТА        X^1            X^3             X^5           X^7           X^9           X^11          X^13
9.4000000000, -0.8472033388, 0.0735249635, -0.0062599658, 0.0003695303, -0.0000133332, 0.0000002612, -0.0000000021,
9.4500000000, -0.8326301419, 0.0718696952, -0.0060785505, 0.0003572399, -0.0000128394, 0.0000002501, -0.0000000020,
9.5000000000, -0.8183117756, 0.0705903198, -0.0059835915, 0.0003534515, -0.0000127879, 0.0000002510, -0.0000000020,
9.5500000000, -0.8033783080, 0.0690001051, -0.0058086504, 0.0003401121, -0.0000122053, 0.0000002383, -0.0000000019,
9.6000000000, -0.7886743381, 0.0677068011, -0.0057066111, 0.0003349782, -0.0000120816, 0.0000002376, -0.0000000019,
9.6500000000, -0.7729790825, 0.0659508589, -0.0055543759, 0.0003275695, -0.0000118824, 0.0000002346, -0.0000000019,
9.7000000000, -0.7580276359, 0.0644122779, -0.0054215873, 0.0003204357, -0.0000116429, 0.0000002300, -0.0000000019,
9.7500000000, -0.7435577007, 0.0631064018, -0.0052945476, 0.0003116507, -0.0000113081, 0.0000002239, -0.0000000018

};


//ТАБЛИЦА ОБЛАСТЕЙ ОПРЕДЕЛЕНИЯ ДЛЯ ПРАВОЙ ВЕТВИ

// Правая ветвь гиперболы. Формат таблицы:
// первый столбец - частота
// второй столбец - tmin -нижняя грантца области определения переменной полинома P(t)
// третий столбец - tmax -верхняя грантца области определения переменной полинома P(t)
// четвертый столбец - TettaMin -нижняя граница области значений полинома P(t)
// пятый столбец - TettaMax - верхняя граница области значений полинома P(t)


const double constArrAreaCoeffHypCaseRight[] = {
//ЧАСТОТА     tmin      tmax      TettaMin  TettaMax
9.4000000000, -5.7008141370, -0.8787808573, 18.2118797579, 36.9306698001,
9.4500000000, -5.6633915460, -0.8557457695, 17.8270038917, 36.8865660803,
9.5000000000, -5.6708054844, -0.8333343298, 17.4406601541, 36.8429266102,
9.5500000000, -5.6241711814, -0.8115185384, 17.0638522006, 36.7997440979,
9.6000000000, -5.6202804841, -0.7902718984, 16.6855076253, 36.7570114035,
9.6500000000, -5.6096878902, -0.7695692246, 16.3110837192, 36.7147215349,
9.7000000000, -5.5914344747, -0.7493864654, 15.9405198534, 36.6728676443,
9.7500000000, -5.5646742990, -0.7297005356, 15.5737566426, 36.6314430244
};


//ТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ОТРЕЗКА РЯДА ЛОРАНА
//ДЛЯ ПРАВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]
const int LEnArrLoranCoeff = 7;
const double constArrLoranCoeffPolinHypCaseRight[] = {
//ЧАСТОТА          1/(X^2)         1/X           1             X^0            X^1             X^2           X^3
9.4000000000, -0.9408551552, -4.0351990713, -0.4177925462, -0.4403854452, -0.0451017264, -0.0017899539, 0.0000047493,
9.4500000000, -0.9122413043, -3.9748155335, -0.4173475122, -0.4351714416, -0.0445591466, -0.0017543139, 0.0000062557,
9.5000000000, -0.8857867623, -3.9201192664, -0.4235113183, -0.4335512656, -0.0450195954, -0.0018571494, 0.0000002624,
9.5500000000, -0.8573047144, -3.8590164829, -0.4231560193, -0.4285486675, -0.0445501924, -0.0018308152, 0.0000013700,
9.6000000000, -0.8306120348, -3.8025549761, -0.4284532606, -0.4267031458, -0.0449810039, -0.0019307977, -0.0000044543,
9.6500000000, -0.8037549392, -3.7449764875, -0.4330968254, -0.4246368188, -0.0453703396, -0.0020259731, -0.0000100206,
9.7000000000, -0.7767736707, -3.6863034906, -0.4371155373, -0.4223631493, -0.0457225343, -0.0021169746, -0.0000153628,
9.7500000000, -0.7497168360, -3.6265869584, -0.4405741099, -0.4199162263, -0.0460479856, -0.0022053037, -0.0000205632
};
class TParAnt
{
public:
 // к-во парциальных диагшрамм в веере
 int mQuantDiagr;
 // апертура антенны
 double   mAppert;
 // длина волны
 double	 mLambda ;
 // угол сдвига парциальных диаграмм  радианы
double  mAngSdvig ;
 //массив с дисперсиями шума в парциальных диаграммах
 double mNoiseDisp;

 // разброс коэффиц усиления k в парциальной диаграмме
 // то есть, коэффиц усиления равен P = (1 + k)P0
double mAmplFactSig;




 TParAnt() ;


// конструктор копирования
 TParAnt(const TParAnt &R) ;
 TParAnt operator=(TParAnt  R2) ;
 // парам констр
 __fastcall TParAnt(const int quantDiagr,const double Appert,const double Lambda
   ,const double AngSdvig, const double NoiseDisp, const double AmplFactSig);

static double findPartDiagrWidth07_GeneralizedAng() ;

double   findPartDiagrWidth07_Rad();

double transformGeneralizedAngToAng (const double  GeneralizedAng  ) ;

double transformAngToGeneralizedAng (const double  valAng  ) ;

void   createGraphsPartial_and_Sum_Diagrams(wchar_t *wchFoldName1) ;

double findSecondZero_For_SumDiagr();

double   findFirstZero_For_SumDiagr();

double   findTet07_For_SumDiagr( );

double   fncSumDiagr( double tetGeneralized);

double   fncDerivSumDiagr( double tetGeneralized);


double   findSumDiagrWidth07_GeneralizedAng();

double   findSumDiagrWidth07_Rad();

double findCrossLevel_For_SumDiagr();

void   ImitateMeasureArrayPartialDiagrams( double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp
		, TComp *cmparrPartS, TComp *cmparrPartSZv)  ;

void   ImitateMeasureArraySumDiagrams(  TComp *cmparrPartS, TComp *cmparrPartSZv,  TComp *cmparrSumS, TComp *cmparrSumSZv);

static bool calcVect_b(TComp *cmparrS, double *arr_b);

double 	    findRootMethChord_For_Fgr_For_3SumDiagr(double *arr_b, double  valX0, double valX1);

 double   fncFGr_SumDiagr(double *arr_b, double valmu);

 int   findRootsFgr_For_3SumDiagr(double *arr_b,  double *arrRoots);

 int   EstGenAngsThreeSumDiagr(TComp *cmparrS
  ,  double  *pvalGenTargEps,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

int   estimateMethThreeSumDiagr( TComp *cmparrS
, int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 ) ;

bool   calcMtrxCorrGenAngs_Meth3SumDiagr(const int NumAnsamble, double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp , double *arrMtrxCorrGenAngs );

void createMtrxCorrForMeasuresForAnsambleOfSumDiagrs(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorrSumDiagr ) ;

void createMtrxCorrForMeasuresForAnsambleOfPartDiagrs(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr )  ;

void   createPtesentGraphs_for_SumDiagr( wchar_t *wchrPresntSumDiagrams, int iNumRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleY );

void   createGraphSumDiagr_from_rad(wchar_t *FileName,double  valSmeshenieGoriz
  , double  valSMeshenieVert,  double valScaleY )  ;

void   createGraphFGreece_from_rad_For_3SumDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void   createGraphFGreece_For_3SumDiagr_from_GenAng(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void   createGraphsSKZ_from_AngDiffer_For_3SumDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp );

void createGraphs_for_Compare_Diagrams(wchar_t *Fold  );

int   estimateMethThreePartDiagr( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

int   EstGenAngsThreePartDiagr(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

 int   findRootsFgr_For_3PartDiagr(double *arr_b,  double *arrRoots);

double   fncFGr_PartDiagr(double *arr_b, double valmu);

double 	    findRootMethChord_For_Fgr_For_3PartDiagr(double *arr_b, double  valX0, double valX1)  ;

bool calcMtrxCorrGenAngs_Meth3PartDiagr(const int NumAnsamble, double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, double *arrMtrxCorrGenAngs );

void   createGraphFGreece_For_3PartDiagr_from_GenAng(wchar_t *Fold, double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 );

void   createPtesentGraphs_for_PartDiagr( wchar_t *wchrPresntPartDiagrams, int iNumPartRayTriple
		,double  alfUMTrg,TComp cmpKTarg,double  alfUMAntp, TComp cmpKAntp
		,double  valEstRadAngTarg,double  valEstRadAngAntp
		  ,TComp cmpKTargEst , TComp cmpKAntpEst, double val_b0, double val_b1, double valSMeshenieVert
		  , double valScaleY ) ;

void   createGraphPartDiagr_from_rad(wchar_t *FileName,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY );

void createGraphPartDiagr_from_Gen(wchar_t *FileName,double  valSmeshenieGoriz, double  valSMeshenieVert,  double valScaleY );

void   createGraphFGreece_from_rad_For_3PartDiagr(wchar_t *Fold,double  valSmeshenieGoriz
	, double  valSMeshenieVert, double val_b0, double val_b1 ) ;

void   createGraphsSKZ_from_AngDiffer_For_3PartDiagr(wchar_t *Fold ,double  valSmeshenieGoriz
	, double  valSMeshenieVert, TComp cmpKTarg, TComp cmpKAntp ) ;

void createGraphs_for_Compare_PartDiagrams(wchar_t *Fold );

double  calcRadAngSdvigaPartDiagr(const int numPartDiagr);

double  calGenAngSdvigaPartDiagr(const int numPartDiagr);

double  calcRadAngSdvigaSumDiagr(const int numSumDiagr);

double  calGenAngSdvigaSumDiagr(const int numSumDiagr);

void createGraphs_For_BearingFuncs( wchar_t *wchFoldName1)  ;

double  calcFirstCoefTalorBearFunc_For_SumDiagr();

double  calcFirstCoefTalorBearFunc_For_PartDiagr();

double  calcSecondDerivBearFunc_For_PartDiagr(const double tet);

double  calcThirdCoefTalorBearFunc_For_PartDiagr();

int estimateARSMethPartDiagr( TComp *cmparrPartSZv
		  ,const int NumPartRayPare,  double *pvalEstPartDiaRadAngTarg
		  ,  TComp *pcmpPartDiaKTarg , double *pvalPartDiaDisp, double *pPsi2);

int estimateARSMethSumDiagr(TComp *cmparrSumSZv ,const int NumSumRayPare
   ,double *pvalEstSumDiaRadAngTarg, TComp *pcmpSumDiaKTarg ,double *pvalSumDiaDisp, double *pPsi2) ;

void createGraphsRezImitARSM(wchar_t *wchrImitRezARSMFold, const double AlfRadUMTrg
	   , const TComp CmpKTarg, const int NumSumRayPare, const int NumPartRayPare, const int QuantIspit) ;

void createGraphsSKZ_ARSM(wchar_t *wchFoldName1, TComp cmpKTarg);

void  createGraphsForTabulation1(wchar_t *wchrTableFold,  double valTargAngRad
  ,double valAntpAngRad, TURPointXY  pntSdvig, double *arrMtxPer_From_b_To_Image_b);

void calcApproxPolinomCoeffArray_ForHyperbolicCase (wchar_t * wchrTableFold, const int NUmValuePoints
  , const int NPolinom,double *arrCoeff
  , double *pvalLeftTabDiapMin,double *pvalLeftTabDiapMax, double *pvalLeftRazreshenieMinRad, double *pvalLeftRazreshenieMaxRad
  , double *pvalRightTabDiapMin,double *pvalRightTabDiapMax, double *pvalRightRazreshenieMinRad, double *pvalRightRazreshenieMaxRad
  , double *pvalGenMaxError);

void CreateValueArraysHyperbolicType( const double VAlDiapMin, const double VAlDiapMax
   ,double *arrArgs,double *arrFunc, const int NUmPoints0,  int *piNumPointsReal );

static bool  CreateRegTab(double *arrArg, double *arrFunc, const int LEnArr
  ,double *arrRegularTabFunc , const int LEnRegTab, const double STepTab);

void calcApproxPolinomCoeffArray_ForEllipticCase (wchar_t * wchrTableFold, const int NUmValuePoints
  , const int NPolinom,double *arrCoeff , double *pvalTabDiap, double *pvalGenMaxError
	, double *pvalRazreshenie, double *valSlop);

//void calcApproxPolinomCoeffArray_ForEllipticCase_v1 (wchar_t * wchrTableFold, const int NUmValuePoints
 // , const int NPolinom,double *arrCoeff , double *pvalTabDiap, double *pvalGenMaxError, double *pvalRazreshenie);

bool  calcTransformationParams( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b);


void CreateValueArraysEllipsType( const double VAlDiap,double *arrArgs,double *arrFunc, const int NUmPoints0 ) ;

static bool  CreateRegularTab(double *arrArg, double *arrFunc, const int LEnArr
  ,double *arrRegularTabFunc , const int LEnRegTab, const double STepTab) ;

static int WriteTXTReportForPolinomCoefArray(const wchar_t*FileName,	double *arrOutMtrxTransf, double *arrOutVect00
 , double *arrOutPolinomCoefEllipticCase, const int numHyperbolic ,const int numElliptic
  , double *arrOutBoundsHyperbCaseLeft,double *arrOutPolinomCoefHyperbCaseLeft
  ,double *arrOutBoundsHyperbCaseRight,double *arrOutPolinomCoefHyperbCaseRight,const int NPolinom);

void transformArrayEnvelopePoints(const double VAlDiap
  , const int NUmPoints0, TURPointXY *arrPoints0, int *numPointsActual, TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b);

void calcEvelopePointsArray( const double VAlDiap, const int NUmPoints0, TURPointXY *arrPoints, int *numPointsActual )  ;

bool calcEnvelopePoint(const double VAlMu, TURPointXY &pntB);

TURPolygon calcEvelopePolyg( const double VAlDiap, const int NUmPoints0 );

TURPolygon transformEnvelopePolygon( const double VAlDiap, const int NUmPoints0
  , TURPolygon &plgEnvelopeRotated, TURPointXY *ppntSdvig, double *arrMtxPer_From_b_To_Image_b);

bool calc_b2( const double VAlMu, const double VAl_b1, double *pval_b2) ;

double 	solvSingularEquationMethChord( double  valX0, double valX1);

double fncFi(const double VAlMu0 );

bool  findEnvelopeSigularPoint (double *pvalMuSingular, double *pvalBSingular);

double fncCurvation();

void ShowTheoreticRotatedEnvelopeEllipseType(wchar_t *FileName );

 void ShowTheoreticRotatedEnvelopeHyperbolaType(wchar_t *FileName ) ;

void CreateValueArraysEllipsType__( const double VAlDiap,double *arrArgs,double *arrFunc, const int NUmPoints0 );

int tabulatedSolution_5P10_03( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

int tabulatedEstimationGenAngs_5P10_03(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

void  createGraphsForTabulation2(wchar_t *wchrTableFold, TComp *cmparrPartS
	 , TURPointXY  pntSdvig, double *arrMtxPer_From_b_To_Image_b);

bool findGenAng_UsingApproximation_ForHyperbolicCase(TURPointXY pntInp, double *pGenAng) ;

bool  calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp ,const int NumAnsamble, double *arrMtrxCorrGenAngs);

static bool findAlfaOverSwitch(wchar_t *wchrAlfaOverStitch, double *pvalAlfa, double *pvalb, double *pvalMu);

static double fncBSingular( double valAlf0,  double  valMuSingular) ;

static double fncFF( double valAlf0,  double  valMuSingular) ;

bool  calcTransformationParams_old( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b);

bool  calcTransformationParamsEllipticCase( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b);

bool  TParAnt::calcTransformationParamsHyperbolicCase( TURPointXY *pntSdvig, double *arrMtxPer_From_b_To_Image_b) ;

 bool calcAngOverSwitchForHyperbBranch (wchar_t * wchrTableFold, double *pvalGenOverSwitch);


};


void CalcStatParams_(const double valX, const int N
	, double &valAver,double  &valAverSquare);

double fncSKO_(const double valXAver ,const double valX2Aver);

void fncFlipDblArr(double *arrFunc , const int NUmPointsTab);

 double fncPolinom(double *arrCoeff, int lenarrCoeff, double valArg);

  double fncOddPolinom(double *arrCoeff, int lenarr, double valArg);

 void CalcPolinomCoefApprox(const int NPolinom, double *arrArg,double *arrF,const int  lenArr
   ,double * arrCoeff,  double *pMaxDelta);

 void CalcPolinomCoefForOddFunc(const int NPolinom, double *arrArg,double *arrF
   ,double * arrCoeff, int lenArr, double *pMaxDelta);

TURPointXY  calcCentrePoint( TURPolygon plgEnvelope);

int findMaxSemiAxe( TURPolygon plgEnvelopeCentred);

void fncPrintTab(FILE *fw, double *parr,  int nrows, int ncols )  ;

double fncPoinomApproximatedEstimationEllipticCase(double valArgTemp ,double valSlop
	 ,double *arrCoeff, int lenArrCoeff);

void fncCorrectArg(double *t);

int findNumRow(const double *arr, int  numRows, int  numCols, double valFreq);

void calc_dB0_po_dS(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS);

void calc_dB1_po_dS(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS);

double fncFTemp1(double valGenAng,double  VAlGenAngSdvig);

#endif
