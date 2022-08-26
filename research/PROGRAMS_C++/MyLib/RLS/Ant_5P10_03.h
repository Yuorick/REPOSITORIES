//---------------------------------------------------------------------------

#ifndef Ant_5P10_03H
#define Ant_5P10_03H
// текст программы вычисления углов цели и антипода
// по измерениям 3 последовательно идущих диаграмм
// угол измеряется относительно центральной оси средней диаграммы

#include <stdio.h>
#include <stdlib.h>
class TComp;


//ЧАСТОТА       A11       A12       A21       A22

const double constArrMtrxTransf__[] = {
//ЧАСТОТА       A11           A12         A21         A22
9.4000000000, -0.0266463565, -0.0266463565, 0.1196834215, -0.1196834215,
9.4500000000, -0.0161557976, -0.0161557976, 0.0922329702, -0.0922329702,
9.5000000000, -0.0056609227, -0.0056609227, 0.0539057230, -0.0539057230,
9.5500000000, -0.0048377238, -0.0048377238, 0.0485743569, -0.0485743569,
9.6000000000, -0.0153500284, -0.0153500284, 0.0852693322, -0.0852693322,
9.6500000000, -0.0258704508, -0.0258704508, 0.1093052736, -0.1093052736,
9.7000000000, -0.0363919298, -0.0363919298, 0.1280859182, -0.1280859182,
9.7500000000, -0.0469081430, -0.0469081430, 0.1437060128, -0.1437060128
};


//ТАБЛИЦА ТОЧЕК СДВИГА НАЧАЛА КООРДИНАТ ПО ЧАСТОТАМ


const double constArrPointsSdvig__[] = {
//ЧАСТОТА         X0            Y0
9.4000000000, -18.2642914374, -18.2642914374,
9.4500000000, -30.4486422282, -30.4486422282,
9.5000000000, -87.8248233460, -87.8248233460,
9.5500000000, 103.8319842140, 103.8319842140,
9.6000000000, 33.0375551896, 33.0375551896,
9.6500000000, 19.7820731032, 19.7820731032,
9.7000000000, 14.1866487375, 14.1866487375,
9.7500000000, 11.0998075514, 11.0998075514
} ;


//       ЭЛЛИПТИЧЕСКИЙ СЛУЧАЙ
const int LEnArrCoefPolinEllipticCase__ = 7;
const int QUantCasesEllipticType__ = 5 ;
//  ТАБЛИЦА КОЭФФИЦИЕНТОВ ПОЛИНОМОВ СТЕПЕНИ 6
 double constArrCoeffPolinEllipticCase__[] = {
//ЧАСТОТА       Диап.арг.     Разреш,мрад      valKSlop            1           X^1          X^2            X^3              X^4            X^5         X^6
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


const double constArrAreaCoeffHypCaseLeft__[] = {
//ЧАСТОТА         tmin           tmax        TettaMin    TettaMax
9.4000000000, 0.0000000000, -3.6258026077, 0.0000000000, 2.9533718143,
9.4500000000, 0.0000000000, -4.0980082263, 0.0000000000, 2.3198014610,
9.5000000000, 0.0000000000, -5.1573813266, 0.0000000000, 1.3875592795
};


//ТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ПОЛИНОМА НЕЧЕТНОЙ СТЕПЕНИ
//ДЛЯ ЛЕВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]

const int LEnArrCoefPolinHypCaseLeft__ = 7;
const int QUantCasesHyperbolicType__ = 3 ;
const double constArrCoeffPolinHypCaseLeft__[] = {
//ЧАСТОТА        X^1             X^3             X^5            X^7          X^9         X^11            X^13
9.4000000000, -0.1841612512, 0.0152180805, -0.0014337667, 0.0001217256, -0.0000080073, 0.0000003138, -0.0000000048,
9.4500000000, -0.1425853561, 0.0118101573, -0.0010923232, 0.0000842397, -0.0000046397, 0.0000001514, -0.0000000021,
9.5000000000, -0.0836714166, 0.0067283213, -0.0005610416, 0.0000350588, -0.0000014197, 0.0000000322, -0.0000000003
};


//ТАБЛИЦА ОБЛАСТЕЙ ОПРЕДЕЛЕНИЯ ДЛЯ ПРАВОЙ ВЕТВИ

// Правая ветвь гиперболы. Формат таблицы:
// первый столбец - частота
// второй столбец - tmin -нижняя грантца области определения переменной полинома P(t)
// третий столбец - tmax -верхняя грантца области определения переменной полинома P(t)
// четвертый столбец - TettaMin -нижняя граница области значений полинома P(t)
// пятый столбец - TettaMax - верхняя граница области значений полинома P(t)


const double constArrAreaCoeffHypCaseRight__[] = {
//ЧАСТОТА       tmin               tmax      TettaMin      TettaMax
9.4000000000, -3.2111389022, -0.1828540080, 2.9533718143, 32.5443143977,
9.4500000000, -2.7082491422, -0.1433506183, 2.3198014610, 32.5108690768,
9.5000000000, -2.0762821729, -0.0855412924, 1.3875592795, 32.4777758120
};


//ТАБЛИЦА КОЭФФИЦИЕНТОВ АППРОКСИМИРУЮЩЕГО ОТРЕЗКА РЯДА ЛОРАНА
//ДЛЯ ПРАВОЙ ВЕТВИ ГИПЕРБОЛЫ. АППРОКСИМАЦИЯ НА ИНТЕРВАЛЕ [tmin;tmax]
const int LEnArrLoranCoeff__ = 7;
const double constArrLoranCoeffPolinHypCaseRight__[] = {
//ЧАСТОТА          1/(X^2)        1/X           1               X^0         X^1              X^2           X^3
9.4000000000, -0.0131285553, -0.7876363561, -0.0856156550, -0.1333386588, -0.0317733028, -0.0058913084, -0.0005164000,
9.4500000000, -0.0049324120, -0.5935287741, -0.0473534282, -0.0950083134, -0.0249526372, -0.0056887408, -0.0006087815,
9.5000000000, -0.0005636879, -0.3390883018, -0.0100444871, -0.0413027052, -0.0092358562, -0.0026860865, -0.0003961098
};
class TAnt_5P10_03
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




 TAnt_5P10_03() ;


// конструктор копирования
 TAnt_5P10_03(const TAnt_5P10_03 &R) ;
 TAnt_5P10_03 operator=(TAnt_5P10_03  R2) ;
 // парам констр
 TAnt_5P10_03(const int quantDiagr,const double Appert,const double Lambda
   ,const double AngSdvig, const double NoiseDisp, const double AmplFactSig);



   // ФУНКЦИИ ЧЛЕНЫ
double transformGeneralizedAngToAng__ (const double  GeneralizedAng  ) ;

double transformAngToGeneralizedAng__ (const double  valAng  ) ;


static bool calcVect_b__(TComp *cmparrS, double *arr_b);


void createMtrxCorrForMeasuresForAnsambleOfPartDiagrs__(double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp, int iNumAnsamble, int lenAnsamble
  , double *arrMtrxCorr )  ;


int   estimateMethThreePartDiagr__( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

int   EstGenAngsThreePartDiagr__(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps, double *pval_b0, double *pval_b1);

 int   findRootsFgr_For_3PartDiagr__(double *arr_b,  double *arrRoots);

double   fncFGr_PartDiagr__(double *arr_b, double valmu);

double 	    findRootMethChord_For_Fgr_For_3PartDiagr__(double *arr_b, double  valX0, double valX1)  ;

double  calcRadAngSdvigaPartDiagr__(const int numPartDiagr);

double  calGenAngSdvigaPartDiagr__(const int numPartDiagr);


int tabulatedSolution_5P10_03__( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr);

int tabulatedEstimationGenAngs_5P10_03__(TComp *cmparrS ,  double  *pvalGenTargEps
   ,  double  *pvalGenAntpEps);


bool findGenAng_UsingApproximation_ForHyperbolicCase__(double pntInpX, double pntInpY,  double *pGenAng) ;

double  calcRadAngSdvigaSumDiagr__(const int numSumDiagr);

double  calGenAngSdvigaSumDiagr__(const int numSumDiagr);

bool calcMtrxCorrGenAngs_AlternativeMeth3PartDiagr__(TComp *cmparrS,  double valGenTargEps,  double valGenAntpEps
  ,TComp cmpKTarg ,TComp cmpKAntp ,const int NumAnsamble, double *arrMtrxCorrGenAngs );

};



 double fncPolinom__(double *arrCoeff, int lenarrCoeff, double valArg);

  double fncOddPolinom__(double *arrCoeff, int lenarr, double valArg);


double fncPoinomApproximatedEstimationEllipticCase__(double valArgTemp ,double valSlop
	 ,double *arrCoeff, int lenArrCoeff);

void fncCorrectArg__(double *t);

int findNumRow__(const double *arr, int  numRows, int  numCols, double valFreq);

int findTangencyHyperbolaPoints(double  pntInp0X,double  pntInp0Y
 ,double &pntOut0X,double &pntOut0Y,  double &pntOut1X,double &pntOut1Y);

bool   InverseMtrx2__(double *arrA, double *arrOut);

void MtrxMultMatrx__(double *parrA,int nRowsA, int nColsA, double * parrB,int nColsB, double *parrRez);

void MtrxMultMatrxTransp__(double *parrA,int nRowsA, int nColsA, double * parrB,int nRowsB, double *parrRez) ;

void calc_dB1_po_dS__(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS);

void calc_dB0_po_dS__(TComp *cmparrS, double valGenAngSdvig, double val_b0, double val_b1, double *arrdB_po_dS) ;

int SolvEq2__(const double a,const double b,const double c,TComp &x1,TComp &x2);





#endif
