//---------------------------------------------------------------------------

#ifndef ParabolicAnt_2DH
#define ParabolicAnt_2DH
class TSincDgr;
class TMeasStand;
class TParabolicAnt_2D
{
public:
 // к-во парциальных диагшрамм в веере
 int mQuantDiagr;
 // ширина диаграммы по уровню 0,5

 // длина волны
 double	 mLambda ;
 // абсолютгого угол сдвига парциальных диаграмм  радианы
double  *mparrAngAbsSdvig ;

 // массив диаграмм в веере
 TSincDgr *mparrDgr;



 TParabolicAnt_2D() ;
 __fastcall ~TParabolicAnt_2D() ;

// конструктор копирования
 TParabolicAnt_2D(const TParabolicAnt_2D &R) ;
 TParabolicAnt_2D operator=(TParabolicAnt_2D  R2) ;
 // парам констр
  __fastcall TParabolicAnt_2D(const int quantDiagr,const double Lambda
   ,double  *parrAngSdvig, TSincDgr *parrDgr);

_fastcall TParabolicAnt_2D(const int quantDiagr,const double Lambda
   , double *arrAngSdvig,  double *arrRadTet05, const double AmplFactSig,  const double NoiseSkz );

void ImitateMeasureArray( double valalfUMTrg, TComp cmpKTarg, double valalfUMAntp,TComp cmpKAntp
		, TComp *cmparrPartS, TComp *cmparrPartSZv);



int estimateMethThreeDiagr( TComp *cmparrS
  , int iNumRayTriple , double *valEstAngTarg, double *valEstAngAntp, TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr
  , double *pval_b0, double *pval_b1 );

 void processMeasureArrFromStandKaurovFile(wchar_t *wchFoldReport, int iNumTriple
	   ,TMeasStand* parrMeas, int iNumMeasures ) ;

void processMeasureArrFromStandKaurovFile_(wchar_t *wchFoldReport, int iNumTriple
	   ,TMeasStand* parrMeas, int iNumMeasures );

void processARSM_ForMeasureArrFromStandKaurovFile(wchar_t *wchFoldName1, int iNumPare
		,TMeasStand* parrMeas, int iNumMeasures, double valTargAngApriori );


};
bool isMeasure(TComp *pcmparrMeas) ;

   #endif
