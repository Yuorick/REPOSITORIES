//---------------------------------------------------------------------------

#ifndef FarH
#define FarH
#include "Faceta.h"
#include "Far_2D.h"
class TComp;
class TFaceta;
class TDiagrSet;
class TFar_2D;
class TFar
{
public:
 // ���������� ������� (�����) ��� ����� � ��������� ��������
 int m_N;
 // ���������� ����� ��������
 double  m_D;
 // ����� �����
 double	 mLambda ;
 // ������
 TFaceta mFaceta;
 //������ � ����������� ���� � ������� (��� �������)
 double *mparrDisp;
 // ������ � ����������� ���������� �� � �������
 double *mparrAmplFactDisp;



 TFar() ;
 ~TFar();

// ����������� �����������
 TFar(const TFar &R) ;
 TFar &operator=(const TFar  &R2) ;
 // ����� ������
 __fastcall TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta);
 // ����� ������
 __fastcall TFar(const int N);
  // ����� ������
 __fastcall TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double *parrDisp);
// ����� ������  3
 __fastcall TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double valFacetaDisp);
// ����� ������  4
 __fastcall TFar(TFar_2D Far_2D, const bool BVertical);
 // ����� ������  5
  __fastcall TFar(const TFar Far, const int NumRows);
  // ����� ������  6
   __fastcall TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double *parrDisp, double *parrAmplFactDisp);
double calcNoiseDisp();

double calcAmplFactDisp();

 double fncPeleng(TComp *pcmpSZv, TComp *pcarr_dPel_po_dSi);

 double calcDispPelengFnc(TComp *pcmpSZv);

 double fncEstimateRsmTetta(TComp *pcmpSZv,  double *pDisp ) ;

 int   fncEstimateMsd(TComp *pcmpSZv,  double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr );

void doArr(TComp cmpInp, double *arrRez);

double calcSmeshMoment( TComp *pcmpSZv,  TComp z1,TComp* cmparr_k1,TComp z2, TComp* cmparr_k2);

void   calcArr_k(TComp *pcmpSZv, TComp z, TComp *cmparr_k);

TComp calcFDeriv_po_z( TComp *pcmpSZv, TComp z);

double calcPhaseDisp( TComp *pcmpSZv, TComp z, TComp *cmparr_k);

void calcMtrxCorr_z1_z2(TComp *pcmpSZv, TComp z1, TComp z2,  double *arrMtrxCorr );

int   solvQuadrEqMeth(TComp *pcmpSZv, TComp * z1, TComp *z2 , double *arrMtrxCorr_fi );

void fncMeasureRSMProcessing(TComp *cmparrS, double *valXi2, double *valEstRSM, double *valRSMDisp  );

double  fncXi2_5P10_UM(double valNoiseSKZ, TComp *pcmpSZv, double valMuZv);

double  fncXi2_5P10_UM_GeneralCase(TComp *pcmpSZv, double valMuZv);

void fncCorrectDisp(TComp *cmparrCorectCoef);

void calcArrayStatCoeff(int quantMeas,TComp * pcmparrMeas, TComp *cmparrCorectCoef, TComp * cmparrCorectCoefStat) ;

TComp calcCorrectCoeff(TComp *pcmparrSZv, int quantMeas,int numArgMax, int numRow);

 void fncMakeArrK(TComp *cmparrCorectCoef,TComp *cmparrRo0
   , TComp *cmparrK) ;

void fncCorrectMeas(TComp *cmparrK,int quantMeas,TComp * pcmparrMeas
   ,TComp * pcmparrCorrectMeas);

void ImitateMeasureArray(double alfUMTrg,TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp, TComp *pcmpS, TComp *pcmpSZv);

void makeMonteCarloGraphs(wchar_t *wchFoldName, const int NIsp, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp);

void makeGraph_SKZ_from_Sig(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp,  double  valSigRangeMax) ;

 void makeGraph_SKZ_from_AngDiffer(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,TComp cmpKAntp) ;

 void makeGraph_SKZ_from_AntpPhAng(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
	   ,double alfUMAntp,TComp cmpKAntp);

 void fncMeasureWhiteProcessing(const double ValAngParam, TComp *pcmpSZv
 , double *valXi2, double *valEstWhite, double *valWhiteDisp  );

 TComp createSquareMeasure(const double ValAngPar, TComp *pcmpSZv);

 bool findRootWhiteEquation(double *valx0,const double ValAngParam
	 ,double valYZv, const double valToler, double *val_dx_po_dY );

double calcSig2Y_MethodWhite(const double ValAngParam, TComp *pcmpSZv);

int   fncEstimateMsd(TComp *pcmpSZv,  double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr , TComp *cmpZ1 , TComp *cmpZ2);

void makeGraph_SKZ_from_dBell(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp,  double  valCoefRange);

void ImitateMeasureArray(int quantTarg, double *arrUMtrg,TComp  *cmparrKTarg
			, TComp *pcmpS, TComp *pcmpSZv);

int   fncPolyn3(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ);

int   fncPolyn2(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ);

void makeMonteCarloGraphs_for_2_and_3_Targs(wchar_t *wchFoldName1, const int NIsp,const int QUantTarg, double *arrUMtrg, TComp *cmparrKTarg );

int   fncPolyn2(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ,  double *valXiSquare);

double calcXiSquare(int quantTarg, double *arrEstUMtrg,TComp  *cmparrEstKTarg , TComp *pcmpSZv) ;

int   fncPolyn3(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ, double *pvalXiSquare);

void makeGraphs_FuncRazr_from_AngTarg2_Hyp3(wchar_t *wchFoldName1, double *arrUMtrg, TComp *cmparrKTarg );

void calcTheoretical_LATarg_SKZ(bool *pbRazreshenie,double *pvalSigTarg
  , double alfUMTrg,  double alfUMAntp, TComp  cmpKTarg, TComp cmpKAntp);

//void Create4RowDiagrSet( TDiagrSet *pDiagrSet0);

//void makeMonteCarloGraphs1(wchar_t *wchFoldName1, const int NIsp, double alfUMTrg, TComp  cmpKTarg
	//		,double  alfUMAntp,TComp cmpKAntp);

double  calc_Guaranted_SKZ_RSM_For_TwoTargs ( double VAlUMTarg, double VAlUMAntp
	,  TComp  CMpKTarg,  TComp CMpKAntp);

static double calcTheorDisp_RSM(const double VAlNWaveCur, const double VAlNAppert
, const double VAlLamb, const double VAlTetta);

double  calcGuarantSystError_RSM_For_2Targs (const double VAlUMTarg, const double VAlUMAntp
	, const double VAlTargAmp, const double VAlAntpAmp, TComp *pcmpCircleCentre, double *pRadius, double *pMu);

double  calcMeanSystError_RSM_For_2Targs (const double VAlUMTarg, const double VAlUMAntp
	, const double VAlTargAmp, const double VAlAntpAmp_Ro, TComp *pcmpCentre
	, double *pvalRadius, double *pvalDeltaFi);

double  calc_Mean_SKZ_RSM_For_TwoTargs ( double valUMTarg, double valUMAntp
	,  TComp  cmpKTarg,  TComp cmpKAntp) ;

double findDiagrWidthApprox();

double  fncFFarApprox (const double valTetta);

	};

double max_d(double a, double b);
double calcDisp_RSM(const double VAlNWaveCur, const double VAlNAppert
, const double VAlLamb, const double VAlTetta);
#endif
