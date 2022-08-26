//---------------------------------------------------------------------------

#ifndef RaznSumMethodH
#define RaznSumMethodH
#include "Faceta.h"
class TComp;
class TFaceta;

class TRaznSumMethod
{
public:


 // ���������� �����(��������)
 int m_NRows;
 // ���������� ����� ��������
 double  m_D;
 // ����� �����
 double	 mLambda ;
 // ������
 TFaceta mFaceta;
 // ������ �-�� ����� � ������
 int *mpiarrNum;
 // ������ ��������� ����� � ������� �� �������
 // mparrDisp[j] - ��������� ���� � ������� ������ � ������� j
 double *mparrDisp;

 ~TRaznSumMethod();
 TRaznSumMethod() ;


// ����������� �����������
 TRaznSumMethod(const TRaznSumMethod &R) ;
 TRaznSumMethod operator=(TRaznSumMethod  R2) ;
 TRaznSumMethod(const int N,const double D,const double Lambda
   ,TFaceta Faceta, int *piarrNum, double *parrDisp);
   /////////////////////////////////////////////////

 double fncPeleng(TComp *pcmpSZv, TComp *pcarr_dPel_po_dSi);

 double calcDispPelengFnc(TComp *pcmpSZv);

 void  imitateVectorIdealDiagrams(const double ValTetta, TComp cmpKTarg, TComp *pcmpS);

 void createIdealPelengFncGraph(wchar_t *wchFileName, const double ValRang );

 TRaznSumMethod correctFirstAndLastRows();

 void  imitateVectorMeasures(const double ValTetta, TComp cmpKTarg, TComp *pcmpS, TComp *pcmpSZv);

 double fncEstimateRsmTetta(TComp *pcmpSZv,  double *pDisp )  ;

 int calcQuantAM();

 double calcTotalDisp()  ;

 double fncEstimateRsmUM_with_ApproxDisp(TComp *pcmpSZv,  double *pDisp );

 double clcSystErr(const double ValTetta);

 void createSKZGraph(wchar_t *wchFileName, const double ValRang );
};

void fncBearingRSM(double valNoiseSKZ, TComp *pcmpSZv, double *pvalPelFunc, double *pvalDispPelFunc);
double  fncXi2_5P10_UM(double valNoiseSKZ, TComp *pcmpSZv, double valMuZv) ;

#endif
