//---------------------------------------------------------------------------

#ifndef FiltKnH
#define FiltKnH
#include "LinFilt.h"
class TLinFilt;
class TMeasure;
class TComp;
class TFiltKn:public TLinFilt
{
public:


	TFiltKn();

	TFiltKn(const double ValT,   double *parrEstX, double *parrMtrxK);

 virtual	void createMtrxA(const TMeasure Measure, const double Val_h, double *pMtrxA) ;

virtual void createMtrxCorrU(const TMeasure Measure,const double Val_h,const double Val_Tau2, double *pMtrxCorrU) ;

virtual void createMtrxC(const TMeasure Meas,const double Val_h, double *pMtrxC) ;

void StabSolutionKalm(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT);

static double  SolvCharactEqKlm(const double c);

virtual void fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT);

static double  SolvCharactEqKlm_Real(const double mu) ;

};
#endif
