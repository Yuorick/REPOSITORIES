//---------------------------------------------------------------------------

#ifndef FiltAlfaH
#define FiltAlfaH
#include "LinFilt.h"
#include "Measure.h"


class TMeasure;
class TLinFilt;
class TComp;
class TFiltAlfa:public TLinFilt
{
public:

~TFiltAlfa() ;

TFiltAlfa() ;

TFiltAlfa(const double ValT,   double *parrEstX, double *parrMtrxK)  ;

virtual	void createMtrxA(const TMeasure Measure, const double Val_h, double *pMtrxA) ;

virtual	void createMtrxCorrU(const TMeasure Measure,const double Val_h,const double Val_Tau2, double *pMtrxCorrU) ;

virtual	void createMtrxC(const TMeasure Meas,const double Val_h, double *pMtrxC) ;

virtual	int fncStep( TMeasure Measure, const double Val_Tau2);

void fncAnalyticSolution(const double Val_h,const double Val_Tau2,const double ValDispKsi
  , double *arrK,double * arrP, TComp *pcmparrRoots
	 , double *arrT);

};


#endif
