//---------------------------------------------------------------------------

#ifndef Article_Rad_El_2016H
#define Article_Rad_El_2016H

class TFar;
class TComp;
void makeMonteCarloGraphs(TFar Far, wchar_t *wchFoldName1, const int NIsp, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp);

int   solveMinMMP_Perebor_Na_Setke(const TFar Far
  , const int iQUANT_STEPS, TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr, double *parrObj , double *pvalMinObj );

double   dblArrMin( double *arr,const int LENArr, int *pinum);

void fncCreateMtrxV(const double valMu0, const double valMu1, TComp *cmparrV);



void fncCreateVectV(const double valMu,  TComp *cmpVectV);
#endif
