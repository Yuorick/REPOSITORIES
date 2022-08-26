//---------------------------------------------------------------------------

#ifndef QuadrEqMethH
#define QuadrEqMethH
class TFaceta;
int   solvQuadrEqMeth(TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , TComp *z1,TComp *z2, double *arrMtrxCorr );




//---------------------------------------------------------------
void calcMtrxCorr_z1_z2(TComp *pcmpSZv, TComp z1, TComp z2,  double *arrMtrxCorr );
double calcPhaseDisp( TComp *pcmpSZv, TComp z, TComp *cmparr_k)  ;
void   calcArr_k(TComp *pcmpSZv, TComp z, TComp *cmparr_k);
double calcSmeshMoment( TComp *pcmpSZv,  TComp z1,TComp* cmparr_k1,TComp z2, TComp* cmparr_k2) ;
void doArr(TComp cmpInp, double *arrRez);
TComp calcFDeriv_po_z( TComp *pcmpSZv, TComp z);
int   solvQuadrEqMeth(TComp *pcmpSZv, TComp * z1, TComp *z2 , double *arrMtrxCorr_fi );
int   fncEstimateMsd(double valAMDist, TFaceta Faceta, TComp *pcmpSZv, double * arrDisp, double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr );


#endif
