//---------------------------------------------------------------------------

#ifndef Kirn_Vers1H
#define Kirn_Vers1H
class TComp;
int   solvNewtonMeth_KirnStage1(const double valSigNoise,TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr ) ;

int   calc_vectG_and_mtrxH_NewtonMeth_Razn (TComp *pcmpSZv, double alfTrg, double alfAnp
		,double* arr_FGreek, double*  arr_dFGreek , TComp*  pcmpZTarg,TComp* pcmpZAnp  );

void  calc_vectF_from_Alfa__( TComp *pcmpSZv, double alfTrg, double alfAnp
  ,double* arr_FGreek, TComp*  pcmpZTarg,TComp* pcmpZAnp  );

bool  find_ZTarg_and_ZAnt__( TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnp, double alfTrg, double alfAnp );

void  calc_ATA_andATS__(int iNumDiagr, TComp cmpSZv, double alfTrg, double alfAnp , double*arrAATCur, double*arrATSCur );

void  calc_F__ ( TComp *pcmpSZv, TComp ZTarg, TComp ZAnt
, double alfTrg, double alfAnp, double* arr_FGreek );

void  calcPartial_F__(int iNumDiagr, TComp cmpSZv,TComp ZTarg, TComp ZAnp
		,double alfTrg,double alfAnp ,double* arr_Part_F);

TComp fncF__ (int iNumDiagr, const double tet);

double fncDiagrSimple__(double valWidthDgr, double valTetRad);

double fncDiagrSinx_div_x__(double tet) ;

TComp dF_po_dTet__(int iNumDiagr, const double tet);

double fncDerivDiagrSimple__(double valWidthDgr, double valTetRad);

double fncDerivDiagrSinx_div_x__(double tet);


#endif
