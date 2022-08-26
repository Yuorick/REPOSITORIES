//---------------------------------------------------------------------------

#ifndef DiagramsH
#define DiagramsH
class TComp;



int find_Targ_and_Antp_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 TComp *pcmpZ1, TComp *pcmpZ2, double *palfTrg, double *palfAnp ) ;

void calcF_and_dF_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF ) ;

 void  calcPartial_F_and_dF_Newton6 (double valWidthDgr,TComp cmpSZv, double valAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF );

 double calcFncObj_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
, TComp cmpZ1, TComp cmpZ2, double alfTrg, double alfAnp ) ;

double calcFncPartialObj_Newton6(double valWidthDgr,TComp cmpSZv, double valAlf,
 double(*pfncF)(double valWidthDgr, double tet), double *arrX );

int find_Targ_and_Antp_GroupRelax6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 TComp *pcmpZ1, TComp *pcmpZ2, double *palfTrg, double *palfAnp )  ;

int Relax_Angles_GroupRelax6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrZ, double *arrMlradAngles );

bool findOtimal_Z1_and_Z2(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet), TComp *pcmpZ, double *arrMilradAngles );

double SIGN__(double a);

double fncDiagrPar(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb, const double valTetScan,const double tet);

double fnc_F3(const double val_gam,const  double val_mu, const double tet);
double fnc_dF3_po_tet(const double val_gam,const  double val_mu, const double tet);
double fnc_F4(const double v,const  double u, const double v1,const  double u1);
double fnc_dDiagrPar_po_dTet(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb, const double valTetScan,const double tet);

double findDiagrWidth(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb);

void createDiagrGraphs_v1(const double val_d,const  double val_D, const int ival_n
  , const  int ival_N,const  double valLamb,const double valTetScan,wchar_t *wchFoldName1 )  ;


 void calcF_and_dF_AnglesRelax_Newton6(const int quantDiagr,double valWidthDgr,TComp *pcmpSZv, double *pAlf,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF );

 void calcPartial_F_and_dF_AnglesRelax_Newton6(double valWidthDgr,TComp cmpSZv, double valAlf0,
 double(*pfncF)(double valWidthDgr, double tet)
 , double(*pfncdF_po_dTet)(double valWidthDgr, double tet), double(*pfncd2F_po_dTet2)(double valWidthDgr, double tet),
 double *arrX, double * arr_F, double * arr_dF );



#endif
