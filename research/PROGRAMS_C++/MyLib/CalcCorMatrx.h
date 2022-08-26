//---------------------------------------------------------------------------

#ifndef CalcCorMatrxH
#define CalcCorMatrxH
void calcMatrP1(const double a, double*matrP);
void calcMatrP2(const double a, double*matrP);
void calcMatrP3(const double a, double*matrP);

void calcMatr_ASK_v_KGSK(double *arrMu,double *matrPereh_ASK_V_KGSK) ;
void calcMatr_ASK_v_PSK(double *arrMu,double *matrPereh_ASK_V_PSK) ;
void calcMatr_PSK_v_KGSK(double *arrMu,double *matrPereh_PSK_V_KGSK) ;



void calc_dF_po_dQ_sq(const double valBet,const double valEps,double *matrRez) ;
void calc_dF_po_dPsi_sq(const double valBet,const double valEps,double *matrRez);
void calc_dF_po_dTet_sq(const double valBet,const double valEps,double *matrRez) ;
void calc_dF_po_dBet_sq(const double valBet,const double valEps,double *matrRez) ;
void calc_dF_po_dEps_sq(const double valBet,const double valEps,double *matrRez) ;

void createExtendMtrx(double *arrInp, double * arrOut);

void calcMatrJ1(const double valBet,const double valEps,double *arrJ1) ;
void calcMatrJ2(const double valBet,const double valEps,double *arrJ2);
void calcMatrJ3(const double valBet,const double valEps,double *arrJ3) ;
void calcMatrJ4(const double valBet,const double valEps,double *arrJ4);

void calcMatrJ5(const double valBet,const double valEps,double *arrJ5) ;

void calcMatrJ1W(const double valBet,const double valEps,double *arrJ1W) ;

void calcMatrJ2W(const double valBet,const double valEps,double *arrJ2W);

void calcMatrJ3W(const double valBet,const double valEps,double *arrJ3W);

void calcMatrJ4W(const double valBet,const double valEps,double *arrJ4W) ;

void calcMatrJ5W(const double valBet,const double valEps,double *arrJ5W);

void calcMatr_dP1_po_da(const double a, double*matrP);

void calcMatr_dP2_po_da(const double a, double*matrP);

void calcMatr_dP3_po_da(const double a, double*matrP);

void calcCorMatrx_PSK_KGSK(const double VAlQ, const double VAlPsi, const double VAlTet
	 , const double VAlDispQ,const double VAlDispPsi	 , const double VAlDispTet
	 , double  *arrVS_PSK,double  *arrVS_KGSK, double *arrMtxCorr);

void calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer);

void calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer);







//---------------------------------------------------------------------------
#endif
