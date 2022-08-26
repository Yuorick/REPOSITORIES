//---------------------------------------------------------------------------

#ifndef CalcCorMatrxH
#define CalcCorMatrxH
void calcMatrP1(const double a, double*matrP);
void calcMatrP2(const double a, double*matrP);
void calcMatrP3(const double a, double*matrP);
void calcMatrP2_AlogWatchArrow(const double a, double*matrP);

void calcMatr_ASK_v_KGSK(double *arrMu,double *matrPereh_ASK_V_KGSK) ;
void calcMatr_ASK_v_PSK(double *arrMu,double *matrPereh_ASK_V_PSK) ;
void calcMatr_PSK_v_KGSK(double *arrMu,double *matrPereh_PSK_V_KGSK) ;
void calcMatr_PSK_v_KGSK_Eilers(double *arrMu,double *matrPereh_PSK_V_KGSK);


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

void calcMatr_dP1_otMinus_a_po_da(const double a, double*matrP);

void calcMatr_d2P1otMinus_a_po_da2(const double a, double*matrP);

void calcMatr_dP3_po_da(const double a, double*matrP);

void calcMatr_d2P1_po_da2(const double a, double*matrP);

void calcMatr_d2P3_po_da2(const double a, double*matrP);

void calcCorMatrx_PSK_KGSK(const double VAlQ, const double VAlPsi, const double VAlTet
     , const double VAlDispQ,const double VAlDispPsi	 , const double VAlDispTet
     , double  *arrVS_PSK,double  *arrVS_KGSK, double *arrMtxCorr);

void calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer);

void calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer);

void recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps);

void recalcSphericalCoord_INTO_Rectangular(const double valR,const  double valBet
                                           ,const  double valEps, double *arrS);

void calc_dSpherical_po_Rectangular(double *arrS, double *arr_dB_po_dS);

void calc_dSpherical_po_Rectangular_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC);

void calc_dSpherical_po_Rectangular_AnglesOnly(double *arrS, double *arr_dB_po_dS);

void calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC);

void calc_dM_ask_to_psk_po_Betta(const double valBet,const double valEps,double *arr_dM_po_dBet);

void calc_dM_ask_to_psk_po_Eps(const double valBet,const double valEps,double *arr_dM_po_dEps);

void calc_dRectangularl_po_dSpherical(const double *arrV, double *arr_dA_po_dV);

void calc_d2B_po_ds_po_dS0(double *arrS, double *arr_d2B_po_dS_po_dS0);

void calc_d2B_po_ds_po_dS1(double *arrS, double *arr_d2B_po_dS_po_dS1);

void calc_d2B_po_ds_po_dS2(double *arrS, double *arr_d2B_po_dS_po_dS2);

void calc_dMtrx3_ASPK_v_PSK_po_dAlf(double *arrMu,double *arr_dM_po_dAlf);

void calc_dMtrx3_ASPK_v_PSK_po_dEps(double *arrMu,double *arr_dM_po_dAlf);

void calc_dMtrx3_ASPK_v_PSK_po_dBet(double *arrMu,double *arr_dM_po_dAlf);

void calc_dM_psk_to_kgsk_po_dQ(double *arrMu,double *arr_dM_po_dBet);

void calc_dM_psk_to_kgsk_po_dPsi(double *arrMu,double *arr_dM_po_dEps);

void calc_dM_psk_to_kgsk_po_dTet(double *arrMu,double *arr_dM_po_dAlf);

void calcMtrx3_ASPK_v_PSK(double *arrMu,double *matrPereh_PSK_V_KGSK);
//---------------------------------------------------------------------------




void calcMatrP1(const long double  a, long double *matrP);
void calcMatrP2(const long double  a, long double *matrP);
void calcMatrP3(const long double  a, long double *matrP);

void calcMatr_ASK_v_KGSK(long double  *arrMu,long double  *matrPereh_ASK_V_KGSK) ;
void calcMatr_ASK_v_PSK(long double  *arrMu,long double  *matrPereh_ASK_V_PSK) ;
void calcMatr_PSK_v_KGSK(long double  *arrMu,long double  *matrPereh_PSK_V_KGSK) ;



void calc_dF_po_dQ_sq(const long double  valBet,const long double  valEps,long double  *matrRez) ;
void calc_dF_po_dPsi_sq(const long double  valBet,const long double  valEps,long double  *matrRez);
void calc_dF_po_dTet_sq(const long double  valBet,const long double  valEps,long double  *matrRez) ;
void calc_dF_po_dBet_sq(const long double  valBet,const long double  valEps,long double  *matrRez) ;
void calc_dF_po_dEps_sq(const long double  valBet,const long double  valEps,long double  *matrRez) ;

void createExtendMtrx(long double  *arrInp, long double  * arrOut);

void calcMatrJ1(const long double  valBet,const long double  valEps,long double  *arrJ1) ;
void calcMatrJ2(const long double  valBet,const long double  valEps,long double  *arrJ2);
void calcMatrJ3(const long double  valBet,const long double  valEps,long double  *arrJ3) ;
void calcMatrJ4(const long double  valBet,const long double  valEps,long double  *arrJ4);

void calcMatrJ5(const long double  valBet,const long double  valEps,long double  *arrJ5) ;

void calcMatrJ1W(const long double  valBet,const long double  valEps,long double  *arrJ1W) ;

void calcMatrJ2W(const long double  valBet,const long double  valEps,long double  *arrJ2W);

void calcMatrJ3W(const long double  valBet,const long double  valEps,long double  *arrJ3W);

void calcMatrJ4W(const long double  valBet,const long double  valEps,long double  *arrJ4W) ;

void calcMatrJ5W(const long double  valBet,const long double  valEps,long double  *arrJ5W);

void calcMatr_dP1_po_da(const long double  a, long double *matrP);

void calcMatr_dP2_po_da(const long double  a, long double *matrP);

void calcMatr_dP1_otMinus_a_po_da(const long double  a, long double *matrP);

void calcMatr_d2P1otMinus_a_po_da2(const long double  a, long double *matrP);

void calcMatr_dP3_po_da(const long double  a, long double *matrP);

void calcMatr_d2P1_po_da2(const long double  a, long double *matrP);

void calcMatr_d2P3_po_da2(const long double  a, long double *matrP);

void calcCorMatrx_PSK_KGSK(const long double  VAlQ, const long double  VAlPsi, const long double  VAlTet
     , const long double  VAlDispQ,const long double  VAlDispPsi	 , const long double  VAlDispTet
     , long double   *arrVS_PSK,long double   *arrVS_KGSK, long double  *arrMtxCorr);

void calcMatrxPer_from_DecartPrSK_To_SSK(long double  *arrV, long double  * parrMtrxPer);

void calcMatrxPer_from_SSK_To_DecartPrSK(long double  *arrV, long double  * parrMtrxPer);

void recalcCoord_INTO_Spherical(long double  *arrInp, long double  &valR, long double  &valBet, long double  &valEps);

void recalcSphericalCoord_INTO_Rectangular(const long double  valR,const  long double  valBet
                                           ,const  long double  valEps, long double  *arrS);

void calc_dSpherical_po_Rectangular(long double  *arrS, long double  *arr_dB_po_dS);

void calc_dSpherical_po_Rectangular_AnglesOnly(long double  *arrS, long double  *arr_dB_po_dS);

void calc_dM_ask_to_psk_po_Betta(const long double  valBet,const long double  valEps,long double  *arr_dM_po_dBet);

void calc_dM_ask_to_psk_po_Eps(const long double  valBet,const long double  valEps,long double  *arr_dM_po_dEps);

void calc_dRectangularl_po_dSpherical(const long double  *arrV, long double  *arr_dA_po_dV);

void calc_d2B_po_ds_po_dS0(long double  *arrS, long double  *arr_d2B_po_dS_po_dS0);

void calc_d2B_po_ds_po_dS1(long double  *arrS, long double  *arr_d2B_po_dS_po_dS1);

void calc_d2B_po_ds_po_dS2(long double  *arrS, long double  *arr_d2B_po_dS_po_dS2);

void calc_dM_psk_to_kgsk_po_dBetta(long double  *arrMu,long double  *arr_dM_po_dBet);

void calc_dM_psk_to_kgsk_po_dEps(long double  *arrMu,long double  *arr_dM_po_dEps);

void calc_dM_psk_to_kgsk_po_dAlf(long double  *arrMu,long double  *arr_dM_po_dAlf);

void calcMtrx3_ASPK_v_PSK(long double  *arrMu,long double  *matrPereh_PSK_V_KGSK);
#endif
