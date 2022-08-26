#ifndef COORDSYSTTRSF_H
#define COORDSYSTTRSF_H

class QPolygonF;
class QCoordSystTrsf
{
public:
    QCoordSystTrsf();

    static void calcMatrP1(const double a, double*matrP);

    static void calcMatrP2(const double a, double*matrP);

    static void calcMatrP3(const double a, double*matrP);

    static void calcMatr_dP1_po_da(const double a, double*matrP);

    static void  calcMatr_dP2_po_da(const double a, double*matrP);

    static void  calcMatr_dP3_po_da(const double a, double*matrP);

    static void calcMatr_PSK_v_KGSK_LeftRot(double *arrMu,double *matrPereh_PSK_V_KGSK);

    static void calcMatr_PSK_v_KGSK_RightRot(double *arrMu,double *matrPereh_PSK_V_KGSK);

    static void   calc_dSpherical_po_Rectangular_AnglesOnly_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC);

    static void   calc_dSpherical_po_Rectangular_AnglesOnly(double *arrS, double *arr_dB_po_dS);

    static void   calc_dSpherical_po_Rectangular_MultCoeff(double valCoeff, double *arrS, double *arr_dB_po_dS_MultC);

    static void   calc_dSpherical_po_Rectangular(double *arrS, double *arr_dB_po_dS);

    static void   recalcSphericalCoord_INTO_Rectangular(const double valR,const  double valBet,const  double valEps
                                   , double *arrS);

    static void   recalcCoord_INTO_Spherical(double *arrInp, double &valR, double &valBet, double &valEps);

    static void calcMtrx3_ASPK_v_PSK(double *arrMu,double *matrPereh_PSK_V_KGSK);

    static void calc_dRectangularl_po_dSpherical(const double *arrV, double *arr_dA_po_dV);

    static void   calc_dMLeft_psk_to_kgsk_po_dQ(double *arrMu,double *arr_dM_po_dBet);

    static void   calc_dMLeft_psk_to_kgsk_po_dPsi(double *arrMu,double *arr_dM_po_dEps);

    static void   calc_dMLeft_psk_to_kgsk_po_dTet(double *arrMu,double *arr_dM_po_dAlf);

    static void   calcMatrxPer_from_DecartPrSK_To_SSK(double *arrV, double * parrMtrxPer);

    static void   calcMatrxPer_from_SSK_To_DecartPrSK(double *arrV, double * parrMtrxPer);

    static void   calcMatrxPer_from_SkSK_To_NormSK(const double psia ,
        const double nua ,const double gamma_a  ,double* parrMtrxPer);

    static void calcMatrxPer_from_SvSK_To_SkSK(const double alf ,
                       const double bet ,double* parrMtrxPer);

    static void calcMatrxPer_from_SvSK_To_NormSK_(double *arrVeloNSK,const double gamma_a
                          ,const double alf,const double bet ,double* parrMtrxPer);

    static void calcMatrxPer_from_NormSK_To_SvSK_(double *arrVeloNSK,const double gamma_a
                                ,const double alf,const double bet ,double* parrMtrxPer);

    static QPolygonF rotate(const QPolygonF plg0, const double ang);

    static void calc_dM_PSK_to_KGSK_po_dt (double *arrEilerCntrKP0, double *arr_dEilers_po_dt0
             , double *arr_dM_PSK_to_KGSK_po_dt);

    static void calc_dMtrx3_ASPK_v_PSK_po_dQ(double *arrMu,double *arr_dM_po_dQ);

    static void calc_dMtrx3_ASPK_v_PSK_po_dPsi(double *arrMu,double *arr_dM_po_dQ);

    static void calc_dMtrx3_ASPK_v_PSK_po_dTet(double *arrMu,double *arr_dM_po_dTet);
};

#endif // COORDSYSTTRSF_H
