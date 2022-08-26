#ifndef ADJUSTMENT_H
#define ADJUSTMENT_H

#define NUM_GADG_PARAMS 6
#define LENX_V1 9

class QAdjustment
{
public:
    QAdjustment();

    static void calc_VIU(double * arrTargS_KGSK , double * arrSins,
         double *arrX,  double &valR, double &valBet, double &valEps);

    static void calc_VIU_and_dVIU_po_dX(double * arrTargS_KGSK ,double * arrSins
                    ,double *arrX, double &valR, double &valBet, double &valEps, double *arr_dVIU_po_dX);

    static void recalc_VIU_from_RLK_to_IU(double * arrVZv_RLK
                                          ,double *arrX_RLK,double *arrX_IU, double &valBet, double &valEps);

    static void recalc_VIU_from_RLK_to_IU_and_dVIU_po_dX(double * arrVZv_RLK
     ,double *arrX_RLK,double *arrX_IU, double &valBet, double &valEps, double *arr_dVIU_po_dX);

    static void calc_dVIU_po_dXj(double *arrS_KGSK ,double * arrSins,  double *arrXGdg
                                              ,const int NUmj, double *arr_dVIU_po_dXj);

    static void calc_dVIU_po_dXj_V1(double *arrVZv_RLK,double *arrX_RLK
                                                 ,double *arrX_IU, const int NUmj, double *arr_dVIU_po_dXj);
};

#endif // ADJUSTMENT_H
