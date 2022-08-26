#ifndef LBLSOLVER_H
#define LBLSOLVER_H
#include "PosSolver.h"


class QLblSolver : public QPosSolver
{
public:

    QLblSolver();

    QLblSolver (const  QLblSolver &R);

    QLblSolver  &operator=( const QLblSolver  &R);


    QLblSolver (const QBigMeasure *parrBigMeasures
                , const int  QntMeas    ,const TTable_1D tblEstPrfl
                , const double Deepth  , const double Toler);

    QLblSolver(QVector<QBigMeasure> VectBigMeasures
               ,const TTable_1D tblEstPrfl, const double Deepth
            , const double Toler );

    virtual void  collectGrlsMeasure(const QBigMeasure &Meas
                                                 ,double *arrMeasure,double *arrD);

    virtual bool calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                           , double *arrNeviaz , double *dArrNeviaz_po_dX);;

    virtual bool calc_2D_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, QBigMeasure Meas
                                                        , double *arrNeviaz , double *dArrNeviaz_po_dX);

    virtual bool calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                                            , double *arrNeviaz );

    void calc_dX_po_dSgps(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dSgps);

    bool calc_FGr_and_dFGr_po_dSgps(double *arrX0,double * parrMtrx_dFGr_po_dSgps);



    void calc_dX_po_dHBeacon(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dHBeacon);

    bool calc_dFGr_po_dHBeacon(double *arrX ,double *arrdfi_po_dHBeacon);

    bool  calc_arrfi_and_arrdfi_po_dHBeacon(double *arrX
                           ,const int NUmBigMeasure,double *dHBeacon);

    bool  calc_Kgps(double *arrXRez,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kgps);

    bool  calc_Kt(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kt);

    bool  calc_dX_po_dAngSins(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_AngSins);



    bool  calc_KsinsQ(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KsinsQ_per_1_rad);

    bool  calc_Ksins_Psi_Tetta(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_Kgps_Psi_Tetta_per_1_rad);

    bool  calc_dX_po_dProfile(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_dX_po_dProfile);

    bool  calc_arrfi_and_arrdfi_po_dProfile(double *arrX
                           ,const int NUmBigMeasure,double *arrfi,double *arrdfi_po_dProfile);

    bool  calc_KProfile(double *arrX0,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KProfile_per_1_m);

    bool  calc_K_ProfilePart(double *arrX0, int numProfilePart,double *arrMtrx_dFGr_po_dX_Inv,double *arr_KPart);

    bool  calc_arrdfi_po_dProfilePart(double *arrX, const int numProfilePart
                           ,const int NUmBigMeasure,double *arrdfi_po_dProfilePart);

    virtual void calcArrPenalty(const double *arrX0,const double valPenalt, double *arrPenalt);

};

#endif // LBLSOLVER_H
