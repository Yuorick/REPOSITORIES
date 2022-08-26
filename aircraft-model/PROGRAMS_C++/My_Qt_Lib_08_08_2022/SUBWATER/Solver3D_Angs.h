#ifndef Solver3D_Angs_H
#define Solver3D_Angs_H
#include "Solver2D_Angs.h"



class QSolver3D_Angs : public QSolver2D_Angs
{
public:
   QSolver3D_Angs();

    QSolver3D_Angs (const  QSolver3D_Angs &R);

    QSolver3D_Angs  &operator=( const QSolver3D_Angs  &R);

    QSolver3D_Angs (const QBigMeasure *parrBigMeasures
                    , const int  QntMeas    ,const TTable_1D tblEstPrfl
                    , const double Deepth
                    , const double Toler,const double *arrPAnt_PSK
                    ,const double *arrSBeacon_GSK);    

    virtual void  collectGrlsMeasure(const QBigMeasure &Meas
                                                 ,double *arrMeasure,double *arrD);


    virtual bool calc_arrNeviaz_and_dArrNeviaz_po_dX_3D( const int NUmBigMeasure
      ,const double e_est,double *arr_dV_po_dJ ,double *arrt7, double *arrNeviaz , double *dArrNeviaz_po_dX);

    virtual bool calc_arrNeviaz_3D( const int NUmBigMeasure ,const double e_est, double *arrNeviaz );

};

#endif // Solver3D_Angs_H
