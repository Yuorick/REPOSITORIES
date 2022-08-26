#ifndef CTRLVELO_H
#define CTRLVELO_H
#include "Ctrl.h"

class QCtrl;
class QMeasurmentImitator;
class QCtrlVelo:public QCtrl
{
public:

    QCtrlVelo();
    QCtrlVelo(const QCtrlVelo &R);
    QCtrlVelo  &operator=( const QCtrlVelo  &R);

    QCtrlVelo( const QElectMotor ElectMotor ,  QLoad*  Load
               , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
               ,const double VAlOmegaStatBegin , const double T0,const double TCur, const double valh,TComp *pCmpArrLamb);

    virtual void CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                                  ,const double *arrC,const double VAlMomOut, TComp* cmpArrEigen);

    virtual void calcGearMtrx(const double *arrObjective, TComp* CmpArrLamb,const double VAlMomOut, double* arrC);


    virtual void calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU);

    void fill_df_po_px_and_mtrxB_Velocity(const double *arrStatinaryPhVect
                     ,double *arr_dF_po_dx,double * arr_dF_po_dW);

    void calc_df_po_dx_plus_BC_Velo(const double VAlOmegaStat
                          ,const double *arrC,const double VAlMomOut, double *arrFGr);


};

#endif // VELOCTRL_H
