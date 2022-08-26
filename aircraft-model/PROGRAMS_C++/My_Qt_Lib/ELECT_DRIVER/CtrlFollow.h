#ifndef CTRLFOLLOW_H
#define CTRLFOLLOW_H


#include "Ctrl.h"
#include "MeasurmentImitator.h"
class QMeasurmentImitator;

class QCtrlFollow:public QCtrl
{
public:
    QCtrlFollow();
    QCtrlFollow(const QCtrlFollow &R);
    QCtrlFollow  &operator=( const QCtrlFollow  &R);

    QCtrlFollow( const QElectMotor ElectMotor ,  QLoad*  Load
                 , const double *arrSpreadParams,const double VAlMomOut, const double VAlTettaBegin
                 ,const double VAlOmegaStatBegin , const double T0, const double TCur, const double valh,TComp *pCmpArrLamb);

    virtual void CalcEigenValues( const double VAlTettaStat, const double VAlOmegaStat
                                  ,const double *arrC,const double VAlMomOut, TComp* cmpArrEigen);

    virtual void calcGearMtrx(const double *arrObjective,  TComp* CmpArrLamb,const double VAlMomOut, double* arrC0);



    virtual void calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU);


    static void calcSubMtrx(const double *arrA, TComp* CmpArrLamb, double* pvalr1, double* pvalr2, double* pvalr3);

   // virtual void calcFx( TEnvironment &Environment, double *arrPhVect, double *arrFx);

   // static void calcPolinom_Coef_4Degree_Vieta( TComp*CmpArrLamb, double *pb, double *pc, double *pd, double *pe);

    void calc_df_po_dx_plus_BC_Follow(const double VAlTettaStat, const double VAlOmegaStat
                          ,const double *arrC,const double VAlMomOut, double *arrFGr);

    void fill_df_po_px_and_mtrxB_Follow(const double *arrStatinaryPhVect
                     ,double *arr_dF_po_dx,double * arr_dF_po_dW);

    virtual void correctMomOut(const double VAlTettaZv
    ,const double VAlDispTetta    , const double VAlW, const double VAlTcur);

    virtual double estimateMomOut(const double VAlTettaZv
       ,const double VAlDispTetta    , const double VAlW, const double VAlTcur);

    virtual double calcEstMomOut();

    virtual int createInputDataReport_Ctrl_Inherited(wchar_t*FileName, const bool bHeader);

    virtual int createInputDataReport_CtrlFollow_Inherited(wchar_t*FileName, const bool bHeader);



    };




#endif // CTRLFOLLOW_H
