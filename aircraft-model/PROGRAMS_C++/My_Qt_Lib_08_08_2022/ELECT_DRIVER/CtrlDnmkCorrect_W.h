#ifndef CTRLDNMKCORRECT_W_H
#define CTRLDNMKCORRECT_W_H
//#include "CtrlFollow.h"
#include "Ctrl_With_Integrators.h"
//class QCtrlFollow;

class QCtrlDnmkCorrect_W : public QCtrl_With_Integrators
{
public:
    double mTPeriodRenew_dU;
    double mTLastRenew_dU;
    double marrW[2];

    //-----------------------------------------
    QCtrlDnmkCorrect_W();

    QCtrlDnmkCorrect_W(const QCtrlDnmkCorrect_W &R);

    QCtrlDnmkCorrect_W  &operator=( const QCtrlDnmkCorrect_W  &R);

    QCtrlDnmkCorrect_W( const QElectMotor ElectMotor ,  QLoad*  Load
                         , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                         ,const double VAlOmegaStatBegin , const double T0,const double TCur,const double valh,TComp *pCmpArrLamb
                         , const double IntegratorTime,const double TPeriodRenew_dU,const double TLastRenew_dU);

    virtual void addW(const double *arrObjective,const double VAlMomOut,double *arrU);

    void calcCorrectingDeltaW(const double *arrObjective,const double VAlMomOut, double *arrW);
};

#endif // CTRLDNMKCORRECT_W_H
