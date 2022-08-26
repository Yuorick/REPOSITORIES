#ifndef CTRLDNMKCORRECT_M_R_H
#define CTRLDNMKCORRECT_M_R_H
//#include "CtrlFollow.h"
#include "Comp.h"
#include "Ctrl_With_Integrators.h"

//class QCtrlFollow;
class TComp;

class QCtrlDnmkCorrect_M_R : public QCtrl_With_Integrators
{
public:
    TComp mCmpArrLamb[4];
    double mTPeriodParamsRenew;
    double mTLastRenew;

    QCtrlDnmkCorrect_M_R();
    QCtrlDnmkCorrect_M_R(const QCtrlDnmkCorrect_M_R &R);
    QCtrlDnmkCorrect_M_R  &operator=( const QCtrlDnmkCorrect_M_R  &R);
    QCtrlDnmkCorrect_M_R(  const QElectMotor ElectMotor ,  QLoad*  Load
         , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
         ,const double VAlOmegaStatBegin, const double T0,const double TCur,const double valh,TComp *pCmpArrLamb
         , const double IntegratorTime, const double TPeriodParamsRenew
         , const double TLastRenew );
    virtual void fncRenewParams_M_R(const double *arrObjective,const double VAlMomOut);

    void calcNewParams(const double *arrObjective,const double VAlMomOut,const int QUantParams, double *arrParams);

    void createMeanPhVect(const double *arrObjective,const double VAlMomOut,double *arrMeanPhVect);
};

#endif // CTRLDNMKCORRECT_M_R_H
