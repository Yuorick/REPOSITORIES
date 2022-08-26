#ifndef CTRL_WITH_INTEGRATORS_H
#define CTRL_WITH_INTEGRATORS_H

#include "CtrlFollow.h"

//class QCtrlFollow;

class QCtrl_With_Integrators : public QCtrlFollow
{
public:
    // фазовый вектор интеграторов отклонений фазового вектора от цели
    double marrIntegrators[QVARS];
     // текущее время привязки
    double mTIntegrCur;
    // постоянная времени
    double mIntegratorTime;



    QCtrl_With_Integrators();

    QCtrl_With_Integrators(const QCtrl_With_Integrators &R);

    QCtrl_With_Integrators  &operator=( const QCtrl_With_Integrators  &R);

    QCtrl_With_Integrators( const QElectMotor ElectMotor ,  QLoad*  Load
                            , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                            ,const double VAlOmegaStatBegin , const double T0,const double TCur, const double valh
                              , TComp *pCmpArrLamb,const double IntegratorTime );

    virtual void calcCurU(const double *arrObjective,const double VAlTObjective,const double VAlMomOut, double *arrU);

    void recalcIntegrators(const double *arrObjective,const double VAlMomOut);

    virtual void addW(const double *arrObjective,const double VAlMomOut,double *arrU);

    virtual void fncRenewParams_M_R(const double *arrObjective,const double VAlMomOut);
};

#endif // CTRL_WITH_INTEGRATORS_H
