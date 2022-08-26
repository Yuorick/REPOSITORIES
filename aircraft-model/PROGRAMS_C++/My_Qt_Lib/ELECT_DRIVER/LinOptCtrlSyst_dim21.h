#ifndef LINOPTCTRLSYST_dim21_H
#define LINOPTCTRLSYST_dim21_H
#include "LinOptCtrlSyst.h"

class QLinOptCtrlSyst;
class QLinOptCtrlSyst_dim21:public QLinOptCtrlSyst
{
public:

    double marr_ekt_b[4];
    double marr_t_ekt_b[4];
    double marr_b0Inv[4];
    double mk1;
    double mk2;
    bool mbInit;
    bool mbMultRoots;
    bool mbZeroRoot;
    QLinOptCtrlSyst_dim21();
    QLinOptCtrlSyst_dim21 (double *pA,
                                       double *pB,
                                       double *pC,
                                       const double T0,
                                       double *arrPhVect0,
                                       const double TCur,
                                       double *arrPhVect);

    QLinOptCtrlSyst_dim21 (const  QLinOptCtrlSyst_dim21 &R);

     QLinOptCtrlSyst_dim21  &operator=( const QLinOptCtrlSyst_dim21  &R);

     bool calcMtrx_ForFundMtrxCalculation();

  //  bool  calcSwitchPointTimes_FixedSignum0_Old(double *arrxEnd,const double VAlUqu0
    //                                        ,double *pvalt0,double *pvalt1,const int isignum0);

    bool calcSwitchPointTimes(double *arrxEnd,const double VAlUqu0
                                          ,double *pvalt0,double *pvalt1, int *isignum0);

    void calcFundMtrx(const double VAlT,double *arrFundMtrx);

    void calcFgr_and_dFgr(double *arrt, double *arrxEnd
                                                 , double *arrFGr, double *arr_dFGr_po_dt);

    void calcFgr_and_dFgr(const double VAlUqu0, const int isignum0
                      , double *arrt, double *arrxEnd, double *arrFGr, double *arr_dFGr_po_dt);

    void  calcIntegralFundMtrx(const double VAlT,double *arrIntegralFundMtrx);

    void calcDefiniteIntegralFundMtrx(const double VAla,const double VAlb,double *arrDefIntegralFundMtrx);

    void calcSpecialIntegralFundMtrx(const double VAlx1,const double VAla,const double VAlb,double *arrDefIntegralFundMtrx);

    void drag( const double VAlIntegrStep, const double VAlU0
              ,const double valt0,const double valt1);

    void calcRightPart(const double VAlU, double *arrRightPart);

    void goToNextTime(const double VAlNextTime, const double VAlU0);

    void calcRightPart(const long double VAlU,  long double *arrPhVect, long  double *arrRightPart);

    void goToNextTimeEiler(const double VAlNextTime, const double VAlU0
                                                  , const double VAlStepIntegr);

    bool calcSwitchPointTimes_FixedSignum0(double *arrxEnd,const double VAlUqu0
                                          ,double *pvalt0,double *pvalt1, const int isignum0);

    void calcFgr_and_dFgr_Full(const double VAlUqu0, const int isignum0
                      , double *arrt, double *arrxEnd, double *arrFGr, double *arr_dFGr_po_dt);


};

#endif // LINOPTCTRLSYST_dim21_H
