#ifndef CTRLFOLLOWADAPT1_H
#define CTRLFOLLOWADAPT1_H
#include "CtrlFollow.h"
#include  "Filtr.h"
class QCtrlFollow;
class QFiltr;

class QCtrlFollowAdapt1:public QCtrlFollow
{
public:
    QCtrlFollowAdapt1();
    QCtrlFollowAdapt1(const QCtrlFollowAdapt1 &R);
    QCtrlFollowAdapt1  &operator=( const QCtrlFollowAdapt1  &R);
    QCtrlFollowAdapt1( const QElectMotor ElectMotor ,  QLoad*  Load
                               , const double *arrSpreadParams, const double MomOut, const double VAlTettaBegin
                               ,const double VAlOmegaStatBegin , const double T0, const double TCur
                               , const double valh,TComp *pCmpArrLamb
                       , const double *ARrFiltK_Begin,const double DispWTochka);

    QFiltr mFiltrW;
    double mDispWTochka;// скз скорости ускорения
    virtual void  correctMomOut(const double VAlTettaZv
       ,const double VAlDispTetta    , const double VAlW, const double VAlTcur);

    virtual double estimateMomOut(const double VAlTettaZv
       ,const double VAlDispTetta    , const double VAlW, const double VAlTcur);

    virtual double calcEstMomOut();

    virtual int createInputDataReport_CtrlFollow_Inherited(wchar_t*FileName, const bool bHeader);

};

#endif // CTRLFOLLOWADAPT1_H
