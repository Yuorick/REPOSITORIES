#ifndef BRLS_H
#define BRLS_H
#include "Load.h"
#include "Segment.h"
#include "Environment.h"
#include "Plane.h"
class QSegment;
class TEnvironment;
class TPlane;

class QBrls:public QLoad
{
public:

    QBrls();
    QBrls (const  QBrls &R);
    QBrls  &operator=( const QBrls  &R);
    QBrls ( const double  Cv, const double   Cx, const double  VAlM, const double  VAlLengthTotal,const double  VAlHeight
            ,const double  Max_dOm_po_dt, const QSegment SgmBarrier, const TPlane Plane);
// полная длина штанги (диаметр вращения)
    double mLength;
    double mHeight;
    //преграда
    QSegment mSgmBarrier;
    // плоскость вращения БРЛС в ГСК (система координат СКБРЛС)
    TPlane mPlane;

    // если нет барьера, то момент аэродинамического сопротивления
    // вычисляется просто и быстро. Поэтому, этот случай следует выделить
    // чтобы иметь возможность эконмить  время

    // признак наличия барьера рядом в БРЛС
    bool mbBarrier;



    virtual void  setLength(const double  VAlLength);
    virtual  double  getLength();

    virtual double  calcAirResistMom( TEnvironment &Environment,const double ValTet, const double VAlOm);

    double  calcCOmegaX ();

    virtual double  calc_dAirResistMom_po_dOmega( TEnvironment &Environment,const double ValTet, const double VAlOm);

    double  calcAirResistFuncForBlade(const TEnvironment Environment
    ,const double ValTet, const double VAlOm
     ,double (*pf)(const double Cx,const double VAlOm,const double  valW,const double  VAla ,const double  VAlb,const double  VAlh));

    double  calcAirResistFunc_With_Limits( const double VAlOm,  const double valW
        ,  const double ValLowLim,  const double ValUpperLim
        ,double (*pf)(const double Cx,const double VAlOm,const double  valW,const double  VAla ,const double  VAlb,const double  VAlh));

    inline double sign__(const double x)
    {
        return (x >=0)?1:-1;
    }

    virtual double calc_dMa_po_Cx(const double VAlOm);

    virtual double  calc_d2Ma_po_dOm_po_Cx(const double VAlOm);

    virtual int createInputDataReport(wchar_t*FileName, const bool bHeader);


};
double calcIntegralAirMom(const double Cx,const double VAlOm,const double  valW,const double  VAla
                                  ,const double  VAlb,const double  VAlh);
double  calcIntegral_dAirMom_po_dOmega(const double Cx, const double VAlOm,const double  valW,const double  VAla
                                  ,const double  VAlb,const double  VAlh);




#endif // BRLS_H
