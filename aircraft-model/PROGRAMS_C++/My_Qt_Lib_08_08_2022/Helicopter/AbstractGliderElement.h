#ifndef ABSTRACTGLIDERELEMENT_H
#define ABSTRACTGLIDERELEMENT_H
#include "LongPlane.h"
#include <math.h>
#include <float.h>

class TLongPlane;

class TAbstractGliderElement
{
public:
    // площадь
    long double mS;

    // связанная система координат элемента планера в осях связанной системы координат вертолета

    TLongPlane mPlaneSvSK;

    // координаты точки приложенияя аэродинамической силы в связанной системе координат элемента планера
    long double marrAirP[3];

    // Коэффициент лобового сопротивления
    long double mCx0;

    // Коэффициент подъемной силы
    long double mCy0;

    // Критичекий угол атаки, при котором коэффициент подъемной силы максималный
    long double mAlfCrit;


    TAbstractGliderElement();

    TAbstractGliderElement (const TAbstractGliderElement &R);

    TAbstractGliderElement operator=(TAbstractGliderElement  R);

    TAbstractGliderElement(const long double S
           , const TLongPlane PlaneSvSK, const long double Cx0, const long double Cy0, const long  double AlfCrit
                                                    ,long  double *arrAirP);
    TAbstractGliderElement(long double *arrInpDataGliders);

    virtual void calMtrxMwave(const long  double VAlAlfa,long  double *arrMwave);

    void calAeroF_and_Mom(long double *arrUa, const long double VAlRo
                          ,const long  double VAlFi_rv,long  double *arrAeroF,long  double *arrAeroMom);
   // static long double DBL_SIGN(const double x)
   // {
   //  if (fabs(x) < DBL_MIN ) return 0.;
   //  return (x > 0.)?1.:0.;
  //  }

    long double fncPodFCoeff(const long  double valAlfEl);

    long double fncDerivPodFCoeff(const long double valAlfEl);

    static long  double fncTay(const long double x);

    static long  double fncDerivTay(const long double x);

   // void calAeroF_and_Mom(long double *arrUa, const long double VAlRo
                                                   //    ,const long double VAlFi_rv, long double *arrAeroF, long double *arrAeroMom);
};
void calcRotMatrX(const long  double a,long  double*matrP);
void calcRotMatrY(const long double a, long double*matrP);
void calcRotMatrZ(const long double a,long  double*matrP);
#endif // ABSTRACTGLIDERELEMENT_H
