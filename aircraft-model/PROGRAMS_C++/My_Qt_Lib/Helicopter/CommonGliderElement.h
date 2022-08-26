#ifndef COMMONGLIDERELEMENT_H
#define COMMONGLIDERELEMENT_H

// класс описывает простой элемент планера, который не может вращаться относительно
// связанной системы координат вертолета
// в классе определяется единственная виртуальная функция
// вычисления матрицы столбцов единичных ортов
// повернутой связанной сиситемы координат руля
// относительно не повернутой
// это единичная матрица
#include "AbstractGliderElement.h"

class TAbstractGliderElement;

class TCommonGliderElement : public TAbstractGliderElement
{
public:
    TCommonGliderElement();

    TCommonGliderElement(const long  double S, const TLongPlane PlaneSvSK, const long  double Cx0, const long  double Cy0, const long   double AlfCrit
                      ,long   double *arrAirP);

    TCommonGliderElement( long double *arrInpDataGliders);

    TCommonGliderElement (const TCommonGliderElement &R);

    TCommonGliderElement operator=(TCommonGliderElement  R);


    virtual void calMtrxMwave(const long   double VAlAlfa,long   double *arrMwave);

};

#endif // COMMONGLIDERELEMENT_H
