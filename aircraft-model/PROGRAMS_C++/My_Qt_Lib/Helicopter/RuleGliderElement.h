#ifndef RULEGLIDERELEMENT_H
#define RULEGLIDERELEMENT_H
// класс описывает элемент планера типа руля
// в классе определяется единственная виртуальная функция
// вычисления матрицы столбцов единичных ортов
// повернутой связанной сиситемы координат руля
// относительно не повернутой

#include "AbstractGliderElement.h"

class TAbstractGliderElement;

class TRuleGliderElement : public TAbstractGliderElement
{
public:
    TRuleGliderElement();

    TRuleGliderElement(const long  double S, const TLongPlane PlaneSvSK, const long  double Cx0, const long double Cy0, const long double AlfCrit
                      ,long  double *arrAirP);

    TRuleGliderElement(long double *arrInpDataGliders);

    TRuleGliderElement (const TRuleGliderElement &R);

    TRuleGliderElement operator=(TRuleGliderElement  R);


    virtual void calMtrxMwave(const long  double VAlAlfa,long  double *arrMwave);

};
#endif // RULEGLIDERELEMENT_H
