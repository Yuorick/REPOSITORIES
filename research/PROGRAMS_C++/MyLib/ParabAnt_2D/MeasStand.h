//---------------------------------------------------------------------------

#ifndef MeasStandH
#define MeasStandH
#include "MeasStand.h"
#include "Comp.h"
class TComp;
class TMeasStand
{
// массив измерений тройки диаграмм со стенда каурова
// mRadAxeAnt - угол наклона центральной оси антенны относительно горизонта
// центральная ось проходит через середину угла, образованного осями 3 и 4 диаграмм
// mpcmparrMeas[3] - массив измерний
public:
 // угол наклона оси антенны
 double	 mRadAxeAnt ;
 // число диаграмм в ансаммбле
 int mCountDgrs;
 // массив измерений диаграмм в веере
TComp  *mpcmparrMeas ;



 __fastcall~TMeasStand();
 TMeasStand() ;


// конструктор копирования
 TMeasStand(const TMeasStand &R) ;
 TMeasStand operator=(TMeasStand  R2) ;
 // парам констр
__fastcall TMeasStand(const double RadAxeAnt, TComp *pcmparrMeas);
__fastcall TMeasStand(const double RadAxeAnt,  int CountDgrs);



};
#endif
