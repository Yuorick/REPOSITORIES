#ifndef DRIVERPARTYDRAW_H
#define DRIVERPARTYDRAW_H

#include "DriveMoveImit.h"

// розыгрыш ситуации со слежением за целью
// в упрощенном виде
// привод стоит на месте
// цель двигается или по прямой или по гармонике с периодом 4-6 с
// внешняя нагрузка - ветер + динамическая нагрузка в виде гармонически изменяющегося внешнего момента
//

class QDriveMoveImit;
class QDriverPartyDraw
{
public:
    // члены класса
    QDriveMoveImit mDriveMoveImit;

    QDriverPartyDraw();
};

#endif // DRIVERPARTYDRAW_H
