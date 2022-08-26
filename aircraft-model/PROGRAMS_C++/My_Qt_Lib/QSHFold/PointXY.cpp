#include "PointXY.h"
#include <QPoint>
#include <math.h>



QPointXY::QPointXY():QPointF()
{

}

QPointXY::QPointXY(const QPoint &point):QPointF(point)
{

}

QPointXY::QPointXY(qreal xpos, qreal ypos):QPointF(xpos, ypos)
{

}

// Конструктор копирования
 QPointXY::QPointXY (const QPointXY &R2)
 {
   this->setX(R2.x());
   this->setY(R2.y());

 }



 // оператор присваивания
 QPointXY &QPointXY::operator=(const QPointXY  &R2)
{
     this->setX(R2.x());
     this->setY(R2.y());
     return *this ;
}
 //-------------------------
 // линейное преобразование полигона  пергруженная
// INPUT:
// valAng - угол поворота
// pntSdvig - точка куда перемещается центр полигона
// valRastigenie - коэффициент растяжения
// OUTPUT:
// возвращает преобразованный полигон
QPointXY   QPointXY::LinTransform(const double  valAng , const QPointXY pntSdvig,const double valRastigenie )
{
  double arrMtxPer[4] = {0.};
  arrMtxPer[0] = cos(valAng);
  arrMtxPer[1] = -sin(valAng);
  arrMtxPer[2] = -arrMtxPer[1];
  arrMtxPer[ 3] = arrMtxPer[0];


  double x0 = this->x() * valRastigenie;
  double x1 = this->y() * valRastigenie;

  QPointXY PntTemp2(arrMtxPer[0] * x0 + arrMtxPer[1] * x1,arrMtxPer[2] * x0 + arrMtxPer[3] * x1 );



  QPointXY PntRez  = PntTemp2;// +pntSdvig;
  PntRez += pntSdvig;


  return PntRez;
}
//-----------------------------
QPointXY   QPointXY::toPntXY(const QPointF& pntF)
{
  QPointXY pntReturn;
  pntReturn.setX(pntF.x());
  pntReturn.setY(pntF.y());
  return pntReturn;
}
//-------------------------------------
// площадь треугольника
double QPointXY::calcVectS(const QPointF& P1, const QPointF& P2)
{

   return (P1.x() * P2.y() - P1.y() * P2.x())/2.;
}
//-------------------------------------
// площадь треугольника
double QPointXY::calc_dist(const QPointF& P1, const QPointF& P2)
{

   return sqrt((P1.x() - P2.x()) * (P1.x() - P2.x())
               + (P1.y() - P2.y()) * (P1.y() - P2.y()));
}
//------------------------------------------------------------------



