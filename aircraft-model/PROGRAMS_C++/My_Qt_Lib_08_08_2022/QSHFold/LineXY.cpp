#include "LineXY.h"
#include <QPointF>
//#include <QLine>
#include <QVector>
#include <math.h>
#include "PointXY.h"
#include <QRectF>
#include "URPolyLine.h"
#include "UrPointXY.h"




QLineXY::QLineXY():QLineF()
{
}
//-------------------------------------

// конструктор копирования
QLineXY :: QLineXY (const  QLineXY &R)
 {
   this->setP1(R.p1());
   this->setP2(R.p2());
 }

 // оператор присваивания
QLineXY  &QLineXY::operator=( const QLineXY  &R)
 {
      if(this == &R)
      {
          return *this;
      }
      this->setP1(R.p1());
      this->setP2(R.p2());

     return *this ;
 }

  // парам конструктор 1
 QLineXY:: QLineXY (const QPointF &p1, const QPointF &p2)
  :QLineF (p1, p2)
 {

 }
//----------------------------------------------
 QLineXY::QLineXY(qreal x1, qreal y1, qreal x2, qreal y2)
     :QLineF (x1,  y1,  x2,  y2)
 {

 }
 //-----------------------------------------
 QLineXY::QLineXY(const QLine &line):QLineF (line)
 {
 }
 //------------------------------------------------
 QPointXY QLineXY::getPXY1()
 {
     return QPointXY(this->p1().x(),this->p1().y());
 }
 //------------------------------------------------
 QPointXY QLineXY::getPXY2()
 {
     return QPointXY(this->p2().x(),this->p2().y());
 }
 //------------------------------
 //создание вектора - стрелки
 void QLineXY::createVectArrow(const double xbegin, const double ybegin
                               , const double e_directx,  const double e_directy
                               ,const double length0,const double valLength
                               , const double valAng, QVector<QLineF> *vctArrow)
 {
    QVector<QLineF> vectRez(3);
    QPointF pointBegin(xbegin,ybegin);
    QPointF pointEnd(xbegin + length0 * e_directx, ybegin + length0 * e_directy);
    QLineF ln0(pointBegin, pointEnd);
    vectRez.replace(0, ln0);

    double d = ln0.length();
    QPointF pnt00 = ln0.pointAt((d -valLength )/d) ;

    QLineF lnt1(QPointF(0.,0.), pointEnd - pointBegin);
    QLineF lnt2 =lnt1.normalVector();
    QLineF lnt3 = lnt2.unitVector();

    QPointF pntDel = lnt3.p2();

    QPointF pnt1 = pnt00 +valLength * tan(valAng) *pntDel;
    QPointF pnt2 = pnt00 -valLength * tan(valAng) *pntDel;

    QLineF ln1(pnt1, pointEnd);
    QLineF ln2(pnt2, pointEnd);

    vectRez.replace(1, ln1);
    vectRez.replace(2, ln2);

    /*
    QPointXY pntTr0(-valLength, valLength*tan(valAng));
    QPointXY pntTr1(-valLength, -valLength*tan(valAng));

    QPointXY point3, point5;
    if(fabs(pointBegin.x() - pointEnd.x() ) < 0.0001)

    {
        point3 = QPointXY(pointEnd.x() -valLength * tan(valAng) , pointEnd.y() -  Sign(pointEnd.y() -pointBegin.y())* valLength);
        point5 = QPointXY(pointEnd.x() +valLength * tan(valAng) , pointEnd.y() +  Sign(pointEnd.y() -pointBegin.y())* valLength);
    }
    else
    {
     double  valAngRot = atan2(pointEnd.y()- pointBegin.y(), pointEnd.x()- pointBegin.x());
     point3 = pntTr0.LinTransform(  valAngRot , pointEnd,1. );
     point5 = pntTr1.LinTransform(  valAngRot , pointEnd,1. );
    }
   vectRez.replace(1, QLineF(QPointF(point3), pointEnd));
   vectRez.replace(2, QLineF(QPointF(point5), pointEnd));
*/
   (*vctArrow) = vectRez;
 }

//-----------------------

 void QLineXY::rastiagenie(QVector<QLineF> *vctLine, const double coef)
 {
     for (int i =0; i < vctLine->size(); ++i)
     {
        QLineF temp = (*vctLine).at(i);
        QPointF pnt = temp.p1();
        pnt *= coef;
        temp.setP1(pnt);

        pnt = temp.p2();
        pnt *= coef;
        temp.setP2(pnt);
        (*vctLine).replace( i, temp);
     }

 }


 //---------------------------------
void QLineXY::stretchVct(QVector<QLineF> *vctLine,const double &scalx,const double &scaly)
{
    for (int i = 0; i < vctLine->size(); ++i)
    {
        QLineF ln0 = vctLine->at(i);
        stretch(ln0, scalx, scaly);
        vctLine->replace(i,ln0);
    }
}
//--------------------------------------------
void QLineXY::stretch(QLineF& ln0,const double &scalx,const double &scaly)
{
    QPointF pnt = ln0.p1();
    QPointF pnt1(pnt.x() *scalx, pnt.y() *scaly);
    ln0.setP1(pnt1);

    pnt = ln0.p2();
    QPointF pnt2(pnt.x() *scalx, pnt.y() *scaly);
    ln0.setP2(pnt2);
}
 //---------------------------------------
 //QVector<QLineF> mvctRelativeVelo
 double QLineXY::Sign(const double x)
 {
     if (x > 0.)
     {
        return 1.;
     }
     if (x < 0.)
     {
        return -1.;
     }
     return 0.;
 }
//-----------------------------
 void QLineXY::translateVct(QVector <QLineF> *vctLine, QPointF &pntOffset)
 {
     for (int i =0; i < (*vctLine).size(); ++i)
     {
       QLineF ln = (*vctLine).at(i);
       ln.translate(pntOffset);
       (*vctLine).replace(i,ln);
     }
 }
 //------------------------------
 //создание вектора линий сетки
 void QLineXY::createVectNet(const double xmin, const double xmax
                               ,const double ymin, const double ymax
                               , const double xstep,  const double ystep
                               , QVector<QLineF> *vctNet)
 {

     int kx = int((xmin - 0.00001)/ xstep);
     double minx = ((double)kx )* xstep;
      kx = int((xmax + 0.00001)/ xstep);
      double maxx = ((double)kx )* xstep;
      int n = int((maxx - minx + 0.000001)/ ystep) +1;
     vctNet->resize(n);


     int numx = 0;
     for (int i =0; i < n; ++i)
     {
         double xcur = minx + ((double )i )*xstep;
         if (xcur > (xmax +0.0000001))
         {
             break;
         }
         QPointF pnt1(xcur, ymin);
         QPointF pnt2(xcur, ymax);
       vctNet->replace(numx, QLineF(pnt1, pnt2));
       numx++;
     }
     vctNet->resize(numx);


     int ky = int((ymin - 0.00001)/ ystep);
     double miny = ((double)ky )* ystep;
      ky = int((ymax + 0.00001)/ ystep);
      double maxy = ((double)ky )* ystep;
      int m = int((maxy - miny + 0.000001)/ ystep) +1;
      vctNet->resize(numx + m);

     int numy =0;
     for (int i =0; i < m; ++i)
     {
         double ycur = miny + ((double )i )*ystep;
         if (ycur > (ymax + 0.0000001))
         {
             break;
         }
         QPointF pnt1(xmin, ycur);
         QPointF pnt2(xmax, ycur);
       vctNet->replace(numx + numy, QLineF(pnt1, pnt2));
       numy++;
     }
     vctNet->resize(numx + numy);

 }
 //------------------------------------------

 //создание вектора линий сетки
 void QLineXY::createVectAxes(const double xmin, const double xmax
                               ,const double ymin, const double ymax
                              ,const double arrowLength, const double valAng
                               , QVector<QLineF> *VctAxes)
 {
      double e_directx = 1.;
      double e_directy = 0.;

     createVectArrow(xmin, 0.,  e_directx,  e_directy
                     ,xmax - xmin, arrowLength,  valAng, VctAxes);

     QVector<QLineF> vctArrowY;

     e_directx = 0.;
     e_directy = 1.;
     createVectArrow(0., ymin, e_directx,  e_directy
                                     ,ymax-ymin, arrowLength
                                     ,  valAng, &vctArrowY);

     //(*VctAxes).append(vctArrowY);
      (*VctAxes) += vctArrowY;


 }
//-----------------------------------------
QRectF  QLineXY::boundBox(QLineF &ln)
{
   double minx =0., maxx = 0., miny = 0., maxy = 0.;
   if((ln.p1().x()<ln.p2().x()))
   {
       minx =ln.p1().x();
       maxx =ln.p2().x();
   }
   else
   {
       minx =ln.p2().x();
       maxx =ln.p1().x();
   }


   if((ln.p1().y()<ln.p2().y()))
   {
       miny =ln.p1().y();
       maxy =ln.p2().y();
   }
   else
   {
       miny =ln.p2().y();
       maxy =ln.p1().y();
   }

   return QRectF(QPointF(minx , miny),QPointF(maxx , maxy));
}
//---------------------------------------------
bool   QLineXY::isIntersected(const QRectF &Rect)
{
    QLineF arrLn[4];
    arrLn[0] =QLineF(Rect.topLeft(), Rect.topRight());
    arrLn[1] =QLineF(Rect.topRight(), Rect.bottomRight());
    arrLn[2] =QLineF(Rect.bottomRight(), Rect.bottomLeft());
    arrLn[3] =QLineF(Rect.bottomLeft(), Rect.topLeft());


    for (int i =0; i < 4; ++i)
    {
    QPointF intersectionPoint;
    if(intersect(arrLn[1], &intersectionPoint)==QLineF::BoundedIntersection)
    return true;
    }
    return false;

}
//-------------------------------------------------
//---------------------------------------------
bool   QLineXY::isIntersected(const QVector<QLineF> &vctLine, QRectF &Rect)
{


    for (int i =0; i < vctLine.size(); ++i)
    {
    QLineXY ln = createLineXY(vctLine.at(i));
    if(ln.isIntersected(Rect))
    return true;
    }
    return false;

}

//----------------------------------------------

QLineXY QLineXY::createLineXY(const QLineF &ln)
{
  QLineXY lnFreturn;
  lnFreturn.setP1(ln.p1());
  lnFreturn.setP2(ln.p2());

  return lnFreturn;
}
//--------------------------------------
TURPolyLine QLineXY::createTURPln(const QLineF &ln)
{
  TURPointXY  pnt1,  pnt2;

  pnt1.X = ln.p1().x();
  pnt1.Y = ln.p1().y();

   pnt2.X = ln.p2().x();
   pnt2.Y = ln.p2().y();
   TURPolyLine lnTreturn( pnt1,  pnt2);
  return lnTreturn;
}

//--------------------------------------
TURPolyLine QLineXY::createTURPln(const QVector<QLineF> &vctLine)
{
    int *iarrParts = new int[vctLine.size()];
    for (int i = 0; i <vctLine.size();++i )
    {
      iarrParts[i] = 2 * i;
    }
    TURPointXY *arrPoints = new TURPointXY[ 2 * vctLine.size()];
    for (int i = 0; i <vctLine.size();++i )
    {
        arrPoints [2 * i].X = vctLine.at(i).p1().x();
        arrPoints [2 * i].Y = vctLine.at(i).p1().y();


        arrPoints [2 * i +1].X = vctLine.at(i).p2().x();
        arrPoints [2 * i +1].Y = vctLine.at(i).p2().y();
    }
    TURPolyLine plnReturn(vctLine.size(),2 * vctLine.size(),iarrParts
               ,arrPoints) ;

    delete []iarrParts;
    delete []arrPoints;


  return  plnReturn;

}
