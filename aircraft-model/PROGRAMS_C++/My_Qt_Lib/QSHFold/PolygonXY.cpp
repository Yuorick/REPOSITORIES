#include "PolygonXY.h"
#include <QLineF>
#include <QPointF>
#include <QPolygonF>
#include <float.h>

#include "URPolygon.h"
#include "UrPointXY.h"
#include "URPolyLine.h"

class QLineXY;
class QPointXY;

QPolygonXY::QPolygonXY():QPolygonF()
{
}


QPolygonXY::QPolygonXY(int size):QPolygonF( size){}

QPolygonXY::QPolygonXY(const QVector<QPointF> &points):QPolygonF(points){}

QPolygonXY::QPolygonXY(QVector<QPointF> &&v):QPolygonF(v){}

QPolygonXY::QPolygonXY(const QRectF &rectangle):QPolygonF(rectangle){}

QPolygonXY::QPolygonXY(const QPolygon &polygon):QPolygonF(polygon){}

QPolygonXY::QPolygonXY(const QPolygonF &polygon):QPolygonF(polygon){}

QPolygonXY::QPolygonXY(QPolygonF &&other):QPolygonF(other){}
//----------------------------------------------------
//----------------------------------------------------

void QPolygonXY::calcBoundBox()
{
     double xMin = 1000000000 ;
     double xMax = -1000000000 ;
     double yMin =  1000000000 ;
     double yMax =  - 1000000000 ;
     for (int i =0; i < size(); i++)
     {
       if (at(i).x() > xMax) xMax =  at(i).x();
       if (at(i).x() < xMin) xMin =  at(i).x();
       if (at(i).y() > yMax) yMax =  at(i).y();
       if (at(i).y() < yMin) yMin =  at(i).y();

     }
    Box[0] =  xMin;
    Box[1] =  yMin;
    Box[2] =  xMax ;
    Box[3] =  yMax;

}
//----------------------------------------------

QPolygonXY::QPolygonXY(const QVector<QPointXY> &points)
{
    this->resize(points.size());
    for (int i =0; i < points.size(); ++i)
    {
    this->replace(i,points.at(i));
    }
}

//------------------------------------
QPolygonXY::QPolygonXY(const QVector<QPointXY> points)
{
    this->resize(points.size());
    for (int i =0; i < points.size(); ++i)
    {
    this->replace(i,points.at(i));
    }
}
//----------------------------------------------

QPolygonF QPolygonXY::createPlgF()
{
  QPolygonF plgFreturn;
  plgFreturn.resize(this->size());
  for (int i =0; i < this->size(); ++i)
  {
  plgFreturn.replace(i,this->at(i));
  }
  return plgFreturn;
}
//----------------------------------------------

QPolygonXY QPolygonXY::createPlgXY(const QPolygonF &plgF)
{
  QPolygonXY plgFreturn;
  plgFreturn.resize(plgF.size());
  for (int i =0; i < plgF.size(); ++i)
  {
  plgFreturn.replace(i,plgF.at(i));
  }
  return plgFreturn;
}
//------------------------------------

QPolygonXY   QPolygonXY::LinTransform(const double  valAng , const QPointXY pntSdvig
                                      ,const double valRastigenie )
{
  QPolygonF plgFreturn = *this;
  for (int i =0; i < this->size(); ++i)
  {
      QPointXY pnt = QPointXY::toPntXY(this->at(i));
      QPointXY pnt1 = pnt.LinTransform( valAng ,  pntSdvig
                       ,valRastigenie );
      plgFreturn.replace(i,pnt1);
  }
  return plgFreturn;

}
//------------------------------------

void   QPolygonXY::LinTransform_(const double  valAng , const QPointXY pntSdvig
                                      ,const double valRastigenie )
{

  for (int i =0; i < this->size(); ++i)
  {
      QPointXY pnt = QPointXY::toPntXY(this->at(i));
      QPointXY pnt1 = pnt.LinTransform( valAng ,  pntSdvig
                       ,valRastigenie );
      this->replace(i,pnt1);
  }


}
//----------------------------------

TURPolygon QPolygonXY::createTURPlg()
{
    TURPointXY *arrpnt = new TURPointXY[this->size()];

    for (int i =0; i < this->size(); ++i)
    {
        QPointF pntF = (*this).at(i);
        arrpnt[i].X = pntF.x();
        arrpnt[i].Y = this->at(i).y();

    }
    TURPolygon plgTreturn(this->size(), arrpnt);
    delete []arrpnt;
    return plgTreturn;
}
//----------------------------------//----------------------------------

TURPolygon QPolygonXY::createTURPlg(QVector <QPolygonF> &vctPlg)
{
    TURPolygon *parrPlg = new TURPolygon [vctPlg.size()];
    for (int i =0; i < vctPlg.size(); ++i)
    {
        QPolygonXY plg = createPlgXY(vctPlg.at(i));
      parrPlg[i] =  plg.createTURPlg();
    }

    TURPolygon plgReturn( parrPlg, vctPlg.size());
    delete []parrPlg;
    return plgReturn;
}
//----------------------------------
TURPolygon QPolygonXY::createTURPlg(const QPolygonF &plg)
{
    TURPointXY *arrpnt = new TURPointXY[plg.size()];

    for (int i =0; i < plg.size(); ++i)
    {
        QPointF pntF = plg.at(i);
        arrpnt[i].X = plg.at(i).x();
        arrpnt[i].Y = plg.at(i).y();

    }
    TURPolygon plgTreturn(plg.size(), arrpnt);
    delete []arrpnt;
    return plgTreturn;
}
//----------------------------------

TURPolyLine QPolygonXY::createTURPln()
{
    TURPointXY *arrpnt = new TURPointXY[this->size()];

    for (int i =0; i < this->size(); ++i)
    {
        QPointF pntF = (*this).at(i);
        arrpnt[i].X = pntF.x();
        arrpnt[i].Y = this->at(i).y();

    }
    TURPolyLine plnTreturn(arrpnt,this->size());
    delete []arrpnt;
    return plnTreturn;
}
//------------------------------------------
// создание полилинии нижнего фронта полигона
QPolygonXY QPolygonXY::createPlnLowFront()
{
    //1
    if(!isClockwise())
    {
       flipPoints();
    }
    ///

    //2
    int iRight = findNumberUtmostRightPoint();

    int iLeft = findNumberUtmostLeftPoint();
    ///

    // 3
    QPolygonXY plgReturn(size());
     int quantPoints = 0;
     if (iLeft >= iRight)
     {
         for (int i = iRight; i <= iLeft; ++i)
         {
           plgReturn.replace(quantPoints, at(i));
           quantPoints++;
         }
     }
     else
     {
        for (int i = iRight; i < size(); ++i)
        {
        plgReturn.replace(quantPoints, at(i));
        quantPoints++;
        }

        for (int i = 1; i <= iLeft; ++i)
        {
        plgReturn.replace(quantPoints, at(i));
        quantPoints++;
        }
     }

     ///

     // 4
     plgReturn.resize(quantPoints);

    return plgReturn;

}

//------------------------------------------
// создание полилинии нижнего фронта полигона
QPolygonXY QPolygonXY::createPlnUpFront()
{
    //1
    if(!isClockwise())
    {
       flipPoints();
    }
    ///

    //2
    int iRight = findNumberUtmostRightPoint();

    int iLeft = findNumberUtmostLeftPoint();
    ///

    // 3
    QPolygonXY plgReturn(size());
     int quantPoints = 0;
     if (iRight >= iLeft)
     {
         for (int i = iLeft; i <= iRight; ++i)
         {
           plgReturn.replace(quantPoints, at(i));
           quantPoints++;
         }
     }
     else
     {
        for (int i = iLeft; i < size(); ++i)
        {
        plgReturn.replace(quantPoints, at(i));
        quantPoints++;
        }

        for (int i = 1; i <= iRight; ++i)
        {
        plgReturn.replace(quantPoints, at(i));
        quantPoints++;
        }
     }

     ///

     // 4
     plgReturn.resize(quantPoints);

    return plgReturn;

}
//-----------------------------------------
double QPolygonXY::calcSq()
{
    double S =0 ;
    if (!(this->isClosed()))
    {
        return 0.;
    }
    for (int i =0; i < this->size() -1; i++)
    {
       S = S + QPointXY::calcVectS(this->at(i),this->at(i +1));
    }
    return -S;
}
//-----------------------------------------
double QPolygonXY::calcSq(const QPolygonF &plg)
{
    double S =0 ;
    if (!(plg.isClosed()))
    {
        return 0.;
    }
    for (int i =0; i < plg.size() -1; i++)
    {
       S = S + QPointXY::calcVectS(plg.at(i),plg.at(i +1));
    }
    return -S;
}
//-------------------------------------------
bool QPolygonXY::isClockwise()
{
    if (this->calcSq() >0.)
    {
        return true;
    }
    return false;
}
//-------------------------------------------

// изменение порядка следования вершин в частичном полигоене с ноером n
void QPolygonXY::flipPoints()
{
    int quanP = this->size();

    for (int i =0; i < (quanP /2); i++)
    {
       QPointF pnt = this->at(quanP -1 - i);
       this->replace(quanP -1 - i, this->at(i));
       this->replace(i, pnt);
    }

}
//--------------------------------------
// поиск номера самой правой вершины
int QPolygonXY::findNumberUtmostRightPoint()
{
    int ireturn =-1;
    double val = -DBL_MAX;
    for (int i =0; i < size(); ++i)
    {
       if (at(i).x() > val)
       {
         val =  at(i).x();
         ireturn = i;
       }
    }
    return ireturn;
}
//
//--------------------------------------
// поиск номера самой верхней вершины
int QPolygonXY::findNumberUtmostUpPoint()
{
    int ireturn =-1;
    double val = -DBL_MAX;
    for (int i =0; i < size(); ++i)
    {
       if (at(i).y() > val)
       {
         val =  at(i).y();
         ireturn = i;
       }
    }
    return ireturn;
}
//--------------------------------------
// поиск номера самой верхней вершины
int QPolygonXY::findNumberUtmostLowPoint()
{
    int ireturn =-1;
    double val = DBL_MAX;
    for (int i =0; i < size(); ++i)
    {
       if (at(i).y() < val)
       {
         val =  at(i).y();
         ireturn = i;
       }
    }
    return ireturn;
}
//--------------------------------------
// поиск номера самой левой вершины
int QPolygonXY::findNumberUtmostLeftPoint()
{
    int ireturn =-1;
    double val = DBL_MAX;
    for (int i =0; i < size(); ++i)
    {
       if (at(i).x() < val)
       {
         val =  at(i).x();
         ireturn = i;
       }
    }
    return ireturn;
}
//---------------------------------------------------------
bool QPolygonXY::calcDistBetweenFrontLines(QPolygonXY& plnUp,QPolygonXY& plnLow, double &valDist)
{
    QRectF UpRect =plnUp.boundingRect();
    QRectF LowRect = plnLow.boundingRect();
    if (UpRect.intersects(LowRect))
    {
        return false;
    }
    valDist = DBL_MAX;
    bool breturn = false;
    // расстояние по вертикали  от вершин верхнего полигона до нижнего
    for (int i = 0; i < plnUp.size(); ++i)
    {
        double x = plnUp.at(i).x();
       if ((x < LowRect.left())  || (x > LowRect.right()) )
       {
           continue;
       }
       double valDistCur = 0.;

       QPointF pntt = plnUp.at(i);
       if(plnLow.calcVertDistFromPoint(pntt, valDistCur))
       {
           if (valDistCur < valDist)
           {
               breturn = true;
               valDist = valDistCur;
           }
       }
      }

    // расстояние по вертикали  от вершин нижнего полигона до верхнего
    for (int i = 0; i < plnLow.size(); ++i)
    {
        double x = plnLow.at(i).x();
       if ((x < UpRect.left())  || (x > UpRect.right()) )
       {
           continue;
       }
       double valDistCur = 0.;

       if(plnUp.calcVertDistFromPoint(plnLow.at(i), valDistCur))
       {
           if (valDistCur < valDist)
           {
               breturn = true;
               valDist = valDistCur;
           }
       }
      }

    return breturn;

}
//-------------------------------------
bool QPolygonXY::calcSdvig(QPolygonXY& plnUp,QPolygonXY& plnLow, double &valSdvig)
{
  valSdvig = DBL_MAX;
  bool breturn = false;
  for (int i = (plnUp.size() -1); i >=0;--i  )
  {
     double valDistCur = -1;
     bool b = false;
     if(plnLow.calcDistNearestRightPoint(plnUp.at(i).x(), valDistCur))
     {
        if(!breturn)
        {
        breturn = true;
        }
        if (valDistCur < valSdvig)
        {
        valSdvig =  valDistCur;
        }
     }
    else
     {
       continue;
     }

  }
 return breturn;
  }
//--------------------------------------
bool QPolygonXY::calcDistNearestRightPoint(const double x, double &valDistCur)
{
    bool breturn = false;
    valDistCur = DBL_MAX;
    for (int i = 0; i < size(); ++i)
    {
      if (x > (at(i).x() - 0.000001))
      {
          continue;
      }
      breturn = true;
      double d = at(i).x() - x;

      if(d < DBL_MIN )
      {
          int yy=0;
      }

      if (d < valDistCur)
      {
          valDistCur = d;
      }
    }
    return  breturn;
}
//---------------------------------------------------------
bool QPolygonXY::calcVertDistFromPoint(const QPointF &pntInp, double &valDist)
{
   bool breturn = false;
   valDist = DBL_MAX;
   double x = pntInp.x();
   for (int j = 0; j < size() -1;++j )
   {
      if ((x - at(j).x())*(x - at(j +1).x()) <=0.)
      {
          breturn = true;
          QLineF ln0(pntInp, QPointF(x, pntInp.y() + 10.));
          QLineF ln1(at(j), at(j +1));
          QPointF intersectionPoint;
          ln1.intersect(ln0, &intersectionPoint);
          QPointF pnt = intersectionPoint - pntInp;
          double disttemp = sqrt(pnt.x() * pnt.x() + pnt.y() * pnt.y());
          if (disttemp < valDist)
          {
              valDist = disttemp;
          }
      }
   }
   return breturn;

}

//---------------------------------------------------
// Вычисление первой точки пересечения линии, заданной начальной точкой
// pntNose и вектором arrVIzd_TargSvSK0
// OUTPUT:
// pntOutput - первая точка пересечения
bool QPolygonXY::findFirstIntersectPoint(const QPointF &pntNose
                                         ,const double* arrVIzd_TargSvSK0, QPointF &pntOutput)
{
    bool breturn  = false;
    // ед вектор направления     arre[2]
    double arre[2] = {0.};
    arre[0] = arrVIzd_TargSvSK0[0];
    arre[1] = arrVIzd_TargSvSK0[1];
    double norm = sqrt( arre[0] *  arre[0]  + arre[1] *  arre[1] );
    if (norm < DBL_MIN)
    {
        return false;
    }
    arre[0] = arre[0]/norm;
    arre[1] = arre[1]/norm;
    ///

    // начальная точка нос
    QPointF pnt0 = pntNose;
    pnt0 += QPointF(-arre[0] * 1.E9, -arre[1] * 1.E9);

    //
    QPointF pnt1 = pntNose;
    pnt1 += QPointF(arre[0] * 1.E9, arre[1] * 1.E9);

    const QLineF lnDirect(pnt1,pnt0);
    ///
    double dist = DBL_MAX;
    for (int i =0; i < size()-1; ++i)
    {
       const QLineF lnCur(at(i), at(i +1));
       QPointF intersectionPoint;
       switch(lnCur.intersect(lnDirect, &intersectionPoint))
       {
       case QLineF::UnboundedIntersection:
       case QLineF::NoIntersection:
           continue;
           break;
       default:
           break;
       }
       QPointF pntDel = pnt0 - intersectionPoint;
       double distcur = sqrt(pntDel.x() * pntDel.x() + pntDel.y() * pntDel.y());
       if (distcur < dist)
       {
         dist = distcur ;
         pntOutput = intersectionPoint;
         breturn  = true;
       }
    }

  return breturn;
}
//---------------------------------------------------
void QPolygonXY::stretch(const double &scalx,const double &scaly)
{
    for (int i = 0; i < size(); ++i)
    {
      QPointF pnt = at(i);
      QPointF pnt1(pnt.x() *scalx, pnt.y() *scaly);
      replace(i,pnt1);
    }
}
//---------------------------------------------------
void QPolygonXY::stretch(QPolygonF &plg,const double &scalx,const double &scaly)
{
    for (int i = 0; i < plg.size(); ++i)
    {
      QPointF pnt = plg.at(i);
      QPointF pnt1(pnt.x() *scalx, pnt.y() *scaly);
      plg.replace(i,pnt1);
    }
}
//--------------------------------
void QPolygonXY::stretch(QVector<QPolygonF> *vctPlg,const double &scalx,const double &scaly)
{
for (int i = 0; i < vctPlg->size(); ++i)
{
    QPolygonF ln0 = vctPlg->at(i);
    stretch(ln0, scalx, scaly);
    vctPlg->replace(i,ln0);
}
}
//--------------------------------
void QPolygonXY::mirrowRect_X(QRectF &rect)
{

    QPointF topleft = rect.topLeft();

   QPointF bottomright = rect. bottomRight();

   QPointF topleft0  =topleft;

   topleft0.setY(-bottomright.y() );

   QPointF bottomright0 = bottomright;


   bottomright0.setY(-topleft.y() );

   rect.setBottomRight(bottomright0);
   rect.setTopLeft(topleft0);


}
//--------------------------------------
// поиск номера наиболее удаленной вершины
int QPolygonXY::findMostDistantVertexNumber(QPolygonF &plg,QPointF &pnt)
{
    int ireturn =-1;
    double val = 0;
    for (int i =0; i < plg.size(); ++i)
    {
        double valCur = QPointXY::calc_dist(plg.at(i), pnt);
       if (valCur > val)
       {
         val =  valCur;
         ireturn = i;
       }
    }
    return ireturn;
}
//-------------------------------------------------
//задана полилиния pln и прямоугольник Rect
//требуется найти номер самщй низ-кой вершины
// в экстенте прямоугольника
bool QPolygonXY::findUtmostLowPoint_Bounded_Rect(const QRectF &Rect, QPointF *intersectionPoint)
{
    QVector <QPointF> vctIntersect;
    findIntersectingRectPointsArray(Rect, &vctIntersect);
    int ireturn = -1;
    double val = DBL_MAX;
    for (int i =0; i < vctIntersect.size(); ++i)
    {

        if (vctIntersect.at(i).y() < val)
        {
          val = vctIntersect.at(i).y() ;
          *intersectionPoint = vctIntersect.at(i);
          ireturn = i;
        }
    }
    if (ireturn > -1)
    {
        return true;;
    }
    if (!Rect.contains(at(0)))
    {
        return false;
    }

    int i = findNumberUtmostLowPoint();
    *intersectionPoint = at(i);
    return true;
}
//-------------------------------------------------
//задана полигон  и прямоугольник Rect
//требуется найти массив точек пересечения строн прямоунгольника
// с полигоном
void QPolygonXY::findIntersectingRectPointsArray(const QRectF &Rect, QVector <QPointF> *vctIntersect)
{
    QLineF arrLn[4];
    arrLn[0] =QLineF(Rect.topLeft(), Rect.topRight());
    arrLn[1] =QLineF(Rect.topRight(), Rect.bottomRight());
    arrLn[2] =QLineF(Rect.bottomRight(), Rect.bottomLeft());
    arrLn[3] =QLineF(Rect.bottomLeft(), Rect.topLeft());

    vctIntersect->resize(0);
    for (int i =0; i < 4; ++i)
    {
        QVector <QPointF> vctIntersectCur;
        findIntersectingLinePointsArray(arrLn[i], &vctIntersectCur);
        //vctIntersect->append(vctIntersectCur);
        (*vctIntersect)+=vctIntersectCur;



    }

}
//---------------------------------------------------------
void  QPolygonXY::findIntersectingLinePointsArray(const QLineF &Line, QVector <QPointF> *vctIntersect)
{
    int len = 100;
    int lencur = len;
   vctIntersect->resize(lencur);
   int num  = 0;
   for (int i =0; i < (size()-1); ++i)
   {
    QLineF lnCur(at(i), at(i +1));

    QPointF intersectionPoint;
        if (Line.intersect(lnCur, &intersectionPoint) != QLineF::BoundedIntersection)
        {
            continue;
        }
        if ((num+1) > lencur)
        {
          lencur += len;
          vctIntersect->resize(lencur);
        }
        vctIntersect->replace(num,intersectionPoint);
        num++;

   }
   vctIntersect->resize(num);
}
//-----------------------------------------
bool   QPolygonXY::isIntersected(const QRectF &rect)
{
    QVector <QPointF> vctIntersect;
    findIntersectingRectPointsArray(rect, &vctIntersect);
    if (vctIntersect.size() > 0)
    {
    return true;
    }
    else {
    return false;
    }

}
//--------------------------------------------------
// создание точечной 3D модели осесиметричного относительно оси OX полигона
// INPUT:
// plg - полигон
// lenghtStep- шаг по оси OX(метры)
// angleStep - шаг по углу, рад
// bNoseOnly - создается единственная самая правая точка (или 2 точки)
// OUTPUT:
// buffVectPoints - вектор координат точек
// длина вектора равна 3 * (к-во точек)
void QPolygonXY::createPoints3DModel(const QPolygonF &plg
      ,const  double &lenghtStep,const  double  &angleStep
      , bool &bNoseOnly,QVector <double> &buffVectPoints)
{
  QPolygonXY plgXY = createPlgXY(plg);
  plgXY.calcBoundBox();

  const int I = (bNoseOnly)? 1:((int)((plgXY.Box[2] - plgXY.Box[0])/lenghtStep));
  const int J = 2. * M_PI / angleStep;

  int numCurPoint = 0.;
  int numPointsReserved = 500;
  buffVectPoints.resize(3 * numPointsReserved);
  for (int i =0; i < I; ++i )
  {
      double xcur = plgXY.Box[2] - ((double)i) * lenghtStep;
      QLineF ln(QPointF(xcur, 0.),QPointF(xcur, 1000000.));
      QVector <QPointF> vctIntersect;
      plgXY.findIntersectingLinePointsArray(ln, &vctIntersect);
      double r = 0.;
      for (int k = 0; k < vctIntersect.size(); ++k)
      {
          if (vctIntersect.at(k).y() > r)
          {
              r = vctIntersect.at(k).y() ;
          }
      }
      for (int j =0 ; j < J; ++j)
      {
         double ficur = ((double)j) * angleStep;
         buffVectPoints.replace( 3 * numCurPoint, xcur);
         buffVectPoints.replace( 3 * numCurPoint +1, r * cos(ficur));
         buffVectPoints.replace( 3 * numCurPoint +1, r * sin(ficur));
         numCurPoint++;
         if (numCurPoint > numPointsReserved )
         {
             numPointsReserved += 500;
             buffVectPoints.resize(numPointsReserved );
         }

      }
  }

  buffVectPoints.resize( 3 * numCurPoint);

}
//--------------------------------------
void   QPolygonXY::flip(QPolygonF &plg)
{
    for (int i = 0; i < plg.size()/2; ++i)
    {
        QPointF pnt = plg.at(i);
        plg.replace(i, plg.at(plg.size() - 1-i));
        plg.replace(plg.size() - 1-i, pnt);
    }
}
//---------------------------------------------------
// Вычисление первой точки пересечения линии, заданной начальной точкой
// pntNose и вектором arrVIzd_TargSvSK0
// OUTPUT:
// pntOutput - первая точка пересечения
bool QPolygonXY::findFirstIntersectPoint(const QPolygonF &plg,const QPointF &pntNose
                                         ,const double* arrVIzd_TargSvSK0, QPointF &pntOutput)
{
    bool breturn  = false;
    // ед вектор направления     arre[2]
    double arre[2] = {0.};
    arre[0] = arrVIzd_TargSvSK0[0];
    arre[1] = arrVIzd_TargSvSK0[1];
    double norm = sqrt( arre[0] *  arre[0]  + arre[1] *  arre[1] );
    if (norm < DBL_MIN)
    {
        return false;
    }
    arre[0] = arre[0]/norm;
    arre[1] = arre[1]/norm;
    ///

    // начальная точка нос
    QPointF pnt0 = pntNose;
    pnt0 += QPointF(-arre[0] * 1.E9, -arre[1] * 1.E9);

    //
    QPointF pnt1 = pntNose;
    pnt1 += QPointF(arre[0] * 1.E9, arre[1] * 1.E9);

    const QLineF lnDirect(pnt1,pnt0);
    ///
    double dist = DBL_MAX;
    for (int i =0; i < plg.size()-1; ++i)
    {
       const QLineF lnCur(plg.at(i), plg.at(i +1));
       QPointF intersectionPoint;
       switch(lnCur.intersect(lnDirect, &intersectionPoint))
       {
       case QLineF::UnboundedIntersection:
       case QLineF::NoIntersection:
           continue;
           break;
       default:
           break;
       }
       QPointF pntDel = pnt0 - intersectionPoint;
       double distcur = sqrt(pntDel.x() * pntDel.x() + pntDel.y() * pntDel.y());
       if (distcur < dist)
       {
         dist = distcur ;
         pntOutput = intersectionPoint;
         breturn  = true;
       }
    }

  return breturn;
}
//---------------------------------------------------
