#include "Segment.h"
#include "UrPointXY.h"
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"
#include "LnSgm.h"

QSegment::QSegment():TURPolyLine(1,2)
{
 Parts[0] = 0;
 Points[0]=Points[1]= TURPointXY(0.,0.);
}
// парм консируктор 1
QSegment::QSegment(const  TURPolyLine pln):TURPolyLine(1,2)
{
    Parts[0] = 0;
    Points[0] = pln.Points[0];
    Points[1] = pln.Points[1];
}

// парм консируктор 2
//QSegment::QSegment(const TURPointXY pnt0,const TURPointXY pnt1):TURPolyLine(1,2)
//{
//    Points[0] = pnt0;
//    Points[1] = pnt1;
//}

QSegment::QSegment(const TURPointXY  pnt1, const TURPointXY  pnt2):TURPolyLine( pnt1, pnt2)
{
}


//-------------------------------------------------------------
// нахождение незатенны участков сегмента segmBase
//Заданы 2 сегмента  segmBase и segmBarrier
//задано направление (вектор)pntV
//пусть свет падает по направлениею pntV
// сегмент segmBarrier создает за собой тень
//требуется найти освещенные участки сегмента  segmBase
// освещенных участков может быть максимум 2
// функция возвращает к-во освещенных участков
// arrSegm[2] -массив освещенных сегментов в сиситеме координат сегмента segmBase
// ось OX направлена от точки segmBase.Points[0] к точке segmBase.Points[1]
// начло координат в точке segmBase
void QSegment::calcOpen_and_Shadow_Segments( QSegment segmBase,  QSegment segmBarrier,  double *arrV
              , QSegment *arrOpenSegm, int *piQuantOpenSegm, QSegment *pShadowSegm, int *piQuantShadowSegm )
{
    const double valLength = segmBase.calcLeng();
    const TURPointXY  pntZero(0.,0.), pntEnd(valLength, 0.);
    if (segmBarrier.calcLeng() <= 0.00000001)
    {
      arrOpenSegm[0] = QSegment(pntZero, pntEnd);
      *piQuantOpenSegm = 1;
      *piQuantShadowSegm =0;
      return;
    }

    if ( NormVect3(arrV) < 0.0000001)
    {
      arrOpenSegm[0] = QSegment(pntZero, pntEnd);
      *piQuantOpenSegm = 1;
      *piQuantShadowSegm =0;
      return;
    }
    // 1. вычисление угла наклона сегмента к оси OX
    double valAngSlope = segmBase.calcSlopeAng();
    ///

    // 2. перевод сегмента segmBarrier в сиситему координат сегмента segmBase
    // ось OX направлена от точки 0 к точке 1, начало коорлиннат в точке 0
    TURPolyLine pln0 = segmBarrier.LinTransform(-valAngSlope , segmBase.Points[0],1. ) ;
    QSegment segmBarrier0 (pln0);
    ///

    // 3.  перевод вектора arrV в сиситему координат сегмента segmBase
    // ось OX направлена от точки 0 к точке 1, начало коорлиннат в точке 0
    TURPointXY pnt(arrV[0], arrV[1]), pntXero(0.,0.);
    TURPointXY pnt0 = pnt.LinTransform(-valAngSlope  , pntXero,1. );
    if(fabs(pnt0.Y) < 2 * DBL_MIN)
    {  // весь сегмент в тени
        *pShadowSegm = QSegment(pntZero, pntEnd);
        *piQuantOpenSegm = 0;
        *piQuantShadowSegm =1;
        return;
    }
    double arrV0[2] = {0.};
    arrV0[0] = pnt0.X;
    arrV0[1] = pnt0.Y;
    ///

    //4. нахождение точек пересечения arrx[0] и arrx[1]линий segmBarrier0.Points[0] + t* arrV0
    // и segmBarrier0.Points[1] + t* arrV0 с осью OX системы координат сегмента segmBase
    // (то есть линией сегмента)

    double arrx[2] = {0.};
    double arrt[2] = {0.};
    arrt[0] = -segmBarrier0.Points[0].Y/ arrV0[1];
    arrx[0] = segmBarrier0.Points[0].X + arrV0[0] * arrt[0] ;
    arrt[1] = -segmBarrier0.Points[1].Y/ arrV0[1];
    arrx[1] = segmBarrier0.Points[1].X + arrV0[0]* arrt[1] ;

    ///


    //сортировка точек    arrx[0]и arrx[1] по возрастанию:
    if (arrx[1] < arrx[0])
    {
      double t = arrx[0];
      arrx[0] = arrx[1];
      arrx[1] = t;
      t = arrt[0];
      arrt[0] = arrt[1];
      arrt[1] = t;
    }
    ///


     const double EPS = 0.000001;

     if (fabs (arrx[1] - arrx[0]) < 10. * DBL_MIN )
     {  // тень очень мала, весь сегменьт открыт
         arrOpenSegm[0] = QSegment(pntZero, pntEnd);
         *piQuantOpenSegm = 1;
         *piQuantShadowSegm =0;
         return;
     }

    // случай, когда одна из точек segmBarrier0 или сразу обе расположена за прямой сегмента
    if(arrt[0] < -0.000000001)// точка 0 лежит за сегментом и не загораживает прямую
    {
        if(arrt[1] < EPS)
        { // точка 1 расположена за прямой сегмента-> барьер не загораживает прямую

            arrOpenSegm[0].Points[0] = TURPointXY(0.,0.);
            arrOpenSegm[0].Points[1] = TURPointXY(valLength ,0.);
            *piQuantOpenSegm = 1;            
            return ;
        }
        else
        { // точка 1  расположена до прямой сегмента

           double valk = ( segmBarrier0.Points[1].Y - segmBarrier0.Points[0].Y) / ( segmBarrier0.Points[1].X - segmBarrier0.Points[0].X);
           arrx[0] = segmBarrier0.Points[0].X -  segmBarrier0.Points[0].Y / valk;

        }
    }
    else
    {
      // точка 0 лежит перед прямой сегмента
        if(arrt[1] > EPS)
        {
            // точка 1 лежит за прямой сегмента
            double valk = ( segmBarrier0.Points[1].Y - segmBarrier0.Points[0].Y) / ( segmBarrier0.Points[1].X - segmBarrier0.Points[0].X);
            arrx[1] = segmBarrier0.Points[0].X -  segmBarrier0.Points[0].Y / valk;
        }
    }

   ///

    // теперь надо найти пересечение отрезка [0;valLength] и отрезка  [arrx[0];arrx[1]]
    QLnSgm LnSgmBase(0., valLength);
    QLnSgm LnSgmShadow(arrx[0],arrx[1]);
    QLnSgm LnSgmProduct, arrLnSgm[2];

    *piQuantShadowSegm =0;
    if ( QLnSgm::productOfSgms( LnSgmBase, LnSgmShadow, &LnSgmProduct))
    {
        *piQuantShadowSegm =1;
        (*pShadowSegm).Points[0] = TURPointXY(LnSgmProduct.ma,0.);
        (*pShadowSegm).Points[0] = TURPointXY(LnSgmProduct.mb,0.);

    }
    *piQuantOpenSegm =  QLnSgm::SegmMinusSegm(  LnSgmBase, LnSgmShadow, arrLnSgm);
    for (int i = 0; i < (*piQuantOpenSegm); ++i )
    {
       arrOpenSegm[i].Points[0] = TURPointXY(arrLnSgm[i].ma,0.);
       (*pShadowSegm).Points[1] = TURPointXY(arrLnSgm[i].mb,0.);
    }

}

//================================================
bool  QSegment::IsVertical()
{
  return(fabs(Points[0].X -  Points[1].X )< 2.* DBL_MIN);
}
//================================================

// нахожление тангенса угла наклона линии сегмента
bool QSegment::calcTang(double *pvalTang)
{
  if (IsVertical())
  {
    return false;
  }
   *pvalTang =  (Points[0].Y -  Points[1].Y )/(Points[0].X -  Points[1].X );
   return true;
}//================================================

// нахожление угла наклона линии сегмента
double QSegment::calcSlopeAng()
{
  if (IsVertical())
  {
    if(Points[1].Y > Points[0].Y)
    {
        return M_PI / 2.;
    }
    return -M_PI / 2.;
  }

   return atan2((Points[1].Y - Points[0].Y),(Points[1].X - Points[0].X));
}
//---------------------------------------------------------------
// нахождение значения по оси OY точки на линии сегмента с координатой по оси OX равной x
bool  QSegment::calcY(const double x, double *py)
{
  if (IsVertical())
  {
    return false;
  }
  double valTang ;
  calcTang(&valTang);
  *py =  valTang * (x -  Points[0].X) + Points[0].Y;
  return true;

}

//---------------------------------------------------------------
// нахождение  точки на линии сегмента с координатой по оси OX равной x
bool  QSegment::calcPointXY(const double x, TURPointXY *pPointXY)
{
   double y;
  if(!calcY( x, &y) )
  {
      return false;
  }
  *pPointXY = TURPointXY(x,y);
   return true;

}
//-----------------------------------------

