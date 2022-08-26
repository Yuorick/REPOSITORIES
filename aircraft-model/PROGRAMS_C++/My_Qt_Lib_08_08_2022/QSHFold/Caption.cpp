#include "Caption.h"
#include <QPainter>
#include "PolygonXY.h"

QCaption::QCaption()
{
     mRect = QRectF();
     mLine = QLineF();

}
//---------------------------------------------------------------------------

// конструктор копирования
  QCaption:: QCaption(const  QCaption&R)
{
      mRect = R.mRect;
      mLine = R.mLine;

}
  //-------------------------------------------------------------------

// оператор присваивания
QCaption &QCaption::operator=( const QCaption &R)
{
    mRect = R.mRect;
    mLine = R.mLine;

  return *this ;
}

//-------------------------------------------------------------------

// парам конструктор 1
QCaption:: QCaption( const QRectF Rect, const QLineF Line)
{
    mRect = Rect;
    mLine = Line;
}

//-------------------------------------------------------------------

// парам конструктор 2
QCaption:: QCaption( const QRectF Rect)
{
    mRect = Rect;
    // центр прямоугольника

    QPointF pntCentre = Rect.center();
    mLine.setP1(pntCentre);
    mLine.setP2(pntCentre);

}

//--------------------------------------

void QCaption::translate(qreal dx, qreal dy)
{
   mRect.translate( dx, dy);
   mLine.setP1(mRect.center());


}
//-----------------------------------------
void QCaption::translate(const QPointF &offset)
{
    mRect.translate( offset);

    mLine.setP1(mRect.center());

}
//-----------------------------------------
void QCaption::draw(QPainter &painter,const QString &str, bool bRect)
{
  painter.drawLine(mLine);
  painter.setBrush(QBrush(Qt::white));
  qreal xRadius=4, yRadius =4;
  if(bRect)
  {
  painter.drawRoundedRect(mRect, xRadius, yRadius);
  }
  painter.drawText(mRect, Qt::AlignCenter|Qt::TextDontClip,str);

}
//-----------------------------------------------
void QCaption::mirrow_X()
{
    QPolygonXY::mirrowRect_X(mRect);
    double scalx =1.,scaly = -1.;
    QLineXY::stretch(mLine,scalx,scaly);
}

