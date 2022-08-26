#ifndef CAPTION_H
#define CAPTION_H

#include <QRectF>
#include <QLineF>

class QRectF;
class QPainter;

class QCaption
{
public:
    QCaption();

    QCaption(const  QCaption&R);

    QCaption &operator=( const QCaption &R);

    QCaption( const QRectF Rect, const QLineF Line);

    QCaption( const QRectF Rect);

    QRectF mRect;
    QLineF mLine;

    void translate(qreal dx, qreal dy);

    void translate(const QPointF &offset);

    void draw(QPainter &painter,const QString &str, bool bRect);

    void mirrow_X();
};

#endif // CAPTION_H
