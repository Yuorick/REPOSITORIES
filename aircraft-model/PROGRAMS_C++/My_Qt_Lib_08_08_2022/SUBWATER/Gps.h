#ifndef GPS_H
#define GPS_H


class QGps
{
public:
    QGps();
    // конструктор копирования
    QGps  (const QGps &R) ;
     // оператор присваивания
    QGps &operator=(const QGps  &R);
    // парам констр
    QGps(  const double SigXY,const double* arrPos) ;

    void imitateMeasure(double *arrTruePosition_GSK, double *arrPosition_GSK_Zv);

    double mSigXY;
    double marrPos[3];

};

#endif // GPS_H
