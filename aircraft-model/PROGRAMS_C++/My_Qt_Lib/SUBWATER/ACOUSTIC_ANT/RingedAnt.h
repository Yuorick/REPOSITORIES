#ifndef RINGEDANT_H
#define RINGEDANT_H
#include <QVector>
#include <QPointF>
class QReceiver;
class TComp;


class QRingedAnt
{
public:
    QRingedAnt();

    QRingedAnt  (const QRingedAnt &R) ;
    // оператор присваивания
    QRingedAnt &operator=(const QRingedAnt  &R);
    // парам констр
    QRingedAnt( const QVector<QReceiver> vctRingedAnt, const QVector<QPointF> vctPntXY);

    QRingedAnt(const int QuantReceivers, const double a_Receiver
                                      , const double radius);

    QRingedAnt(const int QuantReceivers, const double a_Receiver
              , const double radius, const double scatter_a
                           , const double scatter_ang, const double scatter_radius ) ;

    QVector<QReceiver> mvctReceiver;
    QVector<QPointF> mvctPntXY;

    double mr;
    int mQuant;
    double ma;
    double mscatter_a;
    double mscatter_ang;
    double mscatter_radius;

    TComp calcDiagrB(const double lamb,const double q, const double e);

    void calcAngs(const double lamb,TComp *arrQ, double *pq, double *pe);

    static double calcF(double c1, double c2, double lam1, double lam2, double  al0);

    static double calcdF(double c1, double c2, double lam1, double lam2, double  al0);

    void calcVectMeasures(const double lamb,double valSignalAmp, double q, double e, TComp *arrQ );

    TComp calcSumNormalizedDiagr(const double lamb,const double q0, const double e0
                                    ,const double q, const double e);

    void calcMultiplierArr(const double lamb, const double q0, const double e0
                                       ,TComp *cmparrMultipliers);

    TComp calcIzm( const double lamb,TComp *arrQ,const double q0, const double e0 );

    void calcAmpRaznMeth_e(const double lamb,TComp *arrQ, double q0, double e0, double *pe);

    void collectData_for_e_DiagrGraph(const double lamb,const double q0
                      ,const double e0, const int QUantRows,const double step,double *arrBuff);

    void collectData_forDiagrGraph_e(double lamb,double q0,double e0 , int QUantRows
                                        ,double  step,double * arrBuff,double * arrBuffPolar);

    void calcAmpRaznMeth_q(const double lamb,TComp *arrQ, double q0, double e0, double *pq);
};

#endif // RINGEDANT_H
