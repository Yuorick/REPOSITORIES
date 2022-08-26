#ifndef SubWaterAnt_H
#define SubWaterAnt_H
#include <QVector>
#include <QPointF>
class QReceiver;
class TComp;


class QSubWaterAnt
{
public:
    QSubWaterAnt();

    QSubWaterAnt  (const QSubWaterAnt &R) ;
    // оператор присваивания
    QSubWaterAnt &operator=(const QSubWaterAnt  &R);
    // парам констр
    QSubWaterAnt( const QVector<QReceiver> vctSubWaterAnt, const QVector<QPointF> vctPntXY);

    QSubWaterAnt(const int QuantReceivers, const double a_Receiver
                                      , const double radius);





    virtual void fillReceiversLocation(const double radius);

    QVector<QReceiver> mvctReceiver;
    QVector<QPointF> mvctPntXY;


    int mQuant;
    double ma;


    TComp calcDiagrB(const double lamb,const double q, const double e);



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

#endif // SubWaterAnt_H
