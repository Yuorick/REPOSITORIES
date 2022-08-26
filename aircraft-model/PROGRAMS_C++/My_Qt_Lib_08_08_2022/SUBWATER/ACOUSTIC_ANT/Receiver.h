#ifndef RECEIVER_H
#define RECEIVER_H


class QReceiver
{
public:
    QReceiver();
    QReceiver  (const QReceiver &R) ;
    // оператор присваивания
    QReceiver &operator=(const QReceiver  &R);
    // парам констр
    QReceiver(  const double a) ;

    // диаметр модуля
    double ma;

    double fncBessDiagr(const double tetta, const double lamb);

    double calcSignalAmp(const double dist,const double transmitterWaveAmp, const double lamb);

   // double clcAmpBessDiagr(const double tetta,const double dist
     //                                 ,const double transmitterWaveAmp, const double lamb);

    double calcSigmaNoise();

    void collectData_forDiagrGraph(double lamb, int QUantRows
                                              , double step, double *arrBuff, double *arrBuffPolar);
};

double fncBessel1(const double x);

double fnc_dBessel1_po_dx(const double x);

double fncBx(const double x);

double fnc_dBx_po_dx(const double x);



#endif // RECEIVER_H
