#ifndef TRUEMEASPARAMS_H
#define TRUEMEASPARAMS_H


class QTrueMeasParams
{
public:
    QTrueMeasParams();

    QTrueMeasParams (const  QTrueMeasParams &R);

    QTrueMeasParams  &operator=( const QTrueMeasParams  &R);

    QTrueMeasParams (const double  *arrSVessWave
                     ,const double  *arrMuWave,const double  *arrMuWaveZv,
                       const double  *arrSVess
                     ,const double  *arrMu,const double  *arrMuZv,
                       const double  Tzapr,const double  Totv
                       ,const double q, const double e);

    // момент запросного сигнала
    double mTzapr;
    // вектор положения Vess в ГСК
    double marrSVessWave[3];
    // вектор углов эйлера корабля на момент mTzapr
    double marrMuWave[3];
    // вектор измеренных углов Эйлера на момент mTzapr
    double marrMuWaveZv[3];

    // момент ответного сигнала
    double mTotv;
    // вектор положения центра масс корабля на момент mTotv в ГСК
    double marrSVess[3];
    // вектор углов эйлера корабля на момент mTotv
    double marrMu[3];    
    // вектор измеренных углов Эйлера на момент mTotv
    double marrMuZv[3];

    // угловые измерения
    double mq;
    double me;
};

#endif // TRUEMEASPARAMS_H
