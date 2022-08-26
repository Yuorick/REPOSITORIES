#ifndef BIGMEASURE_H
#define BIGMEASURE_H


class QBigMeasure
{
public:
   // ~QBigMeasure() ;

    QBigMeasure();

    QBigMeasure (const  QBigMeasure &R);

    QBigMeasure  &operator=( const QBigMeasure  &R);

    QBigMeasure (const double  *arrSVessWaveZv,const double  *arrMuWaveZv,
                      const double  *arrSVessZv,const double  *arrMuZv,
                      const double  TzaprZv,const double  TotvZv,
                      const double  Tobr,const int  DimMeas
                      ,const double qzv, const double ezv
                      ,const double Sig_t,const double Sig_q,const double Sig_e);


    // момент запросного сигнала
    double mTzaprZv;
    // вектор положения центра масс корабля на момент mTzaprZv
    double marrSVessWaveZv[3];
    // вектор углов эйлера корабля на момент mTzaprZv
    double marrMuWaveZv[3];

    // момент ответного сигнала
    double mTotvZv;
    // вектор положения центра масс корабля на момент mTotvZv
    double marrSVessZv[3];
    // вектор углов эйлера корабля на момент mTotvZv
    double marrMuZv[3];

    //время обработки на маяке
    double mTobr;

    // размерность измерения антенны
    int mDimMeas;

    // угловые измерения
    double mqzv;
    double mezv;
    // СКЗ измерения времени
    double mSig_t;
    // СКЗ измерения КУ
    double mSig_q;
    // СКЗ измерения УМ
    double mSig_e;

    static bool isEqual(QBigMeasure &meas0, QBigMeasure &meas1, const double VAlTolerDist
                  ,const double VAlTolerSinsAngs,  const double VAlTolerT, const double VAlTolerMeasAngs );

    void fill_info_row( double *arrP_Gps, double *arrRow);
};

#endif // BigMEASURE_H
