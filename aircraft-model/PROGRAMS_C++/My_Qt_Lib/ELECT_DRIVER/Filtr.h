#ifndef FILTR_H
#define FILTR_H


#define QUANT_VARS 5
#define QUANT_OBS  3

class QMeasure;
class QDriverMeasure ;

class QFiltr
{
public:
    QFiltr();
    // длина фазового вектора
    int mQntVars;
    // число наблюдаемых переменных
    int mQntObserv;
    double *marrCurEst;
    double *marrCurK;
    double mTCur;
    double mh;
    bool mbInit;

   // матрица наблюдений
    double *marrC;

    virtual ~QFiltr();

    QFiltr (const  QFiltr &R);

    QFiltr  &operator=( const QFiltr  &R);

    QFiltr ( const int QntVars, const int QntObserv,const double *parrCurEst
     , const double *parrCurK, const double TCur, const double h
     , const bool bInit, const double *parrC);

    QFiltr ( const int QntVars, const int QntObserv,const double *parrCurEst
      , const double TCur, const double h  );
    //--------------------------------------------------------

    void init_(const double valTcur,const double valh, const double*arrPhVect);

   // virtual void processMeasure(QMeasure MEasureCur,const double VAlTcur
              //                  ,double *arr_dF_po_dx,double *arr_dF_po_dW,double *arrKww);

    virtual void processMeasure_Type2( double *ARrYZv, double *ARrK,const double VAlTcur
                                      ,double *arrL,double *arrKww);

    void processMeasure(QDriverMeasure MEasureCur,const double VAlTcur
                                ,double *arrA,double *arrBU,double *arrKww);

};

#endif // FILTR_H
