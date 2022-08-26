//---------------------------------------------------------------------------

#ifndef HIDRORLS_H
#define HIDRORLS_H

//---------------------------------------------------------------------------
//enum enumTypeOfRLS{ LBL,USBL_2D,USBL_3D};

#define NODATA -99999

class QTrueMeasParams;
class QBigMeasure;
class QHidroRLS
{
public:

// точность определения времени (шаг дискретизации)
    double mSigT;

 // Ширина диаграммы (рад)
    double mDiagWidth ;


	// конструктор по умолчанию
    QHidroRLS () ;
	// конструктор копирования
    QHidroRLS  (const QHidroRLS  &R) ;
	// оператор присваивания
    QHidroRLS  &operator=(const QHidroRLS   &R2) ;

	// парам конструктор1
    QHidroRLS ( const double SigR, const double DiagWidth);


    void imitateMesure(const QTrueMeasParams trueMeasParams, QBigMeasure &Meas);

    void imitateMesure(const QTrueMeasParams trueMeasParams, double *valTzaprZv
                                  , double *valTotvZv, double *valSig_t, double *val_qzv
                                  , double *val_ezv, double *valSig_q, double *valSig_e);

    virtual void calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                  , double &Sig_q, double &Sig_e);



    virtual int createInputDataReport(wchar_t*FileName, const bool bHeader);

    double calcDeltaT(const QTrueMeasParams trueMeasParams);

    double calcSig_t(const QTrueMeasParams trueMeasParams);



}  ;
#endif
