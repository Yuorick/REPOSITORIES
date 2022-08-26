#ifndef POSSOLVER_H
#define POSSOLVER_H
#include <QVector>

#include "Table_1D.h"

enum TYPE_OF_SOLVER_TASK
{
    LBL,         //
    USBL_2D,   //
    USBL_3D //
};

class QBigMeasure;
class TTable_1D;
class QPosSolver
{
public:
    virtual~QPosSolver();

    QPosSolver();

    // конструктор копирования
    QPosSolver  (const QPosSolver &R) ;
     // оператор присваивания
    QPosSolver &operator=(const QPosSolver  &R);
    // парам констр
    QPosSolver(const QBigMeasure *parrBigMeasures
               , const int  QntMeas    ,const TTable_1D tblEstPrfl
               , const double Deepth
               , const int DimX, const int DimY, const double Toler);

    QPosSolver(QVector<QBigMeasure> VectBigMeasures
               ,const TTable_1D tblEstPrfl, const double Deepth
            , const int DimX, const int DimY, const double Toler );



    TTable_1D mtblEstPrfl;

    // вектор измерений
    //QBigMeasure *mparrBigMeasures;

    //int mQntMeas;
    QVector<QBigMeasure>mVectBigMeasures;

    // вектор параллакса PosSolver априорный
   // double marrGPSAprioriPosParams[3];

    // размерность вектора переменных
    int mDimX;

    // размерность вектора измерений
    int mDimY;

    // глубина маяка (координата Z маяка)
    double mDeepth;

    // точность решения
    double mToler;



    // максимальное к-во итераций
    int mQuantIterMax;


     bool estimateParams(  double *arrX0,double *arrXRez, double &valNeviaz0);


    double calcSumNeviazka(double *arrXCur, int &quantGoodZamers);

    bool calcPartialNeviazkaSquare(double *arrX,const int NUmBigMeasure, double &valNeviaz2);

    bool doOneIteration(double *arrX0
               ,double *arrDel, double *arrMtrx_dFGr_po_dX_Inv);

    bool calc_FGr_and_dFGr_po_dX(double *arrX0
                   ,double * parrF,double * parrMtrx_dF_po_dX, double &valSumVeviaz2);

    bool  calc_arrfi_and_arrdfi_po_dx(double *arrX
             ,const int NUmBigMeasure,double *arrNeviaz,double *arrfi,double *arrMtrx_dfi_po_dX);

    virtual void  collectGrlsMeasure(const QBigMeasure &Meas
                                                 ,double *arrMeasure,double *arrD);


    virtual bool calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                           , double *arrNeviaz , double *dArrNeviaz_po_dX);


    bool  calc_arrfi(double *arrX, const int NUmBigMeasure,double *arrNeviaz,double *arrfi);

    void  calcNeviazArray_(double * arrX0, double *parrBuffNeviaz,double *arrSKZ);

   // void  calcNeviazArray(double * arrX0,double *parrBuffNeviaz,double *arrSKZ);

    void repaireZamerArray(double *arrX0);

    void  calcNeviazArray(double * arrX0,QVector <double> *vectBuffNeviaz,double *arrSKZ);

    static void unloadQvectorToDoubleArray(QVector <double> *vectInp,double *arrOut);

    void repaireZamerArray_(double *arrX0);

    void repaireZamerArray_(double *arrX0, const double VAlSigmaTreshold
                            ,const int NUmVar, QVector<int> *pvctGoodZameresNum);

    virtual bool calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                                                , double *arrNeviaz );

    void fncRoughCleaningMeasures(double *arrX0,double valTreshold, const int NUmPeremennaya,QVector<int>*pvctGoodZamersNum);

    double calcArrNeviaz_Mean_and_Disp(double *arrXCur, int &quantGoodZamers
                                                 , double *arrMean, double *arrDisp, double *arrNeviazSquare);

    bool calc_arrFGr_and_Mtrx_dFGr_po_dX_Inv(double *arrX0,double *arrFGr,double *arrMtrx_dFGr_po_dX_Inv);

    static void createDblVect(const double *arr,const int len,QVector<double> &vct);

    static void createDblArrayFromDblVector(double *arr,QVector<double> vct);

    virtual void calcArrPenalty(const double *arrX0,const double valPenalt, double *arrPenalt);

    void cleanZamersInAccordanceWithAntHeight(double *arrX0, const double VAlMaxHeightTreshold);

    bool calc_darrfi_po_dSvess(double *arrX, const int NUmBigMeasure
                               , const double del,double *parr_dfi_po_dSvess);

    bool  calc_arrfi_and_arrdfi_po_dt(double *arrX
                       ,const int NUmBigMeasure, const double valPrirashen,double *arrfi,double *arrdfi_po_dt);

    bool  calc_arrfi_and_arrdfi_po_dPsiTetta(double *arrX,const int NUmBigMeasure,double *arrPrirashen
                                      ,double *arrdfi_po_dPsiTetta);

    bool  calc_arrfi_and_arrdfi_po_dQ(double *arrX
                                      ,const int NUmBigMeasure, const double valPrirashen,double *arrdfi_po_dQ);

    bool  calc_arrfi_and_arrdfi_po_dSgps(double *arrX
                           ,const int NUmBigMeasure,double *dSgps);

    bool  calc_arrfi_and_arrdfi_po_dAngSins(double *arrX
                           ,const int NUmBigMeasure,double *arrdfi_po_dAngSins);

    bool  calc_dFGr_po_dSins_(double* arrX0,double*  parrMtrx_dFGr_po_dSins);

    bool  calc_dFGr_po_dProfile(double* arrAngX0,double*  parrMtrx_dFGrA_po_dProfile);

};

#endif // POSSOLVER_H
