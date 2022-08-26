#ifndef MYTHREAD_H
#define MYTHREAD_H

#include <QThread>
#include "ImitMod.h"
#include "Table_1D.h"
#include <QString>
#include "PosSolver.h"
#include "BigMeasure.h"
#include "Solver2D_Angs.h"
#include "PntCloud.h"

enum MYTHREAD_RESULT
{
    MYTHREAD_RESULT_OK,         // вычисления завершены успешно (достигли требуемой точности)
    MYTHREAD_RESULT_MAX_STEP,   // вычисления провалены (достигли максимально допустимого кол-ва итераций)
    MYTHREAD_RESULT_USER_CANCEL, // вычисления провалены (прервано пользователем)
    MYTHREAD_RESULT_ERROR, // вычисления провалены (ошибка в исходных данных)
    MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_LBL, // вычисления провалены после чистки по 3 сигма по задержке врем мало замеров
    MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_USBL // вычисления провалены после чистки по 3 сигма по углам мало замеров
};

enum ALGOR_TYPE
{
    PEREBOR,
    NUTON,
};
class QBigMeasure;
class MyThread : public QThread
{
    Q_OBJECT
public:
    explicit MyThread(QObject *parent = 0);   

    void userBreak() { userBreakRequested = true; } // установка флага завершения работы потока

    void setParams(QVector <QBigMeasure> vctMeas, TYPE_OF_SOLVER_TASK TypeOfSolverTask
                   ,QVector <double> vctX0,TTable_1D tblEstPrfl, double *arrGPSAprioriPosParams
                  ,double HBeaconEst,int quantIter_LBL, int quantIter_USBL
                    ,double Gape,double AngStep, ALGOR_TYPE nALGOR_TYPE,double LBLcoeff,double USBLcoeff
                    , TYPE_OF_OUTPUT_FILE Type_of_Output_File, QString qstrOutPutFold);

    MYTHREAD_RESULT estimateAngs(QSolver2D_Angs *pPosSolver,  double *arrX0
         , int numPart,const int MaxQuantIter,double *arrXRez
         , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv);

    MYTHREAD_RESULT estimateParams(QPosSolver *pPosSolver,  double *arrX0
                                   , int numPart,const int MaxQuantIter, double coeff,double *arrXRez
                                   , double &valNeviaz0, double *arrMtrx_dFGr_po_dX_Inv);

    // 1
    QVector <QBigMeasure> mvctMeas;

    //2
    TYPE_OF_SOLVER_TASK mTypeOfSolverTask;

    //3
    QVector <double> mvctX0;

    //4
    TTable_1D mtblEstPrfl;

    //5  вектор параметров позиционироваия GPS априорный
     double marrGPSAprioriPosParams[3];

     // 6
     double mHBeaconEst;

     // 7 к-во итераций процесса LBL
     int mMaxQuantIter_LBL;

     // 8 к-во итераций процесса LBL
     int mMaxQuantIter_USBL;

     // 8 ворота для перебора по углу
     double mGape;

     // 9 шаг перебора по углу
     double mAngStep;

     ALGOR_TYPE mALGOR_TYPE;

     double mLBLcoeff;

     double mUSBLcoeff;

     TYPE_OF_OUTPUT_FILE mType_of_Output_File;

     QString mqstrOutPutFold;



signals:

    void progress(int numIter, QVector<double>vctNevSquare, QVector<double> vctMean
                  , QVector<double> vctDisp, int numpart ,bool bEndPart, QVector<double>vctX0);

    void calcfinished(int result);

    void kleaningEnded(int numpart, QVector<int> vctGoodZamersNum
                       , QVector<double> vctMean,QVector<double> vctDisp);

public slots:

private:
    void run();

private:
    bool userBreakRequested;


};

#endif // MYTHREAD_H
