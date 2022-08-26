#ifndef PNTCLOUD_H
#define PNTCLOUD_H
#include <QVector>
#include <QVector3D>
#include "Table_1D.h"
#include "PosSolver.h"


class TURPointXY ;
class QBigMeasure;
class TTable_1D;
enum TYPE_OF_OUTPUT_FILE
{
    SHP,         //
    CSV
};
class QPntCloud
{
public:
    QPntCloud();

    QPntCloud (const QPntCloud &R);

    QPntCloud &operator=(const QPntCloud  &R);

    QPntCloud(  const QVector<QBigMeasure> vctBigMesure,const TTable_1D tblEstPrfl
                             ,const double* arrSBeacon_XYZ,const double* arrAntPosParams, TYPE_OF_SOLVER_TASK TYPE);

    QVector<QBigMeasure> mvctBigMesure;

     TTable_1D mtblEstPrfl;

     double marrSBeacon_XYZ[3];

     double marrAntPosParams[6];

     TYPE_OF_SOLVER_TASK mTYPE_OF_SOLVER_TASK;



     void createCloud_3D(QVector<QVector3D> *pvctCloud3D, QVector<int> *pvctNumGood, double *parrBuff);

     QVector3D createMeanVct3D(const QVector<QVector3D> &vctCloud3D);

     void createDispMtrx3D(const QVector<QVector3D> &vctCloud3D, const QVector3D &vctMean
                                       ,double *arrDispMtrx);

     void createClouds_XY_YZ_ZX(QVector<QPointF> *pvctCloud_XY
                     ,QVector<QPointF> *pvctCloud_YZ,QVector<QPointF> *pvctCloud_ZY,double *arrMean
                     , double *arrDisp,QVector<int> *pvctNumGood,double * parrBuff);

     static void createPictFiles(wchar_t *wchOutPutFold,const QVector<QBigMeasure> vctBigMeasure
                                 ,const TTable_1D tblEstPrfl,const double* arrSBeacon_XYZ,const double* arrAntPosParams
                               , TYPE_OF_SOLVER_TASK TYPE,double *arrMean, double *arrDisp,QVector<int> *pvctNumGood
                                                        , TYPE_OF_OUTPUT_FILE Type_of_Output_File);

     static void transf_QPointArray_in_TURPointXYArray(TURPointXY *pPntArr,const QVector<QPointF> vctCloud);
};

#endif // PNTCLOUD_H
