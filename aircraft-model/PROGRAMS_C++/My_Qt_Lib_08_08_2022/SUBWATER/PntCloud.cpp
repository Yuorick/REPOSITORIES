#include "PntCloud.h"
#include "BigMeasure.h"
#include <QPointF>
#include "SubWaterBeam.h"
#include "MatrixProccess.h"
#include "UrPointXY.h"
#include "YrWriteShapeFile.h"
#include "URPolygon.h"
#include "URPolyLine.h"
#include "YrWrite.h"


QPntCloud ::QPntCloud ()
{
    mtblEstPrfl = TTable_1D();

    memset(marrSBeacon_XYZ, 0, 3 * sizeof(double));
    memset(marrAntPosParams, 0, 6 * sizeof(double));
    mTYPE_OF_SOLVER_TASK = LBL;
}
// конструктор копирования
 QPntCloud ::QPntCloud (const QPntCloud &R)
 {
     mvctBigMesure = R.mvctBigMesure ;
     mtblEstPrfl = R.mtblEstPrfl;
     memcpy(marrSBeacon_XYZ, R.marrSBeacon_XYZ, 3 * sizeof(double));
     memcpy(marrAntPosParams, R.marrAntPosParams, 6 * sizeof(double));
     mTYPE_OF_SOLVER_TASK = R.mTYPE_OF_SOLVER_TASK;
 }


// оператор присваивания
QPntCloud &QPntCloud::operator=(const QPntCloud  &R)
{
    mvctBigMesure = R.mvctBigMesure ;
    mtblEstPrfl = R.mtblEstPrfl;
    memcpy(marrSBeacon_XYZ, R.marrSBeacon_XYZ, 3 * sizeof(double));
    memcpy(marrAntPosParams, R.marrAntPosParams, 6 * sizeof(double));
    mTYPE_OF_SOLVER_TASK = R.mTYPE_OF_SOLVER_TASK;

  return *this ;
}
// парам констр
QPntCloud :: QPntCloud(  const QVector<QBigMeasure> vctBigMesure,const TTable_1D tblEstPrfl
                         ,const double* arrSBeacon_XYZ,const double* arrAntPosParams, TYPE_OF_SOLVER_TASK TYPE)

{
    mvctBigMesure = vctBigMesure ;
    mtblEstPrfl = tblEstPrfl;
    memcpy(marrSBeacon_XYZ, arrSBeacon_XYZ, 3 * sizeof(double));
    memcpy(marrAntPosParams, arrAntPosParams, 6 * sizeof(double));
    mTYPE_OF_SOLVER_TASK = TYPE;
}
//-------------------------------------
void QPntCloud ::createCloud_3D(QVector<QVector3D> *pvctCloud3D, QVector<int> *pvctNumGood, double *parrBuff)
{

   memset(parrBuff, 0, mvctBigMesure.size());
   pvctCloud3D->resize(mvctBigMesure.size());
   pvctNumGood->resize(mvctBigMesure.size());
   QVector<bool> pvct_b(mvctBigMesure.size());
   pvct_b.fill(true);

   //#pragma omp parallel num_threads(nt) // OMP (главная директива для запуска нескольких потоков, количество потоков nt)
  #pragma omp parallel // OMP (Если не указывать количество потоков nt, то по умолчанию будет использовано максимальное количество потоков)
   { // OMP (начало блока, который выполняется в нескольких потоках
  #pragma omp for // OMP (директива для распределения итераций цикла между потоками)
    for (int i =0; i < mvctBigMesure.size(); ++i)
    {
        QBigMeasure meas = mvctBigMesure.at(i);
        double val_q = 0., val_e =0., val_t = 0.;
        parrBuff[i *10]=meas.mTotvZv;
        if(!transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, marrSBeacon_XYZ, meas.marrSVessZv, meas.marrMuZv
          ,marrAntPosParams,  &val_q, &val_e,  &val_t))
        {
            for (int j =1;j < 8;++j)
            {
              parrBuff[j] = -99999.;
            }
            pvct_b.replace(i, false);
            continue;
        }
        switch(mTYPE_OF_SOLVER_TASK)
        {
        case LBL:
         meas.mqzv =val_q;
         meas.mezv =val_e;
            break;
        case USBL_2D:
            meas.mezv =val_e;
            break;

        default:
            break;
        }


   // формирование буфера невязок
        // время измерения
        // невязка по времени
        // КУ теор
        // КУ Зв
        // нев КУ = КУзв -КУ
        // УМтеор
        // УМзв
        // нев УМ

        double val_qWave = 0., val_eWave= 0.,  val_tWave = 0.;
        if(!transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, marrSBeacon_XYZ, meas.marrSVessWaveZv, meas.marrMuWaveZv
          ,marrAntPosParams,  &val_qWave, &val_eWave,  &val_tWave))
        {
           for (int j =1;j < 10;++j)
           {
             parrBuff[j] = -99999.;
           }
        }
        else
        {
        parrBuff[i *10 +1]= (meas.mTotvZv- meas.mTzaprZv - meas.mTobr - val_t -val_tWave) * 1000.;
        parrBuff[i *10 +2]= val_q * 180./ M_PI;
        parrBuff[i *10 +3] = meas.mqzv* 180./ M_PI;
        parrBuff[i *10 +4] = (meas.mqzv - val_q)* 180./ M_PI;
        parrBuff[i *10 +5] = val_e* 180./ M_PI;
        parrBuff[i *10 +6] = meas.mezv* 180./ M_PI;
        parrBuff[i *10 +7] = (meas.mezv - val_e)* 180./ M_PI;
        parrBuff[i *10 +8] = meas.mTotvZv- meas.mTzaprZv;
        parrBuff[i *10 +9] = meas.mTobr + val_t +val_tWave;
        }
// формирование буфера невязок !


        double valTZv = 0.;
        double arrSBeaconField[3] = {0.};

        if(!transfMeasure_to_GSK(mtblEstPrfl,meas,marrAntPosParams,NULL
                                  ,  arrSBeaconField, &valTZv))
        {
            pvct_b.replace(i, false);
            continue;
        }


       MtrxMinusMatrx(arrSBeaconField, marrSBeacon_XYZ,1, 3,arrSBeaconField) ;
       QVector3D vect3D((float)arrSBeaconField[0], (float)arrSBeaconField[1],(float) arrSBeaconField[2]);
       pvctCloud3D->replace(i,vect3D);
    }
   }

    int icur = 0 ;
    for(int i = (mvctBigMesure.size() - 1); i >=0;--i)
    {
        if(!pvct_b.at(i))
        {
          pvctCloud3D->remove(i) ;
          continue;
        }
       pvctNumGood->replace(icur,i);
       icur++;
    }

    pvctNumGood->resize(pvctCloud3D->size());
}
//-------------------------------------------------------
QVector3D QPntCloud ::createMeanVct3D(const QVector<QVector3D> &vctCloud3D)
{
  QVector3D vctMean(0.,0.,0.);
  for (int i = 0; i < vctCloud3D.size(); ++i)
  {
     vctMean+= vctCloud3D.at(i);
  }
  vctMean *=1. / ((float)vctCloud3D.size());
  return vctMean;
}
//----------------------------------------------------------
void QPntCloud ::createDispMtrx3D(const QVector<QVector3D> &vctCloud3D, const QVector3D &vctMean
                                  ,double *arrDispMtrx)
{
    memset (arrDispMtrx, 0, 9 * sizeof(double));
    double arrMtrx[9] = {0.};
    QVector3D vctCur;
    double arrCur[3]= {0.},arrCur1[3]= {0.};
    for (int i = 0; i < vctCloud3D.size(); ++i)
    {
       vctCur = vctCloud3D.at(i)-vctMean;
       arrCur[0] =(double)vctCur.x();
       arrCur[1] =(double)vctCur.y();
       arrCur[2] =(double)vctCur.z();
       memcpy(arrCur1, arrCur, 3 * sizeof(double));
       MtrxMultMatrxTransp(arrCur,3, 1, arrCur1,3, arrMtrx) ;
       MtrxSumMatrx(arrDispMtrx, arrMtrx,3, 3, arrDispMtrx) ;
    }
    MatrxMultScalar(arrDispMtrx, 3, 3, 1./((double)vctCloud3D.size() ),arrDispMtrx);

}
//-------------------------------------------
void QPntCloud ::createClouds_XY_YZ_ZX(QVector<QPointF> *pvctCloud_XY
                ,QVector<QPointF> *pvctCloud_YZ,QVector<QPointF> *pvctCloud_ZX,double *arrMean
                , double *arrDisp,QVector<int> *pvctNumGood, double *parrBuff)
{
    QVector<QVector3D> vctCloud3D(mvctBigMesure.size());
    createCloud_3D(&vctCloud3D, pvctNumGood, parrBuff);

    QVector3D vctMean =createMeanVct3D(vctCloud3D);

    double arrDispMtrx[9] = {0.};
    createDispMtrx3D(vctCloud3D, vctMean,arrDispMtrx);
    arrMean[0] = (double)(vctMean.x());
    arrMean[1] = (double)(vctMean.y());
    arrMean[2] = (double)(vctMean.z());
    arrDisp[0] = arrDispMtrx[0];
    arrDisp[1] = arrDispMtrx[4];
    arrDisp[2] = arrDispMtrx[8];
    pvctCloud_XY->resize(vctCloud3D.size());
    pvctCloud_YZ->resize(vctCloud3D.size());
    pvctCloud_ZX->resize(vctCloud3D.size());

    for (int i =0; i < vctCloud3D.size(); ++i)
    {
        double x = (double)(vctCloud3D.at(i).x());
        double y = (double)(vctCloud3D.at(i).y());
        double z = (double)(vctCloud3D.at(i).z());
        pvctCloud_XY->replace(i,QPointF(x,y));
        pvctCloud_YZ->replace(i,QPointF(y,z));
        pvctCloud_ZX->replace(i,QPointF(x,z));

    }
}
//----------------------------------------------
void QPntCloud :: createPictFiles(wchar_t *wchOutPutFold,const QVector<QBigMeasure> vctBigMeasure
                                  ,const TTable_1D tblEstPrfl,const double* arrSBeacon_XYZ,const double* arrAntPosParams
                                , TYPE_OF_SOLVER_TASK TYPE,double *arrMean, double *arrDisp,QVector<int> *pvctNumGood
                                                         , TYPE_OF_OUTPUT_FILE Type_of_Output_File)
{
    if(wcslen(wchOutPutFold) < 3)
    {
       return;
    }

    wchar_t wchSuff[20] = {0};
    switch (Type_of_Output_File)
    {
    case SHP:
     wcscpy(  wchSuff,  L".shp");
        break;
    case CSV:
     wcscpy(  wchSuff,  L".csv");
        break;
    default:
        break;
    }

    // время измерения
    // невязка по времени
    // КУ теор
    // КУ Зв
    // нев КУ = КУзв -КУ
    // УМтеор
    // УМзв
    // нев УМ
    double * parrBuff  = new double [vctBigMeasure.size() * 10];
    QPntCloud Cloud(vctBigMeasure,tblEstPrfl
                             ,arrSBeacon_XYZ,arrAntPosParams, TYPE);

    QVector<QPointF> vctCloud_XY ,vctCloud_YZ,vctCloud_ZX;

    Cloud.createClouds_XY_YZ_ZX(&vctCloud_XY,&vctCloud_YZ, &vctCloud_ZX
               ,arrMean,arrDisp, pvctNumGood, parrBuff);

    // точки траектории корабля
    const int QuantGoogMeas = pvctNumGood->size();

    TURPointXY *pnts = new TURPointXY[2 *QuantGoogMeas];
    TURPointXY *pntsWave_ = new TURPointXY[QuantGoogMeas];
    TURPointXY *pnts_ = new TURPointXY[QuantGoogMeas];
    for (int i =0; i< QuantGoogMeas; ++i)
    {
        int j = pvctNumGood->at(i);
        pnts[2 *i].X = vctBigMeasure.at(j).marrSVessWaveZv[0] - arrSBeacon_XYZ[0];
        pnts[2 *i ].Y = vctBigMeasure.at(j).marrSVessWaveZv[1]-arrSBeacon_XYZ[1];

        pnts[2 *i +1 ].X = vctBigMeasure.at(j).marrSVessZv[0] -arrSBeacon_XYZ[0];
        pnts[2 *i +1].Y = vctBigMeasure.at(j).marrSVessZv[1]-arrSBeacon_XYZ[1];

        pntsWave_[i].X = vctBigMeasure.at(j).marrSVessWaveZv[0] -arrSBeacon_XYZ[0];
        pntsWave_[i].Y = vctBigMeasure.at(j).marrSVessWaveZv[1] -arrSBeacon_XYZ[1];

        pnts_[i].X = vctBigMeasure.at(j).marrSVessZv[0]-arrSBeacon_XYZ[0];
        pnts_[i].Y = vctBigMeasure.at(j).marrSVessZv[1]-arrSBeacon_XYZ[1];

    }
    //1. оси координат
    wchar_t wchAxesFileName[300] ={0};
    wcscpy(  wchAxesFileName,  wchOutPutFold);
    wcscat(wchAxesFileName, L"//AxesArr");
    wcscat(wchAxesFileName, wchSuff);
    TURPointXY pnt00(arrSBeacon_XYZ[0], arrSBeacon_XYZ[1]);

    const TURPointXY pointBeginX(-10000., 0.),  pointEndX(10000., 0.)
         , pointBeginY(0., -10000.), pointEndY(0., 10000.);
     const double valLength = 1;
    TURPolyLine plnAxes = TURPolyLine::fncCreateAxes( pointBeginX,  pointEndX
                                           ,pointBeginY, pointEndY
                                        , valLength);
    // 1!

    // 2.  все хорошие точки траектории
    wchar_t wchVessPointsFileName[300] ={0};
    wcscpy(  wchVessPointsFileName,  wchOutPutFold);
    wcscat(wchVessPointsFileName, L"//PntsTraj");
     wcscat(wchVessPointsFileName, wchSuff);

    // 2!

    //  3. все хорошие точки на передачу
    wchar_t wchRequestVessPointsFileName[300] ={0};
    wcscpy(   wchRequestVessPointsFileName,  wchOutPutFold);
    wcscat(  wchRequestVessPointsFileName, L"\\RequestPointsTraj");
    wcscat( wchRequestVessPointsFileName, wchSuff);

    // 3!

    // 4. все хорошие точки на прием
    wchar_t wchReceiveVessPointsFileName[300] ={0};
    wcscpy(  wchReceiveVessPointsFileName,  wchOutPutFold);
    wcscat(wchReceiveVessPointsFileName, L"\\ReceievePointsTraj");
    wcscat( wchReceiveVessPointsFileName, wchSuff);


    // 4!

    // ВЫВОД ТРАЕКТОРНЫХ ТОЧЕК
    double* pZ = NULL;
    int  lenVars = 0;
    switch (Type_of_Output_File)
    {
    case SHP:
     pnts_[0].WriteSetSHPFiles(wchReceiveVessPointsFileName,pnts_,QuantGoogMeas);
     pntsWave_[0].WriteSetSHPFiles(  wchRequestVessPointsFileName,pntsWave_,QuantGoogMeas);
     plnAxes.WriteSetSHPFiles(wchAxesFileName,&plnAxes,1);
     pnts[0].WriteSetSHPFiles(wchVessPointsFileName,pnts,2 *QuantGoogMeas);
        break;
    case CSV:

    TURPointXY::PutPointsToCsvFile_(wchReceiveVessPointsFileName,pnts_,QuantGoogMeas);
    TURPointXY::PutPointsToCsvFile_(wchRequestVessPointsFileName,pnts_,QuantGoogMeas);
    TURPointXY::PutPointsToCsvFile_(wchVessPointsFileName,pnts_,QuantGoogMeas);
   // plnAxes.WriteToASCII(wchAxesFileName);

        break;
    default:
        break;
    }
    delete []pnts;
    delete []pntsWave_;
    delete []pnts_;
//

    // ОБЛАКА
    TURPointXY *pntsXY = new TURPointXY[QuantGoogMeas];
    transf_QPointArray_in_TURPointXYArray(pntsXY,vctCloud_XY);

    TURPointXY *pntsXZ= new TURPointXY[QuantGoogMeas];
    transf_QPointArray_in_TURPointXYArray(pntsXZ,vctCloud_ZX);

    TURPointXY *pntsZY = new TURPointXY[QuantGoogMeas];
    transf_QPointArray_in_TURPointXYArray(pntsZY,vctCloud_YZ);

    // 5 облако XY
    wchar_t wchCloudFile_XY[300] ={0};
    wcscpy(  wchCloudFile_XY,  wchOutPutFold);
    wcscat(wchCloudFile_XY, L"\\Cloud_XY");
    wcscat( wchCloudFile_XY, wchSuff);


    // 5!

    // 6 облако XZ
    wchar_t wchCloudFile_XZ[300] ={0};
    wcscpy(  wchCloudFile_XZ,  wchOutPutFold);
    wcscat(wchCloudFile_XZ, L"\\Cloud_XZ");
    wcscat( wchCloudFile_XZ, wchSuff);

    // 6!

    // 7 облако YZ
    wchar_t wchCloudFile_YZ[300] ={0};
    wcscpy(  wchCloudFile_YZ,  wchOutPutFold);
    wcscat(wchCloudFile_YZ, L"\\Cloud_"
                            "YZ");
    wcscat(wchCloudFile_YZ, wchSuff);

    // 7!

    // 8 эллипс XY
    double arrElK_XY[4 ] ={0.};
    arrElK_XY[0] = arrDisp[0];
    arrElK_XY[3] = arrDisp[1];
    TURPointXY pointCentre(arrMean[0],arrMean[1]);
    TURPolygon   ellXY =TURPolygon::fncCreateEllipse_( pointCentre
                                                       ,2. * sqrt(arrDisp[0]),2. * sqrt(arrDisp[1]), 501);
   wchar_t wchEll_XY[300] ={0};
    wcscpy(  wchEll_XY,  wchOutPutFold);
    wcscat(wchEll_XY, L"\\ElXY");
    wcscat(wchEll_XY, wchSuff);

     // 8 эллипс XZ
    double arrElK_XZ[4 ] ={0.};
    arrElK_XZ[0] = arrDisp[0];
    arrElK_XZ[3] = arrDisp[2];
    pointCentre = TURPointXY (arrMean[0],arrMean[2]);
    TURPolygon   ellXZ =TURPolygon::fncCreateEllipse_( pointCentre
                                                      ,2. * sqrt(arrDisp[0]),2. * sqrt(arrDisp[2]), 501);
     wchar_t wchEll_XZ[300] ={0};
    wcscpy(  wchEll_XZ,  wchOutPutFold);
    wcscat(wchEll_XZ, L"\\ElXZ");
    wcscat(wchEll_XZ, wchSuff);

     // 8 эллипс YZ
    double arrElK_ZY[4 ] ={0.};
    arrElK_ZY[0] = arrDisp[1];
    arrElK_ZY[3] = arrDisp[2];
     pointCentre = TURPointXY (arrMean[1],arrMean[2]);
    TURPolygon   ellZY =TURPolygon::fncCreateEllipse_( pointCentre
                                                       ,2. * sqrt(arrDisp[1]),2. * sqrt(arrDisp[2]), 501);
    wchar_t wchEll_ZY[300] ={0};
    wcscpy(  wchEll_ZY,  wchOutPutFold);
    wcscat(wchEll_ZY, L"\\ElZY");
    wcscat(wchEll_ZY, wchSuff);

    // 8!

    TURPointXY pntMeanXY(arrMean[0], arrMean[1]);
    TURPointXY pntMeanXZ(arrMean[0], arrMean[2]);
    TURPointXY pntMeanZY(arrMean[1], arrMean[2]);

    // 9. средняя точка XY
    wchar_t wchMean_XY[300] ={0};
    wcscpy(  wchMean_XY,  wchOutPutFold);
    wcscat(wchMean_XY, L"\\Mean_XY");
    wcscat(wchMean_XY, wchSuff);

    // 9!

    //10. средняя точка XZ
    wchar_t wchMean_XZ[300] ={0};
    wcscpy(  wchMean_XZ,  wchOutPutFold);
    wcscat(wchMean_XZ, L"\\Mean_XZ");
    wcscat(wchEll_XZ, wchSuff);

    // 10!

    //11. средняя точка ZY
    wchar_t wchMean_ZY[300] ={0};
    wcscpy( wchMean_ZY,  wchOutPutFold);
    wcscat(wchMean_ZY, L"\\Mean_YZ");
    wcscat(wchMean_ZY, wchSuff);

    // 11!
    switch (Type_of_Output_File)
    {
    case SHP:
     pntMeanZY.WriteSetSHPFiles(wchMean_ZY,&pntMeanZY,1);

     pntMeanXZ.WriteSetSHPFiles(wchEll_XZ,&pntMeanXZ,1);

     pntMeanXY.WriteSetSHPFiles(wchMean_XY,&pntMeanXY,1);

     ellXY.WriteSetSHPFiles(wchEll_ZY,&ellZY,1);

     ellXY.WriteSetSHPFiles(wchEll_XY,&ellXY,1);

     ellXZ.WriteSetSHPFiles(wchEll_XZ,&ellXZ,1);

     pntsXY[0].WriteSetSHPFiles(wchCloudFile_YZ,pntsZY,QuantGoogMeas);

     pntsXY[0].WriteSetSHPFiles(wchCloudFile_XZ,pntsXZ,QuantGoogMeas);

     pntsXY[0].WriteSetSHPFiles(wchCloudFile_XY,pntsXY,QuantGoogMeas);

        break;
    case CSV:

    TURPointXY::PutPointsToCsvFile_(wchMean_ZY,&pntMeanZY,1);
    TURPointXY::PutPointsToCsvFile_(wchMean_XY,&pntMeanXY,1);
    TURPointXY::PutPointsToCsvFile_(wchMean_ZY,&pntMeanZY,1);

    TURPointXY::PutPointsToCsvFile_(wchEll_XY,ellXY.Points,ellXY.NumPoints);
    TURPointXY::PutPointsToCsvFile_(wchEll_XZ,ellXZ.Points,ellXZ.NumPoints);
    TURPointXY::PutPointsToCsvFile_(wchEll_ZY,ellZY.Points,ellZY.NumPoints);

    TURPointXY::PutPointsToCsvFile_(wchCloudFile_YZ,pntsZY,QuantGoogMeas);
    TURPointXY::PutPointsToCsvFile_(wchCloudFile_XZ,pntsXZ,QuantGoogMeas);
    TURPointXY::PutPointsToCsvFile_(wchCloudFile_XY,pntsXY,QuantGoogMeas);
        break;

    default:
        break;
    }
    delete []pntsXY;
    delete []pntsXZ;
    delete []pntsZY;

    if (Type_of_Output_File == CSV)
    {
        int iNumCols = 10;
        int iLenName = 30;
        wchar_t wchBuffResiduals[300] ={0};
        wcscpy( wchBuffResiduals,  wchOutPutFold);
        wcscat(wchBuffResiduals, L"\\Residuals.csv");

        wchar_t pwcharrColNames[10 * 30] = {0};
        memset(pwcharrColNames, 0,iNumCols * 30 * sizeof(wchar_t));
        wcscpy( pwcharrColNames,  L"ZamerT");
        wcscpy( &pwcharrColNames[iLenName],  L"ResidT");
        wcscpy( &pwcharrColNames[iLenName *2],  L"q_Theor");
        wcscpy( &pwcharrColNames[iLenName *3],  L"q_Zv");
        wcscpy( &pwcharrColNames[iLenName *4],  L"q_Resid");
        wcscpy( &pwcharrColNames[iLenName *5],  L"e_Theor");
        wcscpy( &pwcharrColNames[iLenName *6],  L"e_Zv");
        wcscpy( &pwcharrColNames[iLenName *7],  L"e_Resid");
        wcscpy( &pwcharrColNames[iLenName *8],  L"T_izm");
        wcscpy( &pwcharrColNames[iLenName *9],  L"T_est");


        TYrWrite::WriteMassiveInFIleSCV(wchBuffResiduals,parrBuff, vctBigMeasure.size(), iNumCols
                                     ,NULL,pwcharrColNames, iLenName);
    }
    else
    {
        int iNumCols = 10;
        int iLenName = 30;

        wchar_t *pwcharrColNames =new wchar_t[iNumCols* iLenName];
        memset(pwcharrColNames, 0, iLenName * iNumCols * sizeof(wchar_t));
        wcscpy( pwcharrColNames,  L"ZamerT.shp");
        wcscpy( &pwcharrColNames[iLenName],  L"ResidT.shp");
        wcscpy( &pwcharrColNames[iLenName *2],  L"q_Theor.shp");
        wcscpy( &pwcharrColNames[iLenName *3],  L"q_Zv.shp");
        wcscpy( &pwcharrColNames[iLenName *4],  L"q_Resid.shp");
        wcscpy( &pwcharrColNames[iLenName *5],  L"e_Theor.shp");
        wcscpy( &pwcharrColNames[iLenName *6],  L"e_Zv.shp");
        wcscpy( &pwcharrColNames[iLenName *7],  L"e_Resid.shp");
        wcscpy( &pwcharrColNames[iLenName *8],  L"T_izm.shp");
        wcscpy( &pwcharrColNames[iLenName *9],  L"T_est.shp");

         for (int i = 1; i < iNumCols; ++i)
         {
           TYrWriteShapeFile::WriteOneReport_Points(wchOutPutFold  // путь к папке
                                              ,parrBuff// массив с информацией - матрица nBuffRows x nBuffCols
                                              ,iNumCols // - к-во переменных о корорых накоплена информация в буфере
                                              ,vctBigMeasure.size() //  - к-во точек
                                              ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                              ,iLenName// максимальная длина имени переменной
                                              ,0  // номер переменной по оси X
                                              ,i // номер переменной по оси Y
                                              ,1.  //  масштаб по оси X
                                              ,1./180.* M_PI* 1000.  // масштаб по оси Y
                                               );
         }
         delete []pwcharrColNames;
    }
    delete []parrBuff;

}
//----------------------------------------
void QPntCloud ::transf_QPointArray_in_TURPointXYArray(TURPointXY *pPntArr,const QVector<QPointF> vctCloud)
{
    for (int i =0;i < vctCloud.size(); ++i)
    {
        pPntArr[i].X =  vctCloud.at(i).x();

        pPntArr[i].Y =  vctCloud.at(i).y();
    }
}
