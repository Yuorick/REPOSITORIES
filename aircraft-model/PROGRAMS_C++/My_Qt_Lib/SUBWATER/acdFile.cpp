#include "acdFile.h"
#include <QFile>
#include <QTextStream>
#include <math.h>
#include <QDebug>
#include "PeaceVess.h"
#include "MatrixProccess.h"


/*
// INPUT:
// wchReportFile - путь к файлу <*.acd>
// QuantMeas - к-во измерений
// parrMeas - указатель на массив мзмерений
// arrAntPosParams[6] - вектор параметров позиционирования антенны
//     arrAntPosParams[0],arrAntPosParams[1],arrAntPosParams[2]  - офсет в метрах
//      arrAntPosParams[3],arrAntPosParams[4],arrAntPosParams[5] -углы в радианах
//arrGPSPosParams[3] - вектор офсетов GPS
//arrBeaconPos[3] - положение маяка в ГСК (преплот)
//TObrabotki - время обработки сигнала на маяке
//OUTPUT:
//arrAntPosParams[3],arrAntPosParams[4],arrAntPosParams[5] - выводятся в файл в градусах!!!!
bool writeExchangeDataFile(QString wchReportFile,
              const int QuantMeas// к-во измерений
              ,QBigMeasure* parrMeas//
              , double* arrAntPosParams
              ,double* arrGPSPosParams
              ,double* arrBeaconPos
              ,const double TObrabotki
              )
{
    qDebug() << "writeDataFile" << wchReportFile;
    QFile file(wchReportFile);
    if(!file.open(QIODevice::WriteOnly))
    {
        qDebug() << "writeDataFile failed open file";
        //QMessageBox::warning(this, tr("Внимание!"), tr("Ошибка открытия файла %1").arg(file.fileName()));
        return false;
    }
    QTextStream stream(&file);
    double rtog = 180.0 / M_PI;
    int N = QuantMeas;
    stream << "N=" << endl;
    stream << QString::number(N) << endl;
    stream << "P antenna:" << endl;
    stream
            << QString::number(arrAntPosParams[0], 'f', 2) << ";"
            << QString::number(arrAntPosParams[1], 'f', 2) << ";"
            << QString::number(arrAntPosParams[2], 'f', 2) << ";"
            << QString::number(arrAntPosParams[3], 'f', 2) << ";"
            << QString::number(arrAntPosParams[4], 'f', 2) << ";"
            << QString::number(arrAntPosParams[5], 'f', 2)
            << endl;
    stream << "P GPS:" << endl;
    stream
            << QString::number(arrGPSPosParams[0], 'f', 2) << ";"
            << QString::number(arrGPSPosParams[1], 'f', 2) << ";"
            << QString::number(arrGPSPosParams[2], 'f', 2)
            << endl;
    stream << "S beacon:" << endl;
    stream
            << QString::number(arrBeaconPos[0], 'f', 2) << ";"
            << QString::number(arrBeaconPos[1], 'f', 2) << ";"
            << QString::number(arrBeaconPos[2], 'f', 2)
            << endl;
    stream << "T beacon:" << endl;
    stream << QString::number(TObrabotki, 'f', 6) << endl;
    stream << "Table:" << endl;

    QStringList strListDataHeader;
    strListDataHeader << "mTzaprZv";
    strListDataHeader << "marrSVessWaveZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "marrMuWaveZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "mTotvZv";
    strListDataHeader << "marrSVessZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "marrMuZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "q";
    strListDataHeader << "e";
    stream << strListDataHeader.join(";") << endl;

    QStringList strListDatum;
    for(int i = 0; i < N; i++)
    {
        QBigMeasure *measure = &parrMeas[i];
        QString strNoData = "-99999";
        QString strBeaconAzimuth = strNoData, strBeaconElevationAngle = strNoData;
        if(measure->mqzv != NODATA)
            strBeaconAzimuth = QString::number(measure->mqzv*rtog, 'f', 2);
        if(measure->mezv != NODATA)
            strBeaconElevationAngle = QString::number(measure->mezv*rtog, 'f', 2);

        QStringList strListData;
        strListData << QString::number(measure->mTzaprZv, 'f', 6);
        strListData << QString::number(measure->marrSVessWaveZv[0], 'f', 2);
        strListData << QString::number(measure->marrSVessWaveZv[1], 'f', 2);
        strListData << QString::number(measure->marrSVessWaveZv[2], 'f', 2);
        strListData << QString::number(measure->marrMuWaveZv[0]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuWaveZv[1]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuWaveZv[2]*rtog, 'f', 4);
        strListData << QString::number(measure->mTotvZv, 'f', 6);
        strListData << QString::number(measure->marrSVessZv[0], 'f', 2);
        strListData << QString::number(measure->marrSVessZv[1], 'f', 2);
        strListData << QString::number(measure->marrSVessZv[2], 'f', 2);
        strListData << QString::number(measure->marrMuZv[0]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuZv[1]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuZv[2]*rtog, 'f', 4);
        strListData << strBeaconAzimuth;
        strListData << strBeaconElevationAngle;
        strListDatum << strListData.join(";");
    }

    for(int i = 0; i < N; i++)
        stream << strListDatum[i] << endl;

    qDebug() << "writeDataFile finished successfully";
    return true;
}
*/
// INPUT:
// wchReportFile - путь к файлу <*.acd>
// QuantMeas - к-во измерений
// parrMeas - указатель на массив мзмерений
// arrAntPosParams[6] - вектор параметров позиционирования антенны
//     arrAntPosParams[0],arrAntPosParams[1],arrAntPosParams[2]  - офсет в метрах
//      arrAntPosParams[3],arrAntPosParams[4],arrAntPosParams[5] -углы в радианах
//arrGPSPosParams[3] - вектор офсетов GPS
//arrBeaconPos[3] - положение маяка в ГСК (преплот)
//TObrabotki - время обработки сигнала на маяке
//OUTPUT:
//arrAntPosParams[3],arrAntPosParams[4],arrAntPosParams[5] - выводятся в файл в градусах!!!!
bool writeExchangeDataFile_(QString wchReportFile,
              const int QuantMeas// к-во измерений
              ,QBigMeasure* parrMeas//
              , double* arrAntPosParams
              ,double* arrGPSPosParams
              ,double* arrBeaconPos
              ,const double TObrabotki
              )
{
    qDebug() << "writeDataFile" << wchReportFile;
    QFile file(wchReportFile);
    if(!file.open(QIODevice::WriteOnly))
    {
        qDebug() << "writeDataFile failed open file";
        //QMessageBox::warning(this, tr("Внимание!"), tr("Ошибка открытия файла %1").arg(file.fileName()));
        return false;
    }
    QTextStream stream(&file);
    double rtog = 180.0 / M_PI;
    int N = QuantMeas;
    stream << "N=" << endl;
    stream << QString::number(N) << endl;
    stream << "P antenna:" << endl;
    stream
            << QString::number(arrAntPosParams[0], 'f', 2) << ";"
            << QString::number(arrAntPosParams[1], 'f', 2) << ";"
            << QString::number(arrAntPosParams[2], 'f', 2) << ";"
            << QString::number(arrAntPosParams[3]*180./M_PI, 'f', 2) << ";"
            << QString::number(arrAntPosParams[4]*180./M_PI, 'f', 2) << ";"
            << QString::number(arrAntPosParams[5]*180./M_PI, 'f', 2)
            << endl;
    stream << "P GPS:" << endl;
    stream
            << QString::number(arrGPSPosParams[0], 'f', 2) << ";"
            << QString::number(arrGPSPosParams[1], 'f', 2) << ";"
            << QString::number(arrGPSPosParams[2], 'f', 2)
            << endl;
    stream << "S beacon:" << endl;
    stream
            << QString::number(arrBeaconPos[0], 'f', 2) << ";"
            << QString::number(arrBeaconPos[1], 'f', 2) << ";"
            << QString::number(arrBeaconPos[2], 'f', 2)
            << endl;
    stream << "T beacon:" << endl;
    stream << QString::number(TObrabotki, 'f', 6) << endl;
    stream << "Table:" << endl;

    QStringList strListDataHeader;
    strListDataHeader << "mTzaprZv";
    strListDataHeader << "marrSVessWaveZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "marrMuWaveZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "mTotvZv";
    strListDataHeader << "marrSVessZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "marrMuZv[3]";
    strListDataHeader << "";
    strListDataHeader << "";
    strListDataHeader << "q";
    strListDataHeader << "-e";
    stream << strListDataHeader.join(";") << endl;

    QStringList strListDatum;
    for(int i = 0; i < N; i++)
    {
        QBigMeasure *measure = &parrMeas[i];
        QString strNoData = "-99999";
        QString strBeaconAzimuth = strNoData, strBeaconElevationAngle = strNoData;
        if(/*measure->mDimMeas >= 2 && */measure->mqzv != NODATA)
            strBeaconAzimuth = QString::number(measure->mqzv*rtog, 'f', 2);
        if(/*measure->mDimMeas >= 3 && */measure->mezv != NODATA)
            strBeaconElevationAngle = QString::number(-(measure->mezv*rtog), 'f', 2);


        // вычисление вектора положения GPS на момент излучения
        double arrS_Gps_Wave_GSK[3] = {0.}, arrS_Gps_Wave_KGSK[3] = {0.};

         QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrGPSPosParams, measure->marrMuWaveZv
                                                 , NULL,arrS_Gps_Wave_KGSK,3 );
         MtrxSumMatrx(measure->marrSVessWaveZv, arrS_Gps_Wave_KGSK,1, 3, arrS_Gps_Wave_GSK);
         // !

         // вычисление вектора положения GPS на момент приема
         double arrS_Gps_GSK[3] = {0.}, arrS_Gps_KGSK[3] = {0.};

          QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrGPSPosParams, measure->marrMuZv
                                                  , NULL,arrS_Gps_KGSK,3 );
          MtrxSumMatrx(measure->marrSVessZv, arrS_Gps_KGSK,1, 3, arrS_Gps_GSK);
          // !

        QStringList strListData;
        strListData << QString::number(measure->mTzaprZv, 'f', 6);
        strListData << QString::number(arrS_Gps_Wave_GSK[0], 'f', 2);
        strListData << QString::number(arrS_Gps_Wave_GSK[1], 'f', 2);
        strListData << QString::number(arrS_Gps_Wave_GSK[2], 'f', 2);
        strListData << QString::number(measure->marrMuWaveZv[0]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuWaveZv[1]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuWaveZv[2]*rtog, 'f', 4);
        strListData << QString::number(measure->mTotvZv, 'f', 6);
        strListData << QString::number(arrS_Gps_GSK[0], 'f', 2);
        strListData << QString::number(arrS_Gps_GSK[1], 'f', 2);
        strListData << QString::number(arrS_Gps_GSK[2], 'f', 2);
        strListData << QString::number(measure->marrMuZv[0]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuZv[1]*rtog, 'f', 4);
        strListData << QString::number(measure->marrMuZv[2]*rtog, 'f', 4);
        strListData << strBeaconAzimuth;
        strListData << strBeaconElevationAngle;
        strListDatum << strListData.join(";");
    }

    for(int i = 0; i < N; i++)
        stream << strListDatum[i] << endl;

    qDebug() << "writeDataFile finished successfully";
    return true;
}

//----------------------------------
//---------------------------------------------------------------
void readQuantMeasures_From_File_acd(QString DataFileName,int &iQuantMeasures)
{
QFile TargFile;
TargFile.setFileName(DataFileName);
TargFile.open(QIODevice::ReadOnly | QIODevice::Text);
QString line0 = TargFile.readLine();
line0 = TargFile.readLine();

iQuantMeasures = line0.toInt();

TargFile.close();
}
//---------------------------------------
//Чтение файла с исх данными и формиорование массивов отсеков по типам воздействия
void readFile_acd(QString TargFileName,double *arrP_Ant
                  , double *arrP_Gps, double *arrS_Beacon, double &valTBeacon
                  , QVector<QBigMeasure> *pvctBigMesure)
{
    int iQuantMeasures = -1;
    readQuantMeasures_From_File_acd(TargFileName,iQuantMeasures);
    pvctBigMesure->resize(iQuantMeasures);


    QFile TargFile;
    TargFile.setFileName(TargFileName);
    TargFile.open(QIODevice::ReadOnly | QIODevice::Text);

    QString line0 = TargFile.readLine();
    line0 = TargFile.readLine();
    line0 = TargFile.readLine();
    line0 = TargFile.readLine();

    QStringList w =line0.split(";");
    for(int i =0; i < 3;++i )
    {
        arrP_Ant[i] = w.at(i).toDouble();
    }

    for(int i =3; i < 6;++i )
    {
        arrP_Ant[i] = w.at(i).toDouble() * M_PI/ 180.;
    }


    line0 = TargFile.readLine();
    line0 = TargFile.readLine();

    w =line0.split(";");
    for(int i =0; i < 3;++i )
    {
        arrP_Gps[i] = w.at(i).toDouble();
    }

    line0 = TargFile.readLine();
    line0 = TargFile.readLine();

    w =line0.split(";");
    for(int i =0; i < 3;++i )
    {
        arrS_Beacon[i] = w.at(i).toDouble();
    }

    line0 = TargFile.readLine();
    line0 = TargFile.readLine();
    valTBeacon = line0.toDouble();


    QString str = "mTzaprZv";


    while( !TargFile.atEnd())
    {
            line0 = TargFile.readLine();
            if( line0.contains(str))
            {
              break;
            }
     }
    ///
    int iChetchik = 0;
   for (int i =0;i < iQuantMeasures; ++i)
   {
       line0 = TargFile.readLine();
       QStringList w =line0.split(";");
       QBigMeasure measure;
       measure.mTzaprZv = w.at(0).toDouble();

       measure.marrSVessWaveZv[0] = w.at(1).toDouble();
       measure.marrSVessWaveZv[1] = w.at(2).toDouble();
       measure.marrSVessWaveZv[2] = w.at(3).toDouble();

       measure.marrMuWaveZv[0] = w.at(4).toDouble() * M_PI/ 180.;
       measure.marrMuWaveZv[1] = w.at(5).toDouble()* M_PI/ 180.;
       measure.marrMuWaveZv[2] = w.at(6).toDouble()* M_PI/ 180.;

       double arrGps_KGSK[3] = {0.};
       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, measure.marrMuWaveZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(measure.marrSVessWaveZv, arrGps_KGSK,1, 3, measure.marrSVessWaveZv);


       measure.mTotvZv = w.at(7).toDouble();

       measure.marrSVessZv[0] = w.at(8).toDouble();
       measure.marrSVessZv[1] = w.at(9).toDouble();
       measure.marrSVessZv[2] = w.at(10).toDouble();

       measure.marrMuZv[0] = w.at(11).toDouble()* M_PI/ 180.;
       measure.marrMuZv[1] = w.at(12).toDouble()* M_PI/ 180.;
       measure.marrMuZv[2] = w.at(13).toDouble()* M_PI/ 180.;

       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, measure.marrMuZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(measure.marrSVessZv, arrGps_KGSK,1, 3, measure.marrSVessZv);
       measure.mTobr = valTBeacon;
       measure.mqzv = w.at(14).toDouble()* M_PI/ 180.;
       measure.mezv = -w.at(15).toDouble()* M_PI/ 180.;

       pvctBigMesure->replace(i,measure );
       ++iChetchik;
   }

    TargFile.close();
}

//----------------------------------
//---------------------------------------
//Чтение файла с исх данными и формиорование массивов отсеков по типам воздействия
void readBigMeasuresOnly_acd(QString TargFileName,const double *arrP_Ant
                  , double *arrP_Gps, double *arrS_Beacon,const  double valTBeacon
                  , QVector<QBigMeasure> *pvctBigMesure)
{
    int iQuantMeasures = -1;
    readQuantMeasures_From_File_acd(TargFileName,iQuantMeasures);
    pvctBigMesure->resize(iQuantMeasures);


    QFile TargFile;
    TargFile.setFileName(TargFileName);
    TargFile.open(QIODevice::ReadOnly | QIODevice::Text);



    QString str = "mTzaprZv";
    QString line0;


    while( !TargFile.atEnd())
    {
            line0 = TargFile.readLine();
            if( line0.contains(str))
            {
              break;
            }
     }
    ///
    int iChetchik = 0;
   for (int i =0;i < iQuantMeasures; ++i)
   {
       line0 = TargFile.readLine();
       QStringList w =line0.split(";");
       QBigMeasure measure;
       measure.mTzaprZv = w.at(0).toDouble();

       measure.marrSVessWaveZv[0] = w.at(1).toDouble();
       measure.marrSVessWaveZv[1] = w.at(2).toDouble();
       measure.marrSVessWaveZv[2] = w.at(3).toDouble();

       measure.marrMuWaveZv[0] = w.at(4).toDouble() * M_PI/ 180.;
       measure.marrMuWaveZv[1] = w.at(5).toDouble()* M_PI/ 180.;
       measure.marrMuWaveZv[2] = w.at(6).toDouble()* M_PI/ 180.;

      // measure.marrMuWaveZv[1] = -w.at(5).toDouble()* M_PI/ 180.;
      // measure.marrMuWaveZv[2] = w.at(6).toDouble()* M_PI/ 180.;

       double arrGps_KGSK[3] = {0.};
       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, measure.marrMuWaveZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(measure.marrSVessWaveZv, arrGps_KGSK,1, 3, measure.marrSVessWaveZv);


       measure.mTotvZv = w.at(7).toDouble();

       measure.marrSVessZv[0] = w.at(8).toDouble();
       measure.marrSVessZv[1] = w.at(9).toDouble();
       measure.marrSVessZv[2] = w.at(10).toDouble();

       measure.marrMuZv[0] = w.at(11).toDouble()* M_PI/ 180.;
       measure.marrMuZv[1] = w.at(12).toDouble()* M_PI/ 180.;
       measure.marrMuZv[2] = w.at(13).toDouble()* M_PI/ 180.;

      // measure.marrMuZv[1] = -w.at(12).toDouble()* M_PI/ 180.;
      // measure.marrMuZv[2] = w.at(13).toDouble()* M_PI/ 180.;

       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( arrP_Gps, measure.marrMuZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(measure.marrSVessZv, arrGps_KGSK,1, 3, measure.marrSVessZv);
       measure.mTobr = valTBeacon;
       measure.mqzv = w.at(14).toDouble()* M_PI/ 180.;
       measure.mezv = -w.at(15).toDouble()* M_PI/ 180.;

       pvctBigMesure->replace(i,measure );
       ++iChetchik;
   }

    TargFile.close();
}

