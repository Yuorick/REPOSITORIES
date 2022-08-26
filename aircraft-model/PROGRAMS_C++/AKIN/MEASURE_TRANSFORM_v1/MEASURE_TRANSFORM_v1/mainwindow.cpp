#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <float.h>
#include <wchar.h>
#include <math.h>
#include <string.h>
#include <QFileDialog>
#include <QMessageBox>
#include <QVector>
#include <stdlib.h>


#include "YrWriteShapeFile.h"
#include "dir.h"
#include "Comp.h"
#include "Equations.h"

#include "MatrixProccess.h"


#include "SubWaterBeam.h"
//#include "BigMeasure.h"
#include "Gps.h"
#include "CoordSystTrsf.h"







#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsReport0;



//extern const bool BEZ_SHUMOV = true;
extern  bool BEZ_SHUMOV = true;
//шаг фильрации




MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

     ui->setupUi(this);
     this->ui->lineEdit_3->setText("D:\\AKIN\\PROFILES\\TEST.000");
     QString strFold = this->ui->lineEdit_3->text();
     wchar_t array [400] = {0};
     strFold.toWCharArray(array);
     array[strFold.length()] = 0;
     wcscpy(mpwchPrflFIle, array);
     //QSubWaterBeam::createProfileTbl(mpwchPrflFIle, &mtblEstPrfl);
    // mDeepthMax = mtblEstPrfl.mparrArg[mtblEstPrfl.mNumCols - 2];
    // ui->doubleSpinBox_32->setValue(mDeepthMax);

    //mZonaR = calcXMCrit( 2., mDeepthMax, mtblEstPrfl);

     ui->doubleSpinBox_33->setValue(mZonaR);


  // установка параметров позиционирования антенны
  memset(marrTruePosParams, 0, 6 * sizeof(double));
  memset(marrAprioriPosParams, 0, 6 * sizeof(double));

  marrTruePosParams[0] = marrAprioriPosParams[0] = 1.;
  marrTruePosParams[1] = marrAprioriPosParams[1] = 10.;
  marrTruePosParams[2] = marrAprioriPosParams[2] = -8.;

  double arr[12] = {0.};
  for (int j =0; j < 6; ++j)
  {
      arr[j] = marrTruePosParams[j];
      arr[6 + j] = marrAprioriPosParams[j];
  }

  //--------------------------------------------------------------

    mbtableWidget_4Init = false;
    for (int i=0; i<  ui->tableWidget_4->rowCount() ; i++)
                   for (int j = 0; j < ui->tableWidget_4->columnCount(); j++)
                   {
                       QTableWidgetItem* ptwi0 = nullptr;
                       ptwi0 = new QTableWidgetItem(QString::number(0.));
                       ui->tableWidget_4->setItem(i,j,ptwi0);
                   }


    for (int j =0; j < ui->tableWidget_4->columnCount(); ++j)
    {
        int ia = 10000.* marrTruePosParams[j];
        double temp = ia/10000.;
     ui->tableWidget_4->item(0,j)->setText(QString::number(temp));

     ia = 10000.* marrAprioriPosParams[j];
     temp = ia/10000.;
     ui->tableWidget_4->item(1,j)->setText(QString::number(temp));
    }

    mbtableWidget_4Init = true;
    /// !

    //--------------------------------------------------------------
    // Таблица маяка
    // установка вектора маяка
    memset(marrTrueBeacPos, 0, 3 * sizeof(double));
    memset(marrTrueBeacVelo, 0, 3 * sizeof(double));


    marrTrueBeacPos[0]  = 0.;
    marrTrueBeacPos[1]  = 0.;
    marrTrueBeacPos[2] = -99.;
    marrTrueBeacVelo[1]  = -50.;

    mbtableWidget_6Init = false;
    for (int i=0; i<  ui->tableWidget_6->rowCount() ; i++)
                   for (int j = 0; j < ui->tableWidget_6->columnCount(); j++)
                   {
                       QTableWidgetItem* ptwi0 = nullptr;
                       ptwi0 = new QTableWidgetItem(QString::number(0.));
                       ui->tableWidget_6->setItem(i,j,ptwi0);
                   }


    for (int j =0; j < ui->tableWidget_6->columnCount(); ++j)
    {
     ui->tableWidget_6->item(0,j)->setText(QString::number(marrTrueBeacPos[j]));
     ui->tableWidget_6->item(1,j)->setText(QString::number(marrTrueBeacVelo[j]));
    }

    mbtableWidget_6Init = true;
    /// !


    // ТАБЛИЦА РЕЗУЛЬТАТОВ ПОЛОЖЕНИЯ КОРАБЛЯ
  ui->tableWidget_5->setColumnCount(7);
  ui->tableWidget_5->setRowCount(2 );
  ui->tableWidget_5->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_5->verticalHeader()->setVisible(true);

  QStringList lst5;
  lst5<<"T"<<"X"<<"Y"<<"Z"<<"Q"<<"Psi"<<"Tet";
  ui->tableWidget_5->setHorizontalHeaderLabels(lst5);
  QStringList lst6;
  lst6<< "Передача"<<"Прием";
  ui->tableWidget_5->setVerticalHeaderLabels(lst6);
  QString qstr00 = QString::number(0.);

  for (int i = 0; i < (ui->tableWidget_5->rowCount()); ++i)
      for (int j = 0; j < (ui->tableWidget_5->columnCount()); ++j)
      {
         QTableWidgetItem* ptwi = nullptr;
         ptwi = new QTableWidgetItem(qstr00);
         ui->tableWidget_5->setItem(i,j,ptwi);
      }


  //---------------------------------------------
      // ТАБЛИЦА РЕЗУЛЬТАТОВ ПО МАЯКУ
  ui->tableWidget_7->setColumnCount(4);
  ui->tableWidget_7->setRowCount(2 );
  ui->tableWidget_7->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_7->verticalHeader()->setVisible(true);

  QStringList lst7;
  lst7<<"T"<<"X"<<"Y"<<"Z";
  ui->tableWidget_7->setHorizontalHeaderLabels(lst7);
  QStringList lst8;
  lst8<< "Ист"<<"Оц";
  ui->tableWidget_7->setVerticalHeaderLabels(lst8);
  QString qstr000 = QString::number(0.);

  for (int i = 0; i < (ui->tableWidget_7->rowCount()); ++i)
      for (int j = 0; j < (ui->tableWidget_7->columnCount()); ++j)
      {
         QTableWidgetItem* ptwi = nullptr;
         ptwi = new QTableWidgetItem(qstr000);
         ui->tableWidget_7->setItem(i,j,ptwi);
      }

  //---------------------------------------------

  //--------------------------------------------------------------
  // Таблица GPS
  // установка вектора маяка
  marrGPSTruePosParams[0] = 2. ;
  marrGPSTruePosParams[1] = 3.;
  marrGPSTruePosParams[2] = 4.;

  marrGPSAprioriPosParams[0] = 2.;//3. ;
  marrGPSAprioriPosParams[1] = 3.;//4.;
  marrGPSAprioriPosParams[2] =4.;// 5.;


  mbtableWidget_8Init = false;
  for (int i=0; i<  ui->tableWidget_8->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget_8->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = nullptr;
                     ptwi0 = new QTableWidgetItem(QString::number(0.));
                     ui->tableWidget_8->setItem(i,j,ptwi0);
                 }


  for (int j =0; j < ui->tableWidget_8->columnCount(); ++j)
  {
   ui->tableWidget_8->item(0,j)->setText(QString::number(marrGPSTruePosParams[j]));
   ui->tableWidget_8->item(1,j)->setText(QString::number(marrGPSAprioriPosParams[j]));
  }

  mbtableWidget_8Init = true;
  /// !

  mType_of_000 = VAR0;
  ui->comboBox_7->setCurrentIndex(0);

  //--------------------------------------------------------------

    mbtableWidget_9Init = false;
    for (int i=0; i<  ui->tableWidget_9->rowCount() ; i++)
                   for (int j = 0; j < ui->tableWidget_9->columnCount(); j++)
                   {
                       QTableWidgetItem* ptwi0 = nullptr;
                       ptwi0 = new QTableWidgetItem(QString::number(0.));
                       ui->tableWidget_9->setItem(i,j,ptwi0);
                   }

    double arrVessTable[5] = {0., 200., 0., 90., 1.};
    for (int j =0; j < ui->tableWidget_9->columnCount(); ++j)
    {
        int ia = 10000.* arrVessTable[j];
        double temp = ia/10000.;
     ui->tableWidget_9->item(0,j)->setText(QString::number(temp));

    }

    mbtableWidget_9Init = true;
    /// !

   // ui->frame_6->hide();
}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------


void MainWindow:: inputData()
{

    ///
    // 0.
    QString strFold = this->ui->lineEdit_3->text();
    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPrflFIle, array);



      // без шумов
   if(ui->checkBox->isChecked()==true)
   {
      BEZ_SHUMOV = true;
   }
   else
   {
    BEZ_SHUMOV = false;
   }
   //создание профиля скорости звука
   QSubWaterBeam::createProfileTbl(mpwchPrflFIle, &mtblEstPrfl,mType_of_000);
  // 3!

   // 4. создание реального профиля звука
   mtblRealPrfl = mtblEstPrfl;
   switch(ui->comboBox_6->currentIndex())
   {
   case 0:
       break;
   case 1:
       mtblRealPrfl.multValue(1.01);
        break;
   case 2:
       mtblRealPrfl.multValue(0.99);
       break;
   case 3:
       mtblRealPrfl.multRandValue(0.01);
       break;
   case 4:
       mtblEstPrfl = mtblRealPrfl.makeMiddleProfile();
       break;
   default:
       break;
   }



   for (int i = 0; i < 6; ++i)
   {
       if (i <3)
       {
           marrTruePosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();
           marrAprioriPosParams[i] =  ui->tableWidget_4->item(1,i)->text().toDouble();
       }
       else
       {
           marrTruePosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;
           marrAprioriPosParams[i] =  ui->tableWidget_4->item(1,i)->text().toDouble()* M_PI/ 180.;
       }
   }
   ///

   // создание корабля
       // атмосфера
       double valWind_V = ui->doubleSpinBox_30->value();

       double valWind_VertV = 0.;

       double valWind_Alf =  0.;
       mEnvironment  = TEnvironment ( valWind_V , valWind_Alf, valWind_VertV );
       ///

       // СИНС
       // СКЗ угловой ошибки СИНС (углов качек)
       double valSigSinsQ = ui->doubleSpinBox_9->value()/ 1000.;

       double valSigSinsPsi = ui->doubleSpinBox_10->value()/ 1000.;

       // СКЗ угловой ошибки СИНСпо скорости углов (углов качек)
        double valSig_d_po_dt_Sins = ui->doubleSpinBox_10->value()/1000.;
         double valMaxSig_Q      =     valSigSinsQ ; //0.000582;
         double valMaxSig_Psi    =      valSigSinsPsi ; //0.00145;
         double valMaxSig_Tet    =    valSigSinsPsi ; // 0.00145;
         double valMaxSig_dQdt   =      valSig_d_po_dt_Sins ;
         double valMaxSig_dPsidt =      valSig_d_po_dt_Sins ;
         double valMaxSig_dTetdt =      valSig_d_po_dt_Sins ;
         double valK1         = 0.01 ;
         double valSigV  = 0.2;
         double valSigH       =      0.1 ;
         double valMaxSig_H =     0.1 ;
         double valMaxSig_VH =     0.05 ;



         mT0 = 0.;
        mTimeTempSins = 1./ 50.;
        mSins = QPeaceSins (mEnvironment,valMaxSig_Q, valMaxSig_Psi, valMaxSig_Tet
                            ,valMaxSig_dQdt, valMaxSig_dPsidt,valMaxSig_dTetdt
                            ,valMaxSig_H,valMaxSig_VH,valK1
                            ,valSigV, mT0,mTimeTempSins );
       ///
            // 3.4  параметры корабля
        double arrVess[5] ={0.};
        read_tableWidget_9(arrVess);
          // константы
       double valVesselWidth = 40; // ширина(м)
       double valVesselLength = 200.; // длина(м)


       // параметры движения
       double valQ0 = arrVess[3] * M_PI/ 180.; // генеральный курс

       double valVVess = arrVess[4];
       mvalMaxQ =    3./180.*M_PI; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
       double valT_Q = 18.; // период рыскания
       mvalMaxPsi =      3./180.*M_PI;// максимальный угол килевой качки(амплитуда)
       double valT_Psi =  12; // период килевой качки
       mvalMaxTet =      12./180.*M_PI; //максимальный угол боротовой качки(амплитуда)
       double valT_Tet =  6.; // период бортовой качки
       double valMaxVert =     1. ;


       //

       mRls_Usbl3D = QRls_Usbl3D((ui->doubleSpinBox_8->value()) * 0.001,0.04
         ,(ui->doubleSpinBox_7->value()) * 0.001,(ui->doubleSpinBox_13->value()) * 0.001);
       mpHidroRLS = &mRls_Usbl3D;
       //mTypeOfRLS = USBL_3D;

       mPlatform = QPlatform(mpHidroRLS, marrTruePosParams);
       //

       //!
       // GPS
       read_tableWidget_8(marrGPSTruePosParams,marrGPSAprioriPosParams);
        mGPS = QGps(ui->doubleSpinBox_31->value(),marrGPSTruePosParams);

       //!

       mVess =  QPeaceVess (mEnvironment,valVesselWidth,valVesselLength
             ,mvalMaxQ ,valT_Q,mvalMaxPsi,valT_Psi ,mvalMaxTet
             ,valT_Tet,valMaxVert, valQ0,valVVess
             , mSins,  mPlatform, mT0,mGPS, arrVess);
//!

     mTObrabotki = ui->doubleSpinBox_34->value() ;

     double arrTrueBeaconVS[6];
     read_tableWidget_6(arrTrueBeaconVS);
     memcpy(marrTrueBeacPos, arrTrueBeaconVS, 3 * sizeof(double));
     memcpy(marrTrueBeacVelo, &arrTrueBeaconVS[3], 3 * sizeof(double));

     mQuantMeas = 1;

     mImitMod = QImitMod(mVess, mtblRealPrfl
                         ,marrTrueBeacPos, marrTrueBeacVelo
                         , mTypeOfVessTraj, mTObrabotki,marrGPSAprioriPosParams
                         ,marrTrueBeacPos, marrTrueBeacVelo);

}


//--------------------------------------------
void MainWindow::on_pushButton_3_clicked()
{
    inputData();


  mVess.createTrueMeasureParams( mtblRealPrfl, mTObrabotki
                             , marrTrueBeacPos, marrTrueBeacVelo,mTrueMeasParams);




    mImitMod.imitateMeasure(mTrueMeasParams, mBigMeasure);

    double arrSZv[3] ={0.}, valTZv = 0.;

    switch(ui->comboBox_8->currentIndex())
    {
    case 0: // в ГСК
    // transfMeasure_to_GSK(mtblEstPrfl,mBigMeasure,marrAprioriPosParams,  arrSZv,   &valTZv);
   transfMeasure_to_GSK(mtblEstPrfl,mBigMeasure,marrAprioriPosParams,marrTrueBeacVelo,  arrSZv,   &valTZv);

        break;

    case 1:
        // пересчет замера 3D в АСПК
        transfMeasure_to_ASPK(mtblEstPrfl,mBigMeasure,marrAprioriPosParams,  arrSZv,   &valTZv);
        break;
    default:

        break;
    }

    // ОТЛАДКА ОБРАТНОГО ПРЕОБРАЗОВАНИЯ
    double val_q = 0., val_e= 0.,  val_t= 0.;
    transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, arrSZv, mBigMeasure.marrSVessZv, mBigMeasure.marrMuZv
      ,marrAprioriPosParams,  &val_q, &val_e,  &val_t);
    int uu = 0;
    // ! ОТЛАДКА



    QBigMeasure BigMeasureTrue (mTrueMeasParams.marrSVessWave,mTrueMeasParams.marrMuWaveZv,
                                mTrueMeasParams.marrSVess,mTrueMeasParams.marrMuZv,
                                mTrueMeasParams.mTzapr, mTrueMeasParams.mTotv,
                                mTObrabotki,3
                                ,mTrueMeasParams.mq, mTrueMeasParams.me
                                ,0.,0.,0.);

    double arrS[3] ={0.}, valT = 0.;

    switch(ui->comboBox_8->currentIndex())
    {
    case 0: // в ГСК
        //transfMeasure_to_GSK(mtblRealPrfl,BigMeasureTrue,marrTruePosParams,  arrS,   &valT);
        transfMeasure_to_GSK(mtblRealPrfl,BigMeasureTrue,marrTruePosParams,marrTrueBeacVelo,  arrS,   &valT);
        break;

    case 1:
        // пересчет замера 3D в АСПК
        transfMeasure_to_ASPK(mtblRealPrfl,BigMeasureTrue,marrTruePosParams,  arrS,   &valT);
        break;
    default:

        break;
    }

    double arrVessPosParams[14]= {0.};
    arrVessPosParams[0] = mTrueMeasParams.mTzapr;
    arrVessPosParams[7] = mTrueMeasParams.mTotv;
    memcpy(&arrVessPosParams[1],mTrueMeasParams.marrSVessWave, 3 * sizeof(double) );
    memcpy(&arrVessPosParams[8],mTrueMeasParams.marrSVess, 3 * sizeof(double) );
    memcpy(&arrVessPosParams[4],mTrueMeasParams.marrMuWaveZv, 3 * sizeof(double) );
    memcpy(&arrVessPosParams[11],mTrueMeasParams.marrMuZv, 3 * sizeof(double) );
    write_tableWidget_5(arrVessPosParams);


    double arrZamerOut[8]= {0.};
    arrZamerOut[0] = valT;
    memcpy(&arrZamerOut[1],arrS, 3 * sizeof(double) );

    arrZamerOut[4] = valTZv;
    memcpy(&arrZamerOut[5],arrSZv, 3 * sizeof(double) );

    write_tableWidget_7(arrZamerOut);


}
//-------------------------------------------------
void MainWindow::read_tableWidget_4(double *arr)
{
for (int i=0; i < ui->tableWidget_4->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_4->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_4->columnCount()) + j] =  ui->tableWidget_4->item(i,j)->text().toDouble();
      }
}
//-------------------------------------------------
void MainWindow::write_tableWidget_4(double *arr)
{
    if (!mbtableWidget_4Init)
    {
        return;
    }
   int ncols = ui->tableWidget_4->columnCount();
    for (int i=0; i < ui->tableWidget_4->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_4->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}
//-------------------------------------------
//------------------------------------------
void MainWindow::on_tableWidget_6_cellChanged(int row, int column)
{
    if (!mbtableWidget_6Init)
    {
        return;
    }
    double arr[3] = {0.};


}
//-------------------------------------------------
void MainWindow::read_tableWidget_6(double *arr)
{
for (int i=0; i < ui->tableWidget_6->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_6->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_6->columnCount()) + j] =  ui->tableWidget_6->item(i,j)->text().toDouble();
      }
}

//-------------------------------------------------
void MainWindow::read_tableWidget_9(double *arr)
{
for (int i=0; i < ui->tableWidget_9->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_9->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_9->columnCount()) + j] =  ui->tableWidget_9->item(i,j)->text().toDouble();
      }
}

//-------------------------------------------------
void MainWindow::read_tableWidget_8(double *arr1,double *arr2)
{
     for (int i=0; i < ui->tableWidget_8->columnCount(); ++i)
    {
          arr1 [i] =  ui->tableWidget_8->item(0,i)->text().toDouble();
          arr2 [i] =  ui->tableWidget_8->item(1,i)->text().toDouble();
    }
}

//-------------------------------------------------
void MainWindow::write_tableWidget_6(double *arr)
{
    if (!mbtableWidget_6Init)
    {
        return;
    }
   int ncols = ui->tableWidget_6->columnCount();
    for (int i=0; i < ui->tableWidget_6->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_6->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}


//--------------------------------------------
/*
int MainWindow::createInputDataReport(wchar_t*FileName, const bool bHeader, QTargTrackingOperation &TargTrackingOperation)
{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {

      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {

     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"  ОБЩИЕ НАСТРОЙКИ АЛГОРИТМА\n");


if(BEZ_SHUMOV )
{
  fprintf(fw,"  без шумов - ДА\n");
}
else
{
  fprintf(fw,"  без шумов - НЕТ\n");
}


fclose(fw);

}
//-------------------------------------------------------------
int MainWindow::createOutputDataReport(wchar_t*FileName, const bool bHeader
               ,const double* arrHorSyst,const double* arrHorDisp
               ,const double*  arrVertSyst,const double* arrVertDisp)
{
    int len = wcslen(FileName) ;

    if ( !( (FileName[len - 1] == 't') && (FileName[len - 2] == 'x') // проверка, что
     && (FileName[len - 3] == 't') ) )  // указанный файл имеет расширение .flt
    {

      return 1 ;
    }

    FILE *fw ;

    if ((fw = _wfopen(FileName,L"a"))== NULL)

    {

     return 1 ;
    }
if (bHeader)
{
   fprintf(fw,"  Дата и время формирования отчета\n");
   time_t t = time(NULL);
   struct tm* aTm = localtime(&t);
   fprintf(fw,"  Год = %04d\n",aTm->tm_year+1900);
   fprintf(fw,"  Mесяц = %02d\n",aTm->tm_mon+1);
   fprintf(fw,"  День = %02d\n",aTm->tm_mday);
   fprintf(fw,"  Время = %02d:%02d:%02d\n",aTm->tm_hour, aTm->tm_min, aTm->tm_sec);
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
   fprintf(fw,"***************************************\n");
}
fprintf(fw,"***************************************\n");
fprintf(fw,"  РЕЗУЛЬТАТЫ РАБОТЫ ИМТАЦИОННОЙ МОДЕЛИ ПРИВОДА ПЛАТФОРМЫ\n");
fprintf(fw,"\n");
fprintf(fw,"\n");
fprintf(fw,"           ОШИБКИ ГОРИЗОНТАЛЬНОГО ПРИВОДА, мрад\n");
fprintf(fw,"                      Сист.            СКЗ флукт.       СКЗ сумм \n");
fprintf(fw,"угловое положение:    %8.7f       %8.7f        %8.7f \n", arrHorSyst[0] * 1000., sqrt(arrHorDisp[0]) * 1000.,sqrt(arrHorDisp[0] + arrHorSyst[0]* arrHorSyst[0]) * 1000.);
fprintf(fw,"угловая скорость :    %8.7f       %8.7f        %8.7f \n", arrHorSyst[1] * 1000., sqrt(arrHorDisp[1]) * 1000.,sqrt(arrHorDisp[1] + arrHorSyst[1]* arrHorSyst[1]) * 1000.);

fprintf(fw,"\n");
fprintf(fw,"\n");
fprintf(fw,"               ОШИБКИ ВЕРТИКАЛЬНОГО ПРИВОДА, мрад\n");
fprintf(fw,"\n");
fprintf(fw,"                      Сист.            СКЗ флукт.       СКЗ сумм \n");
fprintf(fw,"угловое положение:    %8.7f       %8.7f        %8.7f \n", arrVertSyst[0] * 1000., sqrt(arrVertDisp[0]) * 1000.,sqrt(arrVertDisp[0] + arrVertSyst[0]* arrVertSyst[0]) * 1000.);
fprintf(fw,"угловая скорость :    %8.7f       %8.7f        %8.7f \n", arrVertSyst[1] * 1000., sqrt(arrVertDisp[1]) * 1000.,sqrt(arrVertDisp[1] + arrVertSyst[1]* arrVertSyst[1]) * 1000.);
fclose(fw);

}
*/

void MainWindow::on_pushButton_5_clicked()
{
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости","D:\\AKIN\\PROFILES" ,"*.000");
    this->ui->lineEdit_3->setText(strFold);

    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPrflFIle, array);
    QSubWaterBeam::createProfileTbl(mpwchPrflFIle, &mtblEstPrfl, mType_of_000);
    mDeepthMax = mtblEstPrfl.mparrArg[mtblEstPrfl.mNumCols - 2];
    ui->doubleSpinBox_32->setValue(mDeepthMax);

    mZonaR = calcXMCrit( 2., mDeepthMax, mtblEstPrfl);

    ui->doubleSpinBox_33->setValue(mZonaR);

}
//------------------------------------------------------
//void MainWindow::on_tableWidget_7_itemActivated(QTableWidgetItem *item)
//{
//   read_tableWidget_7(marrMayjak_TblData);
//}
//-------------------------------------------------
/*
void MainWindow::read_tableWidget_7(double *arr)
{
for (int i=0; i < ui->tableWidget_7->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_7->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_7->columnCount()) + j] =  ui->tableWidget_7->item(i,j)->text().toDouble();
      }
}
//-------------------------------------------------
void MainWindow::write_tableWidget_7(double *arr)
{

   int ncols = ui->tableWidget_7->columnCount();
    for (int i=0; i < ui->tableWidget_7->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_7->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}
*/
//void MainWindow::on_tableWidget_7_itemSelectionChanged()
//{
 //  write_tableWidget_7(marrMayjak_TblData);
//}
//---------------------------------

//-------------------------------------------------
void MainWindow::read_tableWidget_5(double *arr)
{
for (int i=0; i < ui->tableWidget_5->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_5->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_5->columnCount()) + j] =  ui->tableWidget_5->item(i,j)->text().toDouble();
      }
}
//-------------------------------------------------
void MainWindow::write_tableWidget_5(double *arr)
{

   int ncols = ui->tableWidget_5->columnCount();
    for (int i=0; i < ui->tableWidget_5->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_5->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}


//-------------------------------------------------
void MainWindow::write_tableWidget_7(double *arr)
{

   int ncols = ui->tableWidget_7->columnCount();
    for (int i=0; i < ui->tableWidget_7->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_7->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}

void MainWindow::on_tableWidget_5_itemSelectionChanged()
{
   //write_tableWidget_5(marrGRLS_TblData);
}
/*
void MainWindow::on_pushButton_clicked()
{
    const double VAlC0 = 1400.;
    const double VAlC1 = 1300.;
    const double VAlC2 = 1299.9;
    const double VAlZn1= 1.;
    const double VAlZn2= 50.;
    const double VAlCosTetta = cos (0.10);
  double valIntA =  calcIntegral_I(VAlC0,VAlZn1
                , VAlZn2,VAlC1, VAlC2,VAlCosTetta);
  double atep = 0.001;
  int iSt = (VAlZn2 - VAlZn1)/atep;
  double valn = VAlC0/ VAlC1;
  double sum = 1./ sqrt(valn * valn - VAlCosTetta* VAlCosTetta)/2.;
  double ccur = VAlC1;
  double valk = (VAlC2 - VAlC1)/(VAlZn2 - VAlZn1);
  for (int i = 1; i < (iSt - 1);++i)
  {
      double valZcur = VAlZn1 + ((double) i) * atep;
      double valCCur = 1300. + valk*(valZcur - VAlZn1);
     valn = VAlC0/ valCCur;
    sum += 1./ sqrt(valn * valn - VAlCosTetta* VAlCosTetta);
  }
  valn = VAlC0/ VAlC2;
  sum +=1./ sqrt(valn * valn - VAlCosTetta* VAlCosTetta)/2.;

  sum = sum *atep;
  int iu = 0;

}
*/

void MainWindow::on_comboBox_7_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:
        mType_of_000 = VAR0;
        break;

    case 1:
        mType_of_000 = VAR1;
        break;

    default:

        break;
    }
}
