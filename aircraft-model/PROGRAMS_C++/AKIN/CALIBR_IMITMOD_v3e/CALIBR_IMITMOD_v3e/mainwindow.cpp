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

#include "YrRead.h"
#include "YrWriteShapeFile.h"
#include "dir.h"
#include "Comp.h"
#include "Equations.h"

#include "MatrixProccess.h"
#include "YrWrite.h"

#include "UrPointXY.h"
#include "URPolyLine.h"

#include "SubWaterBeam.h"
#include "BigMeasure.h"
#include "Gps.h"

#include <QDebug>
#include "acdFile.h"
#include "Gauss.h"






#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsReport0;



//extern const bool BEZ_SHUMOV = true;
extern  bool BEZ_SHUMOV = true;
//шаг фильрации




MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
  myThread(0)
{

     ui->setupUi(this);
     qRegisterMetaType<QVector<double> >("QVector<double>");
     qRegisterMetaType<QVector<QBigMeasure>>("QVector<QBigMeasure>");
     //*********************
     // очистка интерфейса
     ui->pushButtonStartTask->setEnabled(true);
     ui->pushButtonStopTask->setEnabled(false);
     ui->progressBar->setValue(0);
     ui->labelTaskInfo1->setText("");
     ui->labelTaskInfo2->setText("");
     ui->labelTaskInfo3->setText("");
     ui->labelTaskInfo4->setText("");

     switch(ui->comboBox_5->currentIndex())
     {
     case 0:
       ui->groupBox_3 ->hide();
       ui->groupBox->setTitle(QString("РЕШЕНИЕ LBL ЗАДАЧИ"));
         break;

     case 1:
         ui->groupBox_3 ->show();
         ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-2D ЗАДАЧИ"));
         break;

     case 2:
         ui->groupBox_3 ->show();
         ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-3D ЗАДАЧИ"));
         break;
     default: break;
     }

     switch(ui->comboBox_8->currentIndex())
     {
     case 0:
         ui->doubleSpinBox_17->setRange(0.,30.);
         ui->doubleSpinBox_17->setValue(20.);
         ui->doubleSpinBox_17->setPrefix(QString("Gape: "));
         ui->doubleSpinBox_17->setSuffix(QString(" гр"));
         ui->doubleSpinBox_17->setDecimals(1);
         ui->doubleSpinBox_16->setRange(0.,5.);
         ui->doubleSpinBox_16->setValue(1.);
         ui->doubleSpinBox_16->setPrefix(QString("Step: "));
         ui->doubleSpinBox_16->setSuffix(QString(" гр"));
         ui->doubleSpinBox_16->setDecimals(3);
         break;

     case 1:
         ui->doubleSpinBox_17->setRange(0.,1.1);
         ui->doubleSpinBox_17->setValue(0.025);
         ui->doubleSpinBox_17->setPrefix(QString("Коеф: "));
         ui->doubleSpinBox_17->setSuffix(QString(""));
         ui->doubleSpinBox_17->setDecimals(3);
         ui->doubleSpinBox_16->setRange(0.,1000);
         ui->doubleSpinBox_16->setValue(200);
         ui->doubleSpinBox_16->setPrefix(QString("Макс ит.: "));
         ui->doubleSpinBox_16->setSuffix(QString(""));
         ui->doubleSpinBox_16->setDecimals(0);
         break;

     default: break;
     }

     mNumOfTask =1;
     //**********************!
      QString str1 = QCoreApplication::applicationDirPath();
   this->ui->lineEdit_3->setText("D:\\AKIN\\PROFILES\\02102021\\FILE48.000");
    // this->ui->lineEdit_3->setText(str1);
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


   this->ui->lineEdit->setText("D:\\AKIN\\ОТЧЕТ_ПОЗИЦИОНИРОВАНИЕ\\REZ_REPORT");
    // this->ui->lineEdit->setText(str1);

     memset(marrMayjak_TblData, 0, 9 * sizeof(double));
     memset(marrGRLS_TblData, 0, 18 * sizeof(double));


  ui->comboBox_4->setCurrentIndex(5);

  ui->doubleSpinBox_14->setValue(1000);
  ui->doubleSpinBox_14->setReadOnly(false);

  // установка параметров позиционирования антенны
  memset(marrTruePosParams, 0, 6 * sizeof(double));
  memset(marrAprioriPosParams, 0, 6 * sizeof(double));

  marrTruePosParams[0] = marrAprioriPosParams[0] = 5.8;
  marrTruePosParams[1] = marrAprioriPosParams[1] = -7.55;
  marrTruePosParams[2] = marrAprioriPosParams[2] = -5.3;
  marrTruePosParams[5] = marrAprioriPosParams[5] = 30./180.* M_PI;;

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
        double coef = 1.;
        if (j > 2)
        {
          coef = 180./M_PI;
        }
        int ia = 10000.* marrTruePosParams[j] * coef;
        double temp = ia/10000.;
     ui->tableWidget_4->item(0,j)->setText(QString::number(temp));

     ia = 10000.* marrAprioriPosParams[j] * coef;
     temp = ia/10000.;
     ui->tableWidget_4->item(2,j)->setText(QString::number(temp));
    }

    mbtableWidget_4Init = true;
    /// !

    //--------------------------------------------------------------
    // Таблица маяка
    // установка вектора маяка
    memset(marrTrueHeadLightPos, 0, 3 * sizeof(double));
    memset(marrAprioriHeadLightPos, 0, 6 * sizeof(double));

    /*
    marrAprioriHeadLightPos[0] = 10000.;
    marrAprioriHeadLightPos[1] = 500.;
    marrAprioriHeadLightPos[2] = -100.;

    marrTrueHeadLightPos[0]  = 10005.;
    marrTrueHeadLightPos[1]  = 495;
    marrTrueHeadLightPos[2] = -99.;
    */
    marrAprioriHeadLightPos[0] = 0.;
    marrAprioriHeadLightPos[1] = 0.;
    marrAprioriHeadLightPos[2] = -99.;//-100.;

    marrTrueHeadLightPos[0]  = 0.;
    marrTrueHeadLightPos[1]  = 0.;
    marrTrueHeadLightPos[2] = -99.;

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
     ui->tableWidget_6->item(0,j)->setText(QString::number(marrTrueHeadLightPos[j]));
     ui->tableWidget_6->item(1,j)->setText(QString::number(marrAprioriHeadLightPos[j] - marrTrueHeadLightPos[j]));
     ui->tableWidget_6->item(2,j)->setText(QString::number(marrAprioriHeadLightPos[j]));
    }

    mbtableWidget_6Init = true;
    /// !


    // ТАБЛИЦА РЕЗУЛЬТАТОВ ПО АНТЕННЕ
  ui->tableWidget_5->setColumnCount(6);
  ui->tableWidget_5->setRowCount(5 );
  ui->tableWidget_5->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_5->verticalHeader()->setVisible(true);

  QStringList lst5;
  lst5<<"X"<<"Y"<<"Z"<<"Bet"<<"Eps"<<"Alf";
  ui->tableWidget_5->setHorizontalHeaderLabels(lst5);
  QStringList lst6;
  lst6<< "Оценка,(м,гр)"<<"Ошибка,(м,мрад)"<<"СКО,(м,мрад)"<<"Сист,(м,мрад)"<<"Сумм,(м,мрад)";
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
  ui->tableWidget_7->setColumnCount(3);
  ui->tableWidget_7->setRowCount(5 );
  ui->tableWidget_7->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_7->verticalHeader()->setVisible(true);

  QStringList lst7;
  lst7<<"X"<<"Y"<<"Z";
  ui->tableWidget_7->setHorizontalHeaderLabels(lst7);
  QStringList lst8;
  lst8<< "Оценка"<<"Ошибка"<<"СКО"<<"Сист"<<"Сумм";
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
  marrGPSTruePosParams[0] = 0.27 ;
  marrGPSTruePosParams[1] = 15.66;
  marrGPSTruePosParams[2] = 17.83;

  marrGPSAprioriPosParams[0] = 0.27 ;
  marrGPSAprioriPosParams[1] = 15.66;
  marrGPSAprioriPosParams[2] = 17.83;


  mbtableWidget_8Init = false;
  for (int i=0; i<  ui->tableWidget_8->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget_8->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = nullptr;
                     ptwi0 = new QTableWidgetItem(QString::number(0.));
                     ui->tableWidget_8->setItem(i,j,ptwi0);
                 }


  for (int j =0; j < ui->tableWidget_8->rowCount(); ++j)
  {
   ui->tableWidget_8->item(j,0)->setText(QString::number(marrGPSTruePosParams[j]));
   ui->tableWidget_8->item(j,1)->setText(QString::number(marrGPSAprioriPosParams[j]));

  }

  mbtableWidget_8Init = true;
  /// !

  mType_of_000 = VAR0;
  ui->comboBox_7->setCurrentIndex(0);

  ui->pushButton->hide();

 // ui->pushButton_3->setVisible(false);
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

    mQstrOutPutFold = this->ui->lineEdit->text();
    mQstrOutPutFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[mQstrOutPutFold.size()] = 0;


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
   mvalProfMean = ui->doubleSpinBox_35->value();
   mvalProfSigma = ui->doubleSpinBox_36->value();
   for (int i = 0; i < mtblRealPrfl.mNumCols; ++i)
   {
      mtblRealPrfl.mparrVal[i] +=  getGauss(mvalProfMean, mvalProfSigma );
   }

   if(ui->radioButton->isChecked())
   {
       mtblRealPrfl = mtblRealPrfl.makeMiddleProfile();
   }


   for (int i = 0; i < 6; ++i)
   {
       if (i <3)
       {
           marrTruePosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();
           marrAprioriPosParams[i] =  ui->tableWidget_4->item(2,i)->text().toDouble();
       }
       else
       {
           marrTruePosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;
           marrAprioriPosParams[i] =  ui->tableWidget_4->item(2,i)->text().toDouble()* M_PI/ 180.;
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

        marrSinsSyst[0] = ui->doubleSpinBox_37->value()/ 1000.;
        marrSinsSyst[1] = ui->doubleSpinBox_38->value()/ 1000.;
        marrSinsSyst[2] = ui->doubleSpinBox_38->value()/ 1000.;

       // СКЗ угловой ошибки СИНС по скорости углов (углов качек)
         double valSig_d_po_dt_Sins = 1.16/1000.;
         double valMaxSig_Q      =    ui->doubleSpinBox_40->value()/ 1000.;
         double valMaxSig_Psi    =      ui->doubleSpinBox_39->value()/ 1000.;
         double valMaxSig_Tet    =    ui->doubleSpinBox_39->value()/ 1000.;
         double valMaxSig_dQdt   =      valSig_d_po_dt_Sins ;
         double valMaxSig_dPsidt =      valSig_d_po_dt_Sins ;
         double valMaxSig_dTetdt =      valSig_d_po_dt_Sins ;
         double valK1         = 0.01 ;
         double valSigV  = 0.1;
         double valSigH       =      0.1 ;
         double valMaxSig_H =     0.1 ;
         double valMaxSig_VH =     0.05 ;

     /*    double valMaxSig_Q      =     0.;//0.001 ; //0.000582;
           double valMaxSig_Psi    =    0;;//  valSigSins ; //0.00145;
           double valMaxSig_Tet    = 0.;//   valSigSins ; // 0.00145;
           double valMaxSig_dQdt   =  0.;//    valSig_d_po_dt_Sins ;
           double valMaxSig_dPsidt =  0.;//   valSig_d_po_dt_Sins ;
           double valMaxSig_dTetdt =  0.;//    valSig_d_po_dt_Sins ;
           double valK1         = 0.;//0.01 ;
           double valSigV  = 0.;//ui->doubleSpinBox_11->value();
           double valSigH       = 0.;//     0.1 ;
           double valMaxSig_H = 0.;//    0.1 ;
           double valMaxSig_VH = 0.;//    0.05 ;
           */

/*
       double valSig_d_po_dt_Sins = 0.;// ui->doubleSpinBox_10->value();
       double valMaxSig_Q      = 0.;//     valSigSins ; //0.000582;
       double valMaxSig_Psi    = 0.;//      valSigSins ; //0.00145;
       double valMaxSig_Tet    = 0.;//    valSigSins ; // 0.00145;
       double valMaxSig_dQdt   = 0.;//      valSig_d_po_dt_Sins ;
       double valMaxSig_dPsidt = 0.;//      valSig_d_po_dt_Sins ;
       double valMaxSig_dTetdt = 0.;//      valSig_d_po_dt_Sins ;
       double valK1         = 0.;// 0.01 ;
       double valSigV  = 0.;// ui->doubleSpinBox_11->value();
       double valSigH       = 0.;//      0.1 ;
       double valMaxSig_H = 0.;//     0.1 ;
       double valMaxSig_VH = 0.;//     0.05 ;
       */

         mT0 = 0.;
        mTimeTempSins = 1./ 50.;
       mSins = QPeaceSins (mEnvironment,valMaxSig_Q, valMaxSig_Psi, valMaxSig_Tet
                            ,valMaxSig_dQdt, valMaxSig_dPsidt,valMaxSig_dTetdt
                            ,valMaxSig_H,valMaxSig_VH,valK1
                            ,valSigV, mT0,mTimeTempSins,marrSinsSyst );
       ///
            // 3.4  параметры корабля
          // константы
       double valVesselWidth = 40; // ширина(м)
       double valVesselLength = 200.; // длина(м)

       // параметры движения
       double valQ0 = 0. ; // генеральный курс

       double valVVess =  ui->doubleSpinBox_15->value();
       mvalMaxQ =    3./180.*M_PI; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
       double valT_Q = 18.; // период рыскания
       mvalMaxPsi =      3./180.*M_PI;// максимальный угол килевой качки(амплитуда)
       double valT_Psi =  12; // период килевой качки
       mvalMaxTet =      12./180.*M_PI; //максимальный угол боротовой качки(амплитуда)
       double valT_Tet =  6.; // период бортовой качки
       double valMaxVert =     1. ;


       //

               switch (ui->comboBox_5->currentIndex())
               {
               case 0:
                   mHidroRLS = QHidroRLS((ui->doubleSpinBox_8->value()) * 0.001,0.04);
                   mpHidroRLS = &mHidroRLS;                  
                   mTypeOfSolverTask = LBL;
                   break;


               case 1:
                   mRls_Usbl2D = QRls_Usbl2D((ui->doubleSpinBox_8->value()) * 0.001,0.04
                                 ,(ui->doubleSpinBox_7->value()) * 0.001);
                   mpHidroRLS = &mRls_Usbl2D;
                   mTypeOfSolverTask = USBL_2D;
                   break;


               case 2:
                   mRls_Usbl3D = QRls_Usbl3D((ui->doubleSpinBox_8->value()) * 0.001,0.04
                     ,(ui->doubleSpinBox_7->value()) * 0.001,(ui->doubleSpinBox_13->value()) * 0.001);
                   mpHidroRLS = &mRls_Usbl3D;
                   mTypeOfSolverTask = USBL_3D;
                   break;
                default:
                   break;
               }

       mPlatform = QPlatform(mpHidroRLS, marrTruePosParams);
       //
       switch(ui->comboBox_4->currentIndex())
       {
       case 0:
           mTypeOfVessTraj = PNT_6;
           mQuantMeas = 6;
           break;
       case 1:
           mTypeOfVessTraj = PNT_7;
           mQuantMeas = 7;
           break;

       case 2:
           mTypeOfVessTraj = QDRT;
           mQuantMeas = (int)(ui->doubleSpinBox_14->value() +0.5);

           break;
       case 3:
           mTypeOfVessTraj = ZIG_ZAG;
           mQuantMeas = (int)(ui->doubleSpinBox_14->value() +0.5);

           break;

       case 4:
           mTypeOfVessTraj = DIAMETRS;
           mQuantMeas = (int)(ui->doubleSpinBox_14->value() +0.5);

           break;

       case 5:
           mTypeOfVessTraj = LINE6;
           mQuantMeas = (int)(ui->doubleSpinBox_14->value() +0.5);

           break;


       default:
           break;
       }
       //!
       // GPS
       read_tableWidget_8(marrGPSTruePosParams,marrGPSAprioriPosParams);
        mGPS = QGps(ui->doubleSpinBox_31->value(),marrGPSTruePosParams);

       //!

       mVess =  QPeaceVess (mEnvironment,valVesselWidth,valVesselLength
             ,mvalMaxQ ,valT_Q,mvalMaxPsi,valT_Psi ,mvalMaxTet
             ,valT_Tet,valMaxVert, valQ0,valVVess
             , mSins,  mPlatform, mT0,mGPS);
//!

     mTObrabotki = ui->doubleSpinBox_34->value() ;

     read_tableWidget_6();

     mImitMod = QImitMod(mVess, mtblRealPrfl
                         ,marrTrueHeadLightPos, marrAprioriHeadLightPos
                         , mTypeOfVessTraj, mTObrabotki,marrGPSAprioriPosParams);


     mMaxQuantIter_LBL = int(ui->doubleSpinBox_32->value()+ 0.1);
     mLBLcoeff= ui->doubleSpinBox_19->value();
     switch (ui->comboBox_8->currentIndex())
     {
     case 0:
     mALGOR_TYPE = PEREBOR;
         break;

     case 1:
     mALGOR_TYPE = NUTON;
     mMaxQuantIter_USBL = int(ui->doubleSpinBox_16->value()+ 0.1);
     mUSBLcoeff= ui->doubleSpinBox_17->value();
         break;

     default:
         break;
     }
  mBeaconDeepthToleran = ui->doubleSpinBox_42->value();
  mGpsToleranInstall = ui->doubleSpinBox_41->value();
}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\AKIN\\ОТЧЕТ_ПОЗИЦИОНИРОВАНИЕ\\REZ_REPORT");
    this->ui->lineEdit->setText(strFold);


}
//--------------------------------------------
void MainWindow::on_pushButton_3_clicked()
{
    inputData();
    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);


    QBigMeasure *parrMeas = (QBigMeasure *)malloc(mQuantMeas * sizeof(QBigMeasure) );

    mImitMod.createMeasures(mQuantMeas, parrMeas);

    QString qstrOutPutData =  mQstrOutPutFold +"\\inputDataExchangeFile.acd";

    writeExchangeDataFile_(qstrOutPutData,
                  mQuantMeas// к-во измерений
                  ,parrMeas//
                  , marrAprioriPosParams
                  ,marrGPSAprioriPosParams
                  ,marrAprioriHeadLightPos
                  ,mTObrabotki
                  );


    TURPointXY *pnts = new TURPointXY[mQuantMeas];
    for (int i =0; i< mQuantMeas; ++i)
    {
       pnts[i].X = parrMeas[i].marrSVessWaveZv[0];
       pnts[i].Y = parrMeas[i].marrSVessWaveZv[1];
       int uu=0;
    }
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TURPointXY pnt00(marrAprioriHeadLightPos[0], marrAprioriHeadLightPos[1]);
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.,pnt00) ;


    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\PntsTraj.shp");
    pnts[0].WriteSetSHPFiles(wchAxesFileName0,pnts,mQuantMeas);
    delete []pnts;


    mLblSolver = QLblSolver(parrMeas, mQuantMeas    ,mtblEstPrfl
                            ,marrEstHeadLightPos[2], 0.015);

    // ОТЛАДКА

            QBigMeasure meas = mLblSolver.mVectBigMeasures.at(0);
            meas.mTotvZv = 93.17058938;
            meas.mTzaprZv = 91.7839599;
            meas.marrMuWaveZv[0] = 1.570796;
            meas.marrMuWaveZv[1] = 0. + 0.001;
            meas.marrMuWaveZv[2] = 0.;

            meas.marrMuZv[0] = 1.570796;
            meas.marrMuZv[1] = + 0.001;
            meas.marrMuZv[2] = 0.;

            meas.marrSVessWaveZv[0] = 26.7196;
            meas.marrSVessWaveZv[1] = 138.713;
            meas.marrSVessWaveZv[2] = 6.8422;

            meas.marrSVessZv[0] = 26.0211;
            meas.marrSVessZv[1] = 140.6387;
            meas.marrSVessZv[2] = 2.98;

            mLblSolver.mVectBigMeasures.replace(0,meas);
            double arrX0_[] = {5.8553, -7.6252,-5.4747,-0.1116, 0.0019};
            double parrArr_df_po_dAngSins[15] = {0.};
            mLblSolver.calc_arrfi_and_arrdfi_po_dAngSins(arrX0_
                         ,0,parrArr_df_po_dAngSins);
            // ОТЛАДКА  !

    double arrX0[9] ={0.},arrXRez [9] = {0.}, valNeviaz0 = -1.;
    arrX0[0]= marrAprioriPosParams[0];
    arrX0[1]= marrAprioriPosParams[1];
    arrX0[2]= marrAprioriPosParams[2];

    arrX0[3]= marrAprioriHeadLightPos[0];
    arrX0[4]= marrAprioriHeadLightPos[1];
   // arrX0[5]= marrAprioriHeadLightPos[2];

    arrX0[5]= marrAprioriPosParams[3];
    arrX0[6]= marrAprioriPosParams[4];
    arrX0[7]= marrAprioriPosParams[5];

 memcpy(arrXRez, arrX0, 8 * sizeof(double));
 //if(!mLblSolver.estimateParams(arrX0,arrXRez, valNeviaz0))
  //  {
  //      int uu =0;
   // }
 // выгрузка таблицы маяка
 for (int j =0; j < 2; ++j)
 {
  int ia = ((int)(arrXRez[3 + j] * 1000.)) ;
  double a = ((double)ia)/ 1000.;
  ui->tableWidget_7->item(0,j)->setText(QString::number(a));

   ia = ((int)((arrXRez[3 + j] - marrTrueHeadLightPos[j]) * 1000.)) ;
   a = ((double)ia)/ 1000.;
  ui->tableWidget_7->item(1,j)->setText(QString::number(a));
  ui->tableWidget_7->item(2,j)->setText(QString("???"));
 }

 int ia = ((int)(marrAprioriHeadLightPos[2] * 1000.)) ;
 double a = ((double)ia)/ 1000.;
 ui->tableWidget_7->item(0,2)->setText(QString::number(a));

 ia = ((int)((marrAprioriHeadLightPos[2] - marrTrueHeadLightPos[2]) * 1000.)) ;
 a = ((double)ia)/ 1000.;
 ui->tableWidget_7->item(1,2)->setText(QString::number(a));
 ui->tableWidget_7->item(2,2)->setText(QString("???"));
 // !

 // выгрузка таблицы антенны
 for (int j =0; j < 3; ++j)
 {
     int ia = ((int)(arrXRez[ j] * 1000.)) ;
     double a = ((double)ia)/ 1000.;
  ui->tableWidget_5->item(0,j)->setText(QString::number(a));

   ia = ((int)((arrXRez[ j] - marrTruePosParams[j]) * 1000.)) ;
   a = ((double)ia)/ 1000.;
  ui->tableWidget_5->item(1,j)->setText(QString::number(a));

  ui->tableWidget_5->item(2,j)->setText(QString("???"));
 }

    bool bContinueCalculations = false;
    double arrSBeacon_GSK[3] = {0.};
    arrSBeacon_GSK[0] = arrXRez[3];
    arrSBeacon_GSK[1] = arrXRez[4];
    arrSBeacon_GSK[2] = marrTrueHeadLightPos[2];
    switch(ui->comboBox_5->currentIndex())
        {

        case 1:
            mSolver2D_Angs = QSolver2D_Angs(parrMeas, mQuantMeas    ,mtblEstPrfl
                 ,marrTrueHeadLightPos[2],0.0001,arrXRez,arrSBeacon_GSK);

            mpPosSolver = &mSolver2D_Angs;
            bContinueCalculations = true;
            break;

        case 2:
            mSolver3D_Angs = QSolver3D_Angs(parrMeas, mQuantMeas    ,mtblEstPrfl
                       ,marrTrueHeadLightPos[2],0.0001,arrXRez,arrSBeacon_GSK);
            mpPosSolver = &mSolver3D_Angs;
            bContinueCalculations = true;

            break;

        default:

            break;
        }


    if (bContinueCalculations)
    {
       mpPosSolver->estimateParams(&arrX0[5],&arrXRez[5], valNeviaz0) ;
    }

    free (parrMeas);

   if(mTypeOfSolverTask != LBL)
   {
       for (int j =0; j < 3; ++j)
       {
           int ia = ((int)(arrXRez[5 + j] * 180./M_PI * 100.)) ;
           double a = ((double)ia)/ 100.;
        ui->tableWidget_5->item(0,3 +j)->setText(QString::number(a));


        ia = ((int)((arrXRez[5 + j] -marrTruePosParams[3 + j]) *1000. * 100.)) ;
         a = ((double)ia)/ 100.;
        ui->tableWidget_5->item(1,3 + j)->setText(QString::number(a));

        ui->tableWidget_5->item(2,3 + j)->setText(QString("???"));
       }

   }
    // !
}

void MainWindow::on_comboBox_4_currentIndexChanged(int index)
{
 switch(index)
 {
    case 0:
    case 1:    
     ui->doubleSpinBox_14->setValue(index+6);
     ui->doubleSpinBox_14->setReadOnly(true);
     break;
    case 2:
    case 3:
     ui->doubleSpinBox_14->setReadOnly(false);
     break;
 default:
     break;

 }
}

void MainWindow::on_tableWidget_4_cellChanged(int row, int column)
{
   /* if (!mbtableWidget_4Init)
    {
        return;
    }
    double arr[18] = {0.};
    read_tableWidget_4(arr);

    // перевод 0 строки в радианы
    for (int i =3; i < 6; ++i)
    {
      arr[i] *= M_PI/ 180.;
    }
    // формирование 3-ей строки
    for (int i =0; i < 6; ++i)
    {
      if (i < 3)
      {
      arr[12 + i ] = arr[ i ] + arr[6 + i ];
      }
      else
      {
       arr[12 + i ] = arr[ i ] + arr[6 + i ]/1000.;
      }
    }
    // перевод 0 и 3 строк в градусы
    // перевод 0 строки в радианы
    for (int i =3; i < 6; ++i)
    {
      arr[i] =  arr[i] * 180./M_PI;
      arr[12 + i] =  arr[12 + i] * 180./M_PI;
    }

   write_tableWidget_4(arr);
    */
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
    /*if (!mbtableWidget_6Init)
    {
        return;
    }
    double arr[9] = {0.};
    read_tableWidget_6(arr);


    // формирование 3-ей строки
    for (int i =0; i < ui->tableWidget_6->columnCount(); ++i)
    {
        arr[(ui->tableWidget_6->columnCount() *2) + i ] = arr[ i ] + arr[ui->tableWidget_6->columnCount() + i ];

    }

   write_tableWidget_6(arr);
*/
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
void MainWindow::read_tableWidget_8(double *arr1,double *arr2)
{
     for (int i=0; i < ui->tableWidget_8->rowCount(); ++i)
    {
          arr1 [i] =  ui->tableWidget_8->item(i,0)->text().toDouble();
          arr2 [i] =  ui->tableWidget_8->item(i,1)->text().toDouble();
    }
}
//-------------------------------------------------
void MainWindow::read_tableWidget_6()
{
     for (int i=0; i < ui->tableWidget_6->columnCount(); ++i)
    {
          marrTrueHeadLightPos[i] =  ui->tableWidget_6->item(0,i)->text().toDouble();
          marrAprioriHeadLightPos[i] =  ui->tableWidget_6->item(2,i)->text().toDouble();
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
    QString str0 = QCoreApplication::applicationDirPath() ;
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости",str0,"*.000");
   // QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости","..\\POSITION_IMIT_MODEL_RAZD_OBRAB_v0" ,"*.000");
    this->ui->lineEdit_3->setText(strFold);

    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPrflFIle, array);
    QSubWaterBeam::createProfileTbl(mpwchPrflFIle, &mtblEstPrfl, mType_of_000);
    mDeepthMax = mtblEstPrfl.mparrArg[mtblEstPrfl.mNumCols - 2];
   // ui->doubleSpinBox_32->setValue(mDeepthMax);

    mZonaR = calcXMCrit( 2., mDeepthMax, mtblEstPrfl);

    ui->doubleSpinBox_33->setValue(mZonaR);

}
//------------------------------------------------------
void MainWindow::on_tableWidget_7_itemActivated(QTableWidgetItem *item)
{
   read_tableWidget_7(marrMayjak_TblData);
}
//-------------------------------------------------
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

void MainWindow::on_tableWidget_7_itemSelectionChanged()
{
   write_tableWidget_7(marrMayjak_TblData);
}
//---------------------------------
//--------------------------------

void MainWindow::on_tableWidget_5_itemActivated(QTableWidgetItem *item)
{
  // read_tableWidget_5(marrGRLS_TblData);
}
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
void MainWindow::read_tableWidget_8(double *arr)
{
for (int i=0; i < ui->tableWidget_8->rowCount(); ++i)
      for (int j=0; j < ui->tableWidget_8->columnCount(); ++j)
      {
       arr [ i * (ui->tableWidget_8->columnCount()) + j] =  ui->tableWidget_8->item(i,j)->text().toDouble();
      }
}
//-------------------------------------------------
void MainWindow::write_tableWidget_8(double *arr)
{

   int ncols = ui->tableWidget_8->columnCount();
    for (int i=0; i < ui->tableWidget_8->rowCount(); ++i)
         for (int j=0; j < ncols; ++j)
         {
             ui->tableWidget_8->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

         }
}

void MainWindow::on_tableWidget_5_itemSelectionChanged()
{
  // write_tableWidget_5(marrGRLS_TblData);
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

// обработчик нажатия кнопки Start (запуск задачи)
void MainWindow::on_pushButtonStartTask_clicked()
{
    inputData();
    ui->pushButton->show();
    ui->labelTaskInfo1->setText("Имитация...");
    ui->labelTaskInfo2->setText("");
    ui->labelTaskInfo3->setText("");
    ui->labelTaskInfo4->setText("");

    // создаем поток и привязываем сигналы к слотам
    myThread = new MyThread;
    QObject::connect(myThread, SIGNAL(finished()), this, SLOT(onTaskFinished()));

    QObject::connect(myThread, SIGNAL(progress(int , QVector<double>, QVector<double>
                           , QVector<double> , int  ,bool , QVector<double>))
            , this, SLOT(onTaskProgress(int , QVector<double>, QVector<double>
                                        , QVector<double> , int  ,bool , QVector<double>)));
    QObject::connect(myThread, SIGNAL(calcfinished(int)), this, SLOT(onCalcFinished(int)));


    QObject::connect(myThread, SIGNAL(LBLSucceeded( bool
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ))
            , this, SLOT(onLBLTaskSucceeded(              bool
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          )));

    QObject::connect(myThread, SIGNAL(USBLSucceeded( bool
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ))
            , this, SLOT(onUSBLTaskSucceeded(              bool
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          ,QVector<double>
                                                          )));
    qDebug() << "MainWindow::on_pushButtonStartTask_clicked";

    // очистка интерфейса
    ui->pushButtonStartTask->setEnabled(false);
    ui->pushButtonStopTask->setEnabled(true);
    ui->progressBar->setValue(0);
    ui->progressBar->setRange(0, mMaxQuantIter_LBL);

    int iiii= ui->progressBar->maximum();





    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);


    QBigMeasure *parrMeas = (QBigMeasure *)malloc(mQuantMeas * sizeof(QBigMeasure) );

     mImitMod.createMeasures(mQuantMeas, parrMeas);
     ui->labelTaskInfo1->setText("Этап 1: Вычисление офсетов...");

    QString qstrOutPutData =  mQstrOutPutFold +"\\inputDataExchangeFile.acd";

    writeExchangeDataFile_(qstrOutPutData,
                  mQuantMeas// к-во измерений
                  ,parrMeas//
                  , marrAprioriPosParams
                  ,marrGPSAprioriPosParams
                  ,marrAprioriHeadLightPos
                  ,mTObrabotki
                  );

/*
    // ОТЛАДКА
    QVector<QBigMeasure> vctBigMesure1(mQuantMeas);
    double arrAprioriPosParams[6]
            , arrGPSAprioriPosParams[3], arrAprioriHeadLightPos[3], TObrabotki;
    readFile_acd(qstrOutPutData,arrAprioriPosParams
                      , arrGPSAprioriPosParams, arrAprioriHeadLightPos, TObrabotki
                      , &vctBigMesure1);
    int QuantMeas1 = vctBigMesure1.size();
    QBigMeasure *parrMeas1 = (QBigMeasure *)malloc(QuantMeas1 * sizeof(QBigMeasure) );
    for (int i = 0; i < QuantMeas1; ++i)
    {
      parrMeas1[i] = vctBigMesure1.at(i);
      parrMeas1[i].mSig_e = ui->doubleSpinBox_13->value() * 0.001;
      parrMeas1[i].mSig_q = ui->doubleSpinBox_7->value() * 0.001;
      parrMeas1[i].mSig_t = ui->doubleSpinBox_8->value() * 0.001;
      vctBigMesure1.replace(i, parrMeas1[i]);
    }
    QString qstrOutPutData1 =  mQstrOutPutFold +"\\inputDataExchangeFile1.acd";
    writeExchangeDataFile_(qstrOutPutData1,
                  QuantMeas1// к-во измерений
                  ,parrMeas1//
                  , arrAprioriPosParams
                  ,arrGPSAprioriPosParams
                  ,arrAprioriHeadLightPos
                  ,TObrabotki
                  );

  for (int i =0; i < QuantMeas1; ++i)
  {
      const double VAlTolerDist = 0.01
                              , VAlTolerSinsAngs = 0.001,   VAlTolerT = 0.00001,  VAlTolerMeasAngs = 0.0001;
     if(! QBigMeasure::isEqual(parrMeas1[i], parrMeas[i],  VAlTolerDist
                        , VAlTolerSinsAngs,   VAlTolerT,  VAlTolerMeasAngs ))
     {
         int uu = 0;
     }
  }

// ! ОТЛАДКА
*/
    mNumOfTask =1;

    TURPointXY *pnts = new TURPointXY[mQuantMeas];
    for (int i =0; i< mQuantMeas; ++i)
    {
       pnts[i].X = parrMeas[i].marrSVessWaveZv[0];
       pnts[i].Y = parrMeas[i].marrSVessWaveZv[1];
       int uu=0;
    }
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TURPointXY pnt00(marrAprioriHeadLightPos[0], marrAprioriHeadLightPos[1]);
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.,pnt00) ;


    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\PntsTraj.shp");
    pnts[0].WriteSetSHPFiles(wchAxesFileName0,pnts,mQuantMeas);
    delete []pnts;

    QVector <QBigMeasure> vctMeas(mQuantMeas);
    for (int i =0; i < mQuantMeas; ++i)
    {
     vctMeas.replace(i, parrMeas[i]);
    }




    // задаём параметры выполнения задачи
    // здесь можно передать массив измерений, профиль скорости звука и пр. параметры
    double arrX0[8] ={0.};
    arrX0[0]= marrAprioriPosParams[0];
    arrX0[1]= marrAprioriPosParams[1];
    arrX0[2]= marrAprioriPosParams[2];

    arrX0[3]= marrAprioriHeadLightPos[0];
    arrX0[4]= marrAprioriHeadLightPos[1];
   // arrX0[5]= marrAprioriHeadLightPos[2];

    arrX0[5]= marrAprioriPosParams[3];
    arrX0[6]= marrAprioriPosParams[4];
    arrX0[7]= marrAprioriPosParams[5];
    QVector <double>vctX0(8);


    for (int i = 0; i < 8; ++i)
    {
      vctX0.replace(i,arrX0[i]);

    }

    myThread->setParams(vctMeas, mTypeOfSolverTask,vctX0, mtblEstPrfl
                        ,marrGPSAprioriPosParams, marrAprioriHeadLightPos[2]
                      ,mMaxQuantIter_LBL, mMaxQuantIter_USBL
                         ,mGape,mAngStep,  mALGOR_TYPE, mLBLcoeff,mUSBLcoeff);


   mNumOfTask = 1;

    // запускаем поток (при этом автоматически запустится функция MyThread::run())
    myThread->start();
}

// обработчик нажатия кнопки Stop (принудительная остановка задачи пользователем)
void MainWindow::on_pushButtonStopTask_clicked()
{
    qDebug() << "MainWindow::on_pushButtonStopTask_clicked";
    myThread->userBreak(); // устанавливаем в потоке флаг завершения
}

// обработчик сигнала о текущем прогрессе выполнения задачи
void MainWindow::onTaskProgress(int step, QVector<double>vctNevSquare, QVector<double> vctMean
                                , QVector<double> vctDisp, int numPart ,bool bEndPart, QVector<double>vctX0)
{
  //int iiii= ui->progressBar->maximum();
    // вывод в интерфейс

    for(int i =0;i < 8; ++i)
    {
      marrXRez[i] = vctX0.at(i);
    }
    if (bEndPart)
    {
       if(numPart == 1)
       {

           // выгрузка таблицы маяка
           for (int j =0; j < 2; ++j)
           {
            int ia = ((int)(marrXRez[3 + j] * 1000.)) ;
            double a = ((double)ia)/ 1000.;
            ui->tableWidget_7->item(0,j)->setText(QString::number(a));

             ia = ((int)((marrXRez[3 + j] - marrTrueHeadLightPos[j]) * 1000.)) ;
             a = ((double)ia)/ 1000.;
            ui->tableWidget_7->item(1,j)->setText(QString::number(a));
            ui->tableWidget_7->item(2,j)->setText(QString("???"));
           }

           int ia = ((int)(marrAprioriHeadLightPos[2] * 1000.)) ;
           double a = ((double)ia)/ 1000.;
           ui->tableWidget_7->item(0,2)->setText(QString::number(a));

           ia = ((int)((marrAprioriHeadLightPos[2] - marrTrueHeadLightPos[2]) * 1000.)) ;
           a = ((double)ia)/ 1000.;
           ui->tableWidget_7->item(1,2)->setText(QString::number(a));
           ui->tableWidget_7->item(2,2)->setText(QString("???"));
           // !

           // выгрузка таблицы антенны
           for (int j =0; j < 3; ++j)
           {
               int ia = ((int)(marrXRez[ j] * 1000.)) ;
               double a = ((double)ia)/ 1000.;
            ui->tableWidget_5->item(0,j)->setText(QString::number(a));

             ia = ((int)((marrXRez[ j] - marrTruePosParams[j]) * 1000.)) ;
             a = ((double)ia)/ 1000.;
            ui->tableWidget_5->item(1,j)->setText(QString::number(a));

            ui->tableWidget_5->item(2,j)->setText(QString("???"));
           }
       }
       else
       {
          for (int j =0; j < 3; ++j)
           {
               int ia = ((int)(marrXRez[5 + j] * 180./M_PI * 100.)) ;
               double a = ((double)ia)/ 100.;
            ui->tableWidget_5->item(0,3 +j)->setText(QString::number(a));


            ia = ((int)((marrXRez[5 + j] -marrTruePosParams[3 + j]) *1000. * 100.)) ;
             a = ((double)ia)/ 100.;
            ui->tableWidget_5->item(1,3 + j)->setText(QString::number(a));

            ui->tableWidget_5->item(2,3 + j)->setText(QString("???"));
           }
       }
    }

    if (mNumOfTask != numPart)
    {
       mNumOfTask++;
       ui->labelTaskInfo1->setText("Этап 2: вычисление углов...");

       if ( mALGOR_TYPE  == PEREBOR)
       {
           double in_v_storonu = mGape/mAngStep;
           int in = 2* in_v_storonu +1;
           int iRange = in * in * in;
           ui->progressBar->setRange(0, iRange);
       }

       if ( mALGOR_TYPE  == NUTON)
       {
           ui->progressBar->setRange(0, mMaxQuantIter_USBL);
       }

       ui->progressBar->setValue(0);
    }
    ui->progressBar->setValue(step);
    if (mNumOfTask == 1)
    {
        ui->labelTaskInfo2->setText(QString("t(мс): N= %1, СКО=  %2, Сред= %3, СКЗ= %4").arg(step).arg(sqrt(vctNevSquare.at(0)) *1000., 0, 'f', 3).arg(vctMean.at(0) *1000., 0, 'f', 3).arg(sqrt(vctDisp.at(0)) *1000., 0, 'f', 3));

    }
    else
    {
        ui->labelTaskInfo3-> setText(QString("N= %1; q(мр):  СКО=  %2, Сред= %3, СКЗ= %4").arg(step).arg(sqrt(vctNevSquare.at(0)) *1000., 0, 'f', 3).arg(vctMean.at(0) *1000., 0, 'f', 3).arg(sqrt(vctDisp.at(0)) *1000., 0, 'f', 3));
        if(USBL_3D == mTypeOfSolverTask)
        {
         ui->labelTaskInfo4->setText(QString("          e(мр):  СКО=  %1, Сред= %2, СКЗ= %3").arg(sqrt(vctNevSquare.at(1)) *1000., 0, 'f', 3).arg(vctMean.at(1) *1000., 0, 'f', 3).arg(sqrt(vctDisp.at(1)) *1000., 0, 'f', 3));
        }
    }
}
//-------------------------------------------------------------------------
// обработчик сигнала об успешном выполнении задачи LBL
void MainWindow::onLBLTaskSucceeded(bool bbegin
                                    ,QVector<double>vct_arr_dX_po_dSgps
                                    ,QVector<double>vct_arr_dX_po_dHBeacon
                                    ,QVector<double>vct_arr_Kgps
                                    ,QVector<double>vct_arr_Kt
                                    ,QVector<double>vct_arr_dX_po_AngSins
                                    ,QVector<double>vct_arr_KsinsQ_per_1_rad
                                    ,QVector<double>vct_arr_Kgps_Psi_Tetta_per_1_rad
                                    ,QVector<double>vct_arr_dX_po_dProfile
                                    ,QVector<double>vct_arr_KProfile_per_1_m
                                    )
{

    if (bbegin)
    {
       ui->labelTaskInfo1->setText("Этап 1: Расчет точности...");
       return;
    }


    // ОТЛАДКА

    double arr_dX_po_dSgps[15]= {0.};


    double arr_dX_po_dHBeacon[5] = {0.};


    double arr_Kgps[25] = {0.};


    double arr_Kt[25] = {0.};


    double arr_dX_po_AngSins[5*3] = {0.};


    double arr_KsinsQ_per_1_rad[5*5] = {0.};

    double arr_Kgps_Psi_Tetta_per_1_rad[25] ={0.};


    double arr_dX_po_dProfile[5] = {0.};


    double arr_KProfile_per_1_m[25] ={0.};

    QPosSolver:: createDblArrayFromDblVector(arr_dX_po_dSgps,vct_arr_dX_po_dSgps);//

    QPosSolver:: createDblArrayFromDblVector(arr_dX_po_dHBeacon,vct_arr_dX_po_dHBeacon);//

   QPosSolver:: createDblArrayFromDblVector(arr_Kgps,vct_arr_Kgps);//+

   QPosSolver:: createDblArrayFromDblVector(arr_Kt,vct_arr_Kt);//

   QPosSolver:: createDblArrayFromDblVector(arr_dX_po_AngSins,vct_arr_dX_po_AngSins);//

   QPosSolver:: createDblArrayFromDblVector(arr_KsinsQ_per_1_rad,vct_arr_KsinsQ_per_1_rad);//

   QPosSolver:: createDblArrayFromDblVector(arr_Kgps_Psi_Tetta_per_1_rad,vct_arr_Kgps_Psi_Tetta_per_1_rad);//

   QPosSolver:: createDblArrayFromDblVector(arr_dX_po_dProfile,vct_arr_dX_po_dProfile);

   // QPosSolver:: createDblArrayFromDblVector(arr_KProfile_per_1_m,vct_arr_KProfile_per_1_m);//

    // вычисление коррел матрицы
    double arrK[25] = {0.};
    MatrxMultScalar(arr_Kgps, 5, 5, mGPS.mSigXY * mGPS.mSigXY,arr_Kgps);
    MtrxSumMatrx(arrK, arr_Kgps, 5, 5, arrK) ;

    MatrxMultScalar(arr_Kt, 5, 5, mpHidroRLS->mSigT *mpHidroRLS->mSigT,arr_Kt);
    MtrxSumMatrx(arrK,arr_Kt, 5, 5, arrK) ;

    MatrxMultScalar(arr_KsinsQ_per_1_rad, 5, 5, mSins.mSig_Q *mSins.mSig_Q,arr_KsinsQ_per_1_rad);
    MtrxSumMatrx(arrK,arr_KsinsQ_per_1_rad, 5, 5, arrK) ;

    MatrxMultScalar(arr_Kgps_Psi_Tetta_per_1_rad, 5, 5, mSins.mSig_Psi *mSins.mSig_Psi,arr_Kgps_Psi_Tetta_per_1_rad);
    MtrxSumMatrx(arrK,arr_Kgps_Psi_Tetta_per_1_rad, 5, 5, arrK) ;

    MatrxMultScalar(arr_KProfile_per_1_m, 5, 5, mvalProfSigma *mvalProfSigma,arr_KProfile_per_1_m);
    MtrxSumMatrx(arrK,arr_KProfile_per_1_m, 5, 5, arrK) ;

    double arrSKZ[5] ={0.};
    for (int i =0; i < 5;++i)
    {
       arrSKZ[i] =sqrt(arrK[ 5 * i +i]);
    }
    ///
   //систематическая ошибка
    double arrSyst[5] ={0.};
    for (int i = 0; i < 5; ++i)
    {
        double temp0 = arr_dX_po_dSgps[i * 3]*arr_dX_po_dSgps[i * 3]
                +arr_dX_po_dSgps[i * 3+1]*arr_dX_po_dSgps[i * 3+1]
                +arr_dX_po_dSgps[i * 3+2]*arr_dX_po_dSgps[i * 3+2];
        temp0 *= mGpsToleranInstall*mGpsToleranInstall;

       double temp1 = arr_dX_po_AngSins[i * 3] *arr_dX_po_AngSins[i * 3]* marrSinsSyst[0]* marrSinsSyst[0]
               +arr_dX_po_AngSins[i * 3 +1] *arr_dX_po_AngSins[i * 3+1]* marrSinsSyst[1]* marrSinsSyst[1]
               +arr_dX_po_AngSins[i * 3 +2] *arr_dX_po_AngSins[i * 3+2]* marrSinsSyst[2]* marrSinsSyst[2];

       arrSyst[i] =  temp0
               + arr_dX_po_dHBeacon[i]*arr_dX_po_dHBeacon[i] * mBeaconDeepthToleran* mBeaconDeepthToleran
               + temp1
               + arr_dX_po_dProfile[i]*arr_dX_po_dProfile[i] *mvalProfMean *mvalProfMean;
       arrSyst[i] = sqrt(arrSyst[i]);
    }
    //

    // суммарная ошибка
    double arrSum[5] ={0.};
    for (int i = 0; i < 5; ++i)
    {
       arrSum[i] =sqrt(arrSyst[i] * arrSyst[i] +arrSKZ[i] * arrSKZ[i]);
    }

    //
    // выгрузка таблицы маяка
    for (int j =0; j < 2; ++j)
    {
     int ia = ((int)(arrSKZ[3 + j] * 1000.)) ;
     double a = ((double)ia)/ 1000.;
     ui->tableWidget_7->item(2,j)->setText(QString::number(a));

      ia = ((int)(arrSyst[3 + j] * 1000.)) ;
      a = ((double)ia)/ 1000.;
     ui->tableWidget_7->item(3,j)->setText(QString::number(a));

     ia = ((int)(arrSum[3 + j] * 1000.)) ;
     a = ((double)ia)/ 1000.;
    ui->tableWidget_7->item(4,j)->setText(QString::number(a));
    }
    ui->tableWidget_7->item(2,2)->setText(QString::number(0.));
    ui->tableWidget_7->item(3,2)->setText(QString::number(mBeaconDeepthToleran));;
    ui->tableWidget_7->item(4,2)->setText(QString::number(mBeaconDeepthToleran));




    // !

    // выгрузка таблицы антенны
    for (int j =0; j < 3; ++j)
    {
        int ia = ((int)(arrSKZ[ j] * 1000.)) ;
        double a = ((double)ia)/ 1000.;
     ui->tableWidget_5->item(2,j)->setText(QString::number(a));

     ia = ((int)(arrSyst[ j]  * 1000.)) ;
     a = ((double)ia)/ 1000.;
    ui->tableWidget_5->item(3,j)->setText(QString::number(a));

    ia = ((int)(arrSum[ j]  * 1000.)) ;
    a = ((double)ia)/ 1000.;
   ui->tableWidget_5->item(4,j)->setText(QString::number(a));


    }

}

//-------------------------------------------------------------------------
// обработчик сигнала об успешном выполнении задачи USBL
void MainWindow::onUSBLTaskSucceeded(bool bbegin
                                    ,QVector<double>vct_arr_dAngs_po_dSgps
                                    ,QVector<double>vct_arr_dAngs_po_dH //
                                    ,QVector<double>vct_arr_K_dqe
                                    ,QVector<double>vct_arr_dAngs_po_dAngSins
                                    ,QVector<double>vct_arr_dAngs_po_dProfile
                                    ,QVector<double>vct_arr_calc_AngKgps
                                    ,QVector<double>vct_arr_calc_AngKsinsQ
                                    ,QVector<double>vct_arr_calc_AngK_Psi_Tetta_per_1_rad
                                    ,QVector<double>vct_arr_calc_arr_AngKt_per_1sec
                                    )

{

    if (bbegin)
    {
       ui->labelTaskInfo1->setText("Этап 2: Расчет точности...");
       return;
    }

    double   arr_dAngs_po_dSgps[9  ] = {0.};
    double   arr_dAngs_po_dH[3  ] = {0.};
    double   arr_K_dqe[9  ] = {0.};
    double   arr_dAngs_po_dAngSins[9  ] = {0.};
    double   arr_dAngs_po_dProfile[3  ] = {0.};
    double   arr_calc_AngKgps[9  ] = {0.};
    double   arr_calc_AngKsinsQ[9  ] = {0.};
    double   arr_calc_AngK_Psi_Tetta_per_1_rad[9  ] = {0.};
    double   arr_calc_arr_AngKt_per_1sec[9  ] = {0.};

    QPosSolver:: createDblArrayFromDblVector(arr_dAngs_po_dSgps,vct_arr_dAngs_po_dSgps);//++

    QPosSolver:: createDblArrayFromDblVector(arr_dAngs_po_dH,vct_arr_dAngs_po_dH);//++

   QPosSolver:: createDblArrayFromDblVector(arr_K_dqe,vct_arr_K_dqe);//+

   QPosSolver:: createDblArrayFromDblVector(arr_dAngs_po_dAngSins,vct_arr_dAngs_po_dAngSins);//++

   QPosSolver:: createDblArrayFromDblVector(arr_dAngs_po_dProfile,vct_arr_dAngs_po_dProfile);//

   QPosSolver:: createDblArrayFromDblVector(arr_calc_AngKgps,vct_arr_calc_AngKgps);//++

   QPosSolver:: createDblArrayFromDblVector(arr_calc_AngKsinsQ,vct_arr_calc_AngKsinsQ);//++

   QPosSolver:: createDblArrayFromDblVector(arr_calc_AngK_Psi_Tetta_per_1_rad,vct_arr_calc_AngK_Psi_Tetta_per_1_rad);//++

   QPosSolver:: createDblArrayFromDblVector(arr_calc_arr_AngKt_per_1sec,vct_arr_calc_arr_AngKt_per_1sec );//++

    // вычисление коррел матрицы
    double arrK[9] = {0.};
    MatrxMultScalar(arr_calc_AngKgps, 3, 3, mGPS.mSigXY * mGPS.mSigXY,arr_calc_AngKgps);
    MtrxSumMatrx(arrK, arr_calc_AngKgps, 3, 3, arrK) ;

    MatrxMultScalar(arr_calc_arr_AngKt_per_1sec, 3, 3, mpHidroRLS->mSigT *mpHidroRLS->mSigT,arr_calc_arr_AngKt_per_1sec);
    MtrxSumMatrx(arrK,arr_calc_arr_AngKt_per_1sec, 3, 3, arrK) ;

    MatrxMultScalar(arr_calc_AngKsinsQ, 3, 3, mSins.mSig_Q *mSins.mSig_Q,arr_calc_AngKsinsQ);
    MtrxSumMatrx(arrK, arr_calc_AngKsinsQ, 3, 3, arrK) ;

    MatrxMultScalar(arr_calc_AngK_Psi_Tetta_per_1_rad, 3, 3, mSins.mSig_Psi *mSins.mSig_Psi, arr_calc_AngK_Psi_Tetta_per_1_rad);
    MtrxSumMatrx(arrK,arr_calc_AngK_Psi_Tetta_per_1_rad, 3, 3, arrK) ;



    MtrxSumMatrx(arrK,arr_K_dqe, 3, 3, arrK) ;

    double arrSKZ[3] ={0.};
    for (int i =0; i < 3;++i)
    {
       arrSKZ[i] =sqrt(arrK[ 3 * i +i]);
    }
    //
   //систематическая ошибка
    double arrSyst[3] ={0.};
    for (int i = 0; i < 3; ++i)
    {
        double temp0 = arr_dAngs_po_dSgps[i * 3]* arr_dAngs_po_dSgps[i * 3];
              temp0 += arr_dAngs_po_dSgps[i * 3+1] * arr_dAngs_po_dSgps[i * 3+1];
              temp0 += arr_dAngs_po_dSgps[i * 3+2] * arr_dAngs_po_dSgps[i * 3+2];
        temp0 = temp0 *mGpsToleranInstall*mGpsToleranInstall;

       double temp1 = arr_dAngs_po_dAngSins[i * 3] *arr_dAngs_po_dAngSins[i * 3]* marrSinsSyst[0]* marrSinsSyst[0];

              temp1 +=arr_dAngs_po_dAngSins[i * 3 +1] *arr_dAngs_po_dAngSins[i * 3+1]* marrSinsSyst[1]* marrSinsSyst[1];

              temp1 +=arr_dAngs_po_dAngSins[i * 3 +2] *arr_dAngs_po_dAngSins[i * 3+2]* marrSinsSyst[2]* marrSinsSyst[2];

       arrSyst[i] =  temp0
               + arr_dAngs_po_dH[i]*arr_dAngs_po_dH[i] * mBeaconDeepthToleran* mBeaconDeepthToleran
               + temp1
               + arr_dAngs_po_dProfile[i]*arr_dAngs_po_dProfile[i] *mvalProfMean *mvalProfMean;
       arrSyst[i] = sqrt(arrSyst[i]);
    }

    // суммарная ошибка
    double arrSum[3] ={0.};
    for (int i = 0; i < 3; ++i)
    {
       arrSum[i] =sqrt(arrSyst[i] * arrSyst[i] +arrSKZ[i] * arrSKZ[i]);
    }

    // выгрузка таблицы антенны
    for (int j =0; j < 3; ++j)
    {
     int ia = ((int)(arrSKZ[j] * 1000.* 10.)) ;
     double a = ((double)ia)/ 10.;
     ui->tableWidget_5->item(2,3 +j)->setText(QString::number(a));

      ia = ((int)(arrSyst[j] * 1000.* 10.)) ;
      a = ((double)ia)/ 10.;
     ui->tableWidget_5->item(3,3 +j)->setText(QString::number(a));

     ia = ((int)(arrSum[ j] * 1000.*10.)) ;
     a = ((double)ia)/ 10.;
    ui->tableWidget_5->item(4,3 + j)->setText(QString::number(a));
    }

}
//----------------------------------------------
// обработчик сигнала о факте завершении потока (этот сигнал поток автоматически отправляет при выходе из функции MyThread::run())
void MainWindow::onTaskFinished()
{
    qDebug() << "MainWindow::onTaskFinished";
    // интерфейс
    ui->pushButtonStartTask->setEnabled(true);
    ui->pushButtonStopTask->setEnabled(false);

    // удаляем завершившийся поток, обнуляем ссылку
    delete myThread;
    myThread = 0;
}

// обработчик сигнала о факте завершения и результате вычислений
void MainWindow::onCalcFinished(int result)
{
    switch(result)
    {
    case MYTHREAD_RESULT_OK:
        ui->labelTaskInfo1->setText("Task finished successfully");
        ui->progressBar->setValue(ui->progressBar->maximum());
        break;
    case MYTHREAD_RESULT_MAX_STEP:
        ui->labelTaskInfo1->setText("Task failed (maximum iterations reached)");
        break;
    case MYTHREAD_RESULT_USER_CANCEL:
        ui->labelTaskInfo1->setText("Task failed (cancelled by user)");
        break;
    case MYTHREAD_RESULT_ERROR:
        ui->labelTaskInfo1->setText("Task failed (cancelled by error data)");
        break;
    default:
        ui->labelTaskInfo1->setText("Task failed (unknown reason)");
        break;
    }
}

void MainWindow::on_tableWidget_6_itemChanged(QTableWidgetItem *item)
{

}

void MainWindow::on_comboBox_5_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:
      ui->groupBox_3 ->hide();
      ui->groupBox->setTitle(QString("РЕШЕНИЕ LBL ЗАДАЧИ"));
        break;

    case 1:
        ui->groupBox_3 ->show();
        ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-2D ЗАДАЧИ"));
        break;

    case 2:
        ui->groupBox_3 ->show();
        ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-3D ЗАДАЧИ"));
        break;
    default: break;
    }
}

void MainWindow::on_comboBox_5_activated(const QString &arg1)
{

}

void MainWindow::on_comboBox_8_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:
        ui->doubleSpinBox_17->setRange(0.,30.);
        ui->doubleSpinBox_17->setValue(20.);
        ui->doubleSpinBox_17->setPrefix(QString("Gape: "));
        ui->doubleSpinBox_17->setSuffix(QString(" гр"));
        ui->doubleSpinBox_17->setDecimals(1);
        ui->doubleSpinBox_17->setSingleStep(1.);
        ui->doubleSpinBox_16->setRange(0.,5.);
        ui->doubleSpinBox_16->setValue(1.);
        ui->doubleSpinBox_16->setPrefix(QString("Step: "));
        ui->doubleSpinBox_16->setSuffix(QString(" гр"));
        ui->doubleSpinBox_16->setSingleStep(0.01);
        ui->doubleSpinBox_16->setDecimals(2);
        break;

    case 1:
        ui->doubleSpinBox_17->setRange(0.,1.1);
        ui->doubleSpinBox_17->setPrefix(QString("Коеф: "));
        ui->doubleSpinBox_17->setSuffix(QString(""));
        ui->doubleSpinBox_17->setDecimals(3);
        ui->doubleSpinBox_17->setSingleStep(0.001);
        ui->doubleSpinBox_17->setValue(0.025);
        ui->doubleSpinBox_16->setRange(0.,1000);
        ui->doubleSpinBox_16->setPrefix(QString("Макс ит.: "));
        ui->doubleSpinBox_16->setSuffix(QString(""));
        ui->doubleSpinBox_16->setDecimals(0);
        ui->doubleSpinBox_17->setSingleStep(50);
        ui->doubleSpinBox_16->setValue(200);

        break;

    default: break;
    }
}

void MainWindow::on_comboBox_8_activated(const QString &arg1)
{

}

void MainWindow::on_pushButton_clicked()
{
    //откуда
    double arr7[50] = {0.};
    read_tableWidget_7(arr7);

    double arr6[50] = {0.};
    read_tableWidget_6(arr6);
    arr6[6] = arr7[0];
    arr6[7] = arr7[1];
    arr6[8] = arr7[2];
    write_tableWidget_6(arr6);
    // !

    double arr5[50] = {0.}, arr4[50] = {0.};
    read_tableWidget_5(arr5);
    read_tableWidget_4(arr4);
    for (int i =0; i <6; ++i)
    {
    arr4[12 + i] = arr5[i];
    }
    write_tableWidget_4(arr4);

}

void MainWindow::on_checkBox_clicked()
{

}

void MainWindow::on_tableWidget_6_cellActivated(int row, int column)
{

}
