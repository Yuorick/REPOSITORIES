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

     this->ui->lineEdit_3->setText("D:\\AKIN\\PROFILES\\02102021\\FILE48.000");
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

  double arr[18] = {0.};
  for (int j =0; j < 6; ++j)
  {
      arr[j] = marrTruePosParams[j];
      arr[12 + j] = marrAprioriPosParams[j];
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
  ui->tableWidget_5->setRowCount(3 );
  ui->tableWidget_5->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_5->verticalHeader()->setVisible(true);

  QStringList lst5;
  lst5<<"X"<<"Y"<<"Z"<<"Bet"<<"Eps"<<"Alf";
  ui->tableWidget_5->setHorizontalHeaderLabels(lst5);
  QStringList lst6;
  lst6<< "Оценка,(м,гр)"<<"Ошибка,(м,мрад)"<<"СКО,(м,мрад)";
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
  ui->tableWidget_7->setRowCount(3 );
  ui->tableWidget_7->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_7->verticalHeader()->setVisible(true);

  QStringList lst7;
  lst7<<"X"<<"Y"<<"Z";
  ui->tableWidget_7->setHorizontalHeaderLabels(lst7);
  QStringList lst8;
  lst8<< "Оценка"<<"Ошибка"<<"СКО";
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

    strFold = this->ui->lineEdit->text();
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;


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
       // СКЗ угловой ошибки СИНС (углов качек)
       double valSigSins = ui->doubleSpinBox_9->value()/ 1000.;

       // СКЗ угловой ошибки СИНСпо скорости углов (углов качек)
        double valSig_d_po_dt_Sins = ui->doubleSpinBox_10->value()/1000.;
         double valMaxSig_Q      =     valSigSins ; //0.000582;
         double valMaxSig_Psi    =      valSigSins ; //0.00145;
         double valMaxSig_Tet    =    valSigSins ; // 0.00145;
         double valMaxSig_dQdt   =      valSig_d_po_dt_Sins ;
         double valMaxSig_dPsidt =      valSig_d_po_dt_Sins ;
         double valMaxSig_dTetdt =      valSig_d_po_dt_Sins ;
         double valK1         = 0.01 ;
         double valSigV  = ui->doubleSpinBox_11->value();
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
        mTimeTempSins = 1./ ui->doubleSpinBox_12->value();
       mSins = QPeaceSins (mEnvironment,valMaxSig_Q, valMaxSig_Psi, valMaxSig_Tet
                            ,valMaxSig_dQdt, valMaxSig_dPsidt,valMaxSig_dTetdt
                            ,valMaxSig_H,valMaxSig_VH,valK1
                            ,valSigV, mT0,mTimeTempSins );
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
                   mTypeOfRLS = LBL;
                   break;


               case 1:
                   mRls_Usbl2D = QRls_Usbl2D((ui->doubleSpinBox_8->value()) * 0.001,0.04
                                 ,(ui->doubleSpinBox_7->value()) * 0.001);
                   mpHidroRLS = &mRls_Usbl2D;
                   mTypeOfRLS = USBL_2D;
                   break;


               case 2:
                   mRls_Usbl3D = QRls_Usbl3D((ui->doubleSpinBox_8->value()) * 0.001,0.04
                     ,(ui->doubleSpinBox_7->value()) * 0.001,(ui->doubleSpinBox_13->value()) * 0.001);
                   mpHidroRLS = &mRls_Usbl3D;
                   mTypeOfRLS = USBL_3D;
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
                            ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2], 0.015);
    /*switch(ui->comboBox_5->currentIndex())
    {
    case 0:
        mLblSolver = QLblSolver(parrMeas, mQuantMeas    ,mtblEstPrfl
                                ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2]);
        mpPosSolver = &mLblSolver;
        break;

    case 1:
        mUsbl_2D_Solver = QUsbl_2D_Solver(parrMeas, mQuantMeas    ,mtblEstPrfl
                                ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2]);
        mpPosSolver = &mUsbl_2D_Solver;
        break;

    case 2:
        mUsbl_3D_Solver = QUsbl_3D_Solver(parrMeas, mQuantMeas    ,mtblEstPrfl
                                ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2]);
        mpPosSolver = &mUsbl_3D_Solver;
        break;

    default:

        break;
    }
*/
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


 if(!mLblSolver.estimateParams(arrX0,arrXRez, valNeviaz0))
    {
        int uu =0;
    }
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
                 ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2],0.0001,arrXRez,arrSBeacon_GSK);

            mpPosSolver = &mSolver2D_Angs;
            bContinueCalculations = true;
            break;

        case 2:
            mSolver3D_Angs = QSolver3D_Angs(parrMeas, mQuantMeas    ,mtblEstPrfl
                       ,marrGPSAprioriPosParams,marrTrueHeadLightPos[2],0.0001,arrXRez,arrSBeacon_GSK);
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

   if(mTypeOfRLS != LBL)
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
    if (!mbtableWidget_4Init)
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
    /*int ncols = ui->tableWidget_4->columnCount();
     for (int i=0; i < ui->tableWidget_4->rowCount(); ++i)
          for (int j=0; j < ncols; ++j)
          {
              ui->tableWidget_4->item(i,j)->setText(QString::number( arr [ i * ncols + j]));

          }*/
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
    double arr[9] = {0.};
    read_tableWidget_6(arr);


    // формирование 3-ей строки
    for (int i =0; i < ui->tableWidget_6->columnCount(); ++i)
    {
        arr[(ui->tableWidget_6->columnCount() *2) + i ] = arr[ i ] + arr[ui->tableWidget_6->columnCount() + i ];

    }

   write_tableWidget_6(arr);

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
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости","D:\\AKIN\\PROFILES" ,"*.000");
   // QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости","..\\POSITION_IMIT_MODEL_RAZD_OBRAB_v0" ,"*.000");
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
   read_tableWidget_5(marrGRLS_TblData);
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

void MainWindow::on_tableWidget_5_itemSelectionChanged()
{
   write_tableWidget_5(marrGRLS_TblData);
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
