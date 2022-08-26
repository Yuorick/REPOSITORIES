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
#include "acdFile.h"
#include "URPolygon.h"




#include "SubWaterBeam.h"
#include "BigMeasure.h"
#include "Gps.h"

#include <QDebug>
#include "Solver3D_Angs.h"

#include <QVector3D>
#include <QDateTime>
#include <QFile> // Подключаем класс QFile
#include <QTextStream> // Подключаем класс QTextStream









#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsReport0;

extern bool BEZ_SHUMOV = false;
#define nullptr 0




MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    myThread(0)
  {

    ui->setupUi(this);

    ui->comboBox_5->setCurrentIndex(2);
    switch(ui->comboBox_5->currentIndex())
    {
    case 0:
        mTypeOfSolverTask = LBL;
     ui->groupBox_tuneUSBL->hide();
        break;

    case 1:
        mTypeOfSolverTask = USBL_2D;
     ui->groupBox_tuneUSBL->show();
        break;

    case 2:
        mTypeOfSolverTask = USBL_3D;
     ui->groupBox_tuneUSBL->show();
        break;

    default:
        break;
    }

   qRegisterMetaType<QVector<double> >("QVector<double>");// для проогресс бара !!
   qRegisterMetaType<QVector<int> >("QVector<int>");

    //******* для прогресс бара **************
       // очистка интерфейса
       ui->pushButtonStartTask->setEnabled(true);
       ui->pushButtonStopTask->setEnabled(false);
       ui->progressBar->setValue(0);
       ui->labelTaskInfo1->setText("");
       ui->labelTaskInfo2->setText("");
       ui->labelTaskInfo3->setText("");
       ui->labelTaskInfo4->setText("");
       ui->labelTaskInfo6->setText("");
       ui->label_2->setText("");


       ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-3D ЗАДАЧИ"));
       mNumOfTask =1;

       mbtableWidget_5Init = false;
       mbtableWidget_7Init = false;
       //**********************!
       QString str0 = QCoreApplication::applicationDirPath() ;
      // str0 = str0 +"//V000002.TXT.conv.000";
      // this->ui->lineEdit_3->setText("D:\\REPOSITORIES\\aircraft-model\PROGRAMS_C++\\AKIN\\SOLVER_POS_v6\\build-SOLVER_POS_v6-Desktop_Qt_5_12_0_MinGW_64_bit-Debug\\debug\\V000002.TXT.conv.000");
       this->ui->lineEdit_3->setText(str0);

     mqstrPrflFIle = this->ui->lineEdit_3->text();
     wchar_t wchPrflFIle [400] = {0};
     mqstrPrflFIle.toWCharArray(wchPrflFIle);
     wchPrflFIle[mqstrPrflFIle.length()] = 0;

     mType_of_000 = VAR0;

     // графики
     ui->lineEdit->setText("D:\\AKIN\\08-06-2022\\Angles");

     ui->lineEdit->setText(str0);
     mqstrOutPutFold = ui->lineEdit->text();


     str0 = QCoreApplication::applicationDirPath() ;
    // str0 = str0 + "\\inputDataExchangeFile.acd";
     this->ui->lineEdit_5->setText("D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\AKIN\\SOLVER_POS_v6\\build-SOLVER_POS_v6-Desktop_Qt_5_12_0_MinGW_64_bit-Debug\\debug\\beacon-3163 _yr_v2.acd");
     this->ui->lineEdit_5->setText(str0);
     mqstrDataFIle = this->ui->lineEdit_5->text();
     wchar_t wchDataFIle [400] = {0};
     mqstrDataFIle.toWCharArray(wchDataFIle);
     wchDataFIle[mqstrDataFIle.length()] = 0;


  // установка параметров позиционирования антенны

  memset(marrAntPosParams, 0, 6 * sizeof(double));

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
        int ia = 10000.* marrAntPosParams[j];
        double temp = ia/10000.;
        ui->tableWidget_4->item(0,j)->setText(QString::number(temp));
    }

    mbtableWidget_4Init = true;
    /// !

    //--------------------------------------------------------------
    // Таблица маяка
    // установка вектора маяка

    memset(marrBeaconPos, 0, 3 * sizeof(double));

    mbtableWidget_6Init = false;

   for (int j = 0; j < ui->tableWidget_6->columnCount(); j++)
   {
       QTableWidgetItem* ptwi0 = nullptr;
       ptwi0 = new QTableWidgetItem(QString::number(0.));

       ui->tableWidget_6->setItem(0,j,ptwi0);
   }


    for (int j =0; j < 3; ++j)
    {
     ui->tableWidget_6->item(0,j)->setText(QString::number(marrBeaconPos[j]));

    }

 mbtableWidget_6Init = true;
    /// !


    // ТАБЛИЦА РЕЗУЛЬТАТОВ ПО АНТЕННЕ
  ui->tableWidget_5->setColumnCount(6);
 // ui->tableWidget_5->setRowCount(2 );
  ui->tableWidget_5->setRowCount(1 );
  ui->tableWidget_5->horizontalHeader()->setVisible(true);

  ui->tableWidget_5->verticalHeader()->setVisible(true);

  QStringList lst5;
  lst5<<"X"<<"Y"<<"Z"<<"Bet"<<"Eps"<<"Alf";
  ui->tableWidget_5->setHorizontalHeaderLabels(lst5);
  QStringList lst6;
  lst6<< "Оценка,(м,гр)"<<"СКО,(м,мрад)";
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
  ui->tableWidget_7->setRowCount(1 );
  ui->tableWidget_7->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_7->verticalHeader()->setVisible(true);

  QStringList lst7;
  lst7<<"X"<<"Y"<<"Z";
  ui->tableWidget_7->setHorizontalHeaderLabels(lst7);
  QStringList lst8;
  lst8<< "Оценка"<<"СКО";
  lst8<< "Оценка";
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

  mbtableWidget_8Init = false;
    for (int i=0; i<  ui->tableWidget_8->columnCount() ; i++)
    {
     QTableWidgetItem* ptwi0 = nullptr;
     ptwi0 = new QTableWidgetItem(QString::number(0.));
     ui->tableWidget_8->setItem(0,i,ptwi0);
    }


 for (int j =0; j < ui->tableWidget_8->columnCount(); ++j)
 {

   ui->tableWidget_8->item(0,j)->setText(QString::number(0.));

 }
 ui->tableWidget_8->setHorizontalHeaderLabels(lst7);
 ui->tableWidget_8->horizontalHeader()->setVisible(true);

 mbtableWidget_8Init = true;  
 mbtableWidget_4Init = true;
 mbtableWidget_6Init = true;
 mbtableWidget_5Init = true;
 mbtableWidget_7Init = true;

 //ui->doubleSpinBox_9->setReadOnly(true);


// ui->frame_2->hide();
// ui->pushButton_3->hide();


//ui->groupBox_5->hide();

 mbOpenDialog_000 = false;
 mbOpenDialog_acd = false;

 ui->comboBox_9->setCurrentIndex(1);

 tuneGroupBoxUSBL(ui->comboBox_9->currentIndex());



   for (int i=0; i<  ui->tableWidget->rowCount() ; i++)
       for (int j=0; j<  ui->tableWidget->columnCount() ; j++)
   {
    QTableWidgetItem* ptwi0 = nullptr;
    ptwi0 = new QTableWidgetItem(QString::number(0.));
    ui->tableWidget->setItem(i,j,ptwi0);
   }


 ui->pushButton_3->hide();
 ui->doubleSpinBox_15->hide();

 // рабочий журнал
 QDateTime datatime  = QDateTime::currentDateTime();
 ui->textEdit->clear();
 ui->textEdit->append("Начало сессии: " + datatime.toString());
 ui->textEdit->append(" Параметры по умолчанию: ");
 QString str = createInputParamsConstructorString();
 ui->textEdit->append(str);

 mNumCkickStart = 0;
 //ui->textEdit->hide();
}
//------------------------------------
void MainWindow::tuneGroupBoxUSBL(int i)
{
    switch(i)
    {
    case 0:
     mALGOR_TYPE = PEREBOR;

     ui->groupBox_9->hide();
     ui->groupBox_7->show();
     ui->groupBox_7->setGeometry(10,60, 291,72);
        break;

    case 1:
     mALGOR_TYPE = NUTON;
     ui->groupBox_7->hide();
     ui->groupBox_9->show();
     ui->groupBox_9->setGeometry(10,60, 291,72);
        break;
    default: break;
    }
}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------


void MainWindow:: inputData()
{
    // 1. файл профиля звука
    mqstrPrflFIle = this->ui->lineEdit_3->text();
    wchar_t wchPrflFIle [400] = {0};
    mqstrPrflFIle.toWCharArray(wchPrflFIle);
    wchPrflFIle[mqstrPrflFIle.length()] = 0;

    // !1

    // 2. директория с результатами (графиками)
    mqstrOutPutFold = this->ui->lineEdit->text();
    wchar_t wchOutPutFold[400];
    mqstrOutPutFold.toWCharArray(wchOutPutFold);
    wchOutPutFold[mqstrOutPutFold.size()] = 0;
    // !2

    // 3.файл с исходными данными
    mqstrDataFIle = this->ui->lineEdit_5->text();
    wchar_t wchDataFIle [400] = {0};
    mqstrDataFIle.toWCharArray(wchDataFIle);
    wchDataFIle[mqstrDataFIle.length()] = 0;

    // !3


   // 4. создание профиля скорости звука
   QSubWaterBeam::createProfileTbl(wchPrflFIle, &mtblEstPrfl,mType_of_000);
  // 4!
  if (ui->comboBox_6->currentIndex() == 1)
  {
    TTable_1D  tblPrfl = mtblEstPrfl.makeMiddleProfile();
    mtblEstPrfl = tblPrfl;
  }



  for (int i = 0; i < 6; ++i)
  {
      if (i <3)
      {
          marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();

      }
      else
      {
          marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;

      }
  }

  read_tableWidget_8(marrGPSPosParams);

  read_tableWidget_6(marrBeaconPos);
  ///

   mMaxQuantIter_LBL =(int)(ui->doubleSpinBox_32->value()+0.0001);
   mMaxQuantIter_USBL =(int)(ui->doubleSpinBox_37->value()+0.0001);



   switch(ui->comboBox_5->currentIndex())
   {
   case 0:
       mTypeOfSolverTask = LBL;
       break;

   case 1:
       mTypeOfSolverTask = USBL_2D;
       break;

   case 2:
       mTypeOfSolverTask = USBL_3D;
       break;

   default:

       break;

   }


   mTetta_min = ui->doubleSpinBox_35->value() * M_PI/180.;

   mTetta_max = ui->doubleSpinBox_36->value()* M_PI/180.;

   mGape = ui->doubleSpinBox_20->value() ;
   mAngStep = ui->doubleSpinBox_19->value() ;

   switch (ui->comboBox_9->currentIndex())
   {
   case 0:
   mALGOR_TYPE = PEREBOR;
       break;

   case 1:
   mALGOR_TYPE = NUTON;
       break;

   default:
       break;
   }

   mLBLcoeff = ui->doubleSpinBox_18->value() ;

   mUSBLcoeff = ui->doubleSpinBox_22->value() ;

   mWorkR =  ui->doubleSpinBox_40->value() ;

   switch (ui->comboBox->currentIndex())
   {
   case 0:
    mTYPE_OF_OUTPUT_FILE = SHP;
       break;

   case 1:
    mTYPE_OF_OUTPUT_FILE = CSV;
       break;

   default:
       break;
   }

}
//-------------------------------------------------------
void MainWindow::on_pushButton_2_clicked()
{
    QString str0 = QCoreApplication::applicationDirPath() ;
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", str0);
    this->ui->lineEdit->setText(strFold);
    QString str("PushButton ВЫБОР ДИРЕКТОРИИ С ГРАФИКОЙ:\n " +strFold) ;
    ui->textEdit->append(str);
    ui->textEdit->append("// !");
}
//--------------------------------------------
void MainWindow::on_pushButton_3_clicked()
{

    inputData();

    double arrSBeacon_XYZ[3] = {0.};
    arrSBeacon_XYZ[0]= marrXRez [3];
    arrSBeacon_XYZ[1]= marrXRez [4];
    arrSBeacon_XYZ[2]= marrBeaconPos [2];
    double arrAntPosParams [6] ={0.};
    memcpy(arrAntPosParams, marrXRez,3 * sizeof(double));
    memcpy(&arrAntPosParams[3], &marrXRez[5],3 * sizeof(double));
    wchar_t wchOutPutFold0[400] = {0},wchOutPutFold[400] = {0};
    mqstrOutPutFold.toWCharArray(wchOutPutFold0);
    wchOutPutFold0[mqstrOutPutFold.size()] = 0;
    wcscpy(wchOutPutFold, wchOutPutFold0);

    
    QString qstr = this->ui->lineEdit_5->text();

   // readHeaderDataFile(qstr,mQuantMeas, marrAntPosParams
                              //          , marrGPSPosParams, marrBeaconPos, mTObrabotki);
   /* for (int i = 0; i < 6; ++i)
    {
        if (i <3)
        {
            marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();

        }
        else
        {
            marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;

        }
    }
*/
    read_tableWidget_8(marrGPSPosParams);

   // read_tableWidget_6(marrBeaconPos);
   // marrBeaconPos[0] =  ui->tableWidget_7->item(0,0)->text().toDouble();
   // marrBeaconPos[1] =  ui->tableWidget_7->item(0,1)->text().toDouble();
  // marrBeaconPos[2] =  ui->tableWidget_7->item(0,2)->text().toDouble();

    QVector<QBigMeasure> vctBigMesure(mQuantMeas0);
    readBigMeasuresOnly_acd(qstr,arrAntPosParams
                                 , marrGPSPosParams, arrSBeacon_XYZ,mTObrabotki
                                 , &vctBigMesure);

    int QuantMeas = vctBigMesure.size();
    QBigMeasure *parrMeas = (QBigMeasure *)malloc(QuantMeas * sizeof(QBigMeasure) );
    for(int i =0;i <QuantMeas; ++i)
    {
      parrMeas[i] =  vctBigMesure.at(i);
    }
    TURPointXY *pnts = new TURPointXY[2 *QuantMeas];
    TURPointXY *pntsWave_ = new TURPointXY[QuantMeas];
    TURPointXY *pnts_ = new TURPointXY[QuantMeas];
    for (int i =0; i< QuantMeas; ++i)
    {
        pnts[2 *i].X = parrMeas[i].marrSVessWaveZv[0] - arrSBeacon_XYZ[0];
        pnts[2 *i ].Y = parrMeas[i].marrSVessWaveZv[1]-arrSBeacon_XYZ[1];

        pnts[2 *i +1 ].X = parrMeas[i].marrSVessZv[0] -arrSBeacon_XYZ[0];
        pnts[2 *i +1].Y = parrMeas[i].marrSVessZv[1]-arrSBeacon_XYZ[1];

        pntsWave_[i].X = parrMeas[i].marrSVessWaveZv[0]-arrSBeacon_XYZ[0];
        pntsWave_[i].Y = parrMeas[i].marrSVessWaveZv[1]-arrSBeacon_XYZ[1];

        pnts_[i].X = parrMeas[i].marrSVessZv[0]-arrSBeacon_XYZ[0];
        pnts_[i].Y = parrMeas[i].marrSVessZv[1]-arrSBeacon_XYZ[1];
       int uu=0;
    }
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TURPointXY pnt00(marrBeaconPos[0], marrBeaconPos[1]);
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.,pnt00) ;

    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\PntsTraj.shp");
    pnts[0].WriteSetSHPFiles(wchAxesFileName0,pnts,2 *mQuantMeas0);

    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\PntsTrajWave_.shp");
    pntsWave_[0].WriteSetSHPFiles(wchAxesFileName0,pntsWave_,mQuantMeas0);

    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\PntsTraj_.shp");
    pnts_[0].WriteSetSHPFiles(wchAxesFileName0,pnts_,mQuantMeas0);




    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\arrAprioriBeaconPos.shp");
    TURPointXY arrAprioriBeaconPos(marrBeaconPos[0], marrBeaconPos[1]);
    arrAprioriBeaconPos.WriteSetSHPFiles(wchAxesFileName0,&arrAprioriBeaconPos,1);

    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\pntFirst.shp");
    pnts[0].WriteSetSHPFiles(wchAxesFileName0,&pnts[0],1);

    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\pntEnd.shp");

    pnts[mQuantMeas0-1].WriteSetSHPFiles(wchAxesFileName0,&pnts[mQuantMeas0-1],1);

    TURPolyLine pln(pnts, 2 *mQuantMeas0) ;
    wcscpy(  wchAxesFileName0,  wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\plnTraj.shp");

    pln.WriteSetSHPFiles(wchAxesFileName0,&pln,1);
    delete []pnts;
    delete []pntsWave_;
    delete []pnts_;

    QLblSolver LblSolver(parrMeas, QuantMeas   ,mtblEstPrfl
                              ,arrSBeacon_XYZ[2], 0.015);

      double valNeviaz0 = -1.;



      QPosSolver *pPosSolver = &LblSolver;
      double sigmaTreshold = 2.5;
      pPosSolver->repaireZamerArray_(marrXRez, sigmaTreshold,0, NULL);



//***************************************
      QSolver3D_Angs Solver3D_Angs;
      QSolver2D_Angs Solver2D_Angs;
      QSolver2D_Angs *pSolver2D_Angs;
  switch(mTypeOfSolverTask)
  {
  case USBL_2D:
    Solver2D_Angs =   QSolver2D_Angs (parrMeas, QuantMeas,mtblEstPrfl
           ,arrSBeacon_XYZ[2],0.009,marrXRez,arrSBeacon_XYZ);
    pSolver2D_Angs = &Solver2D_Angs;
      break;
  case USBL_3D:
      Solver3D_Angs =   QSolver3D_Angs (parrMeas, QuantMeas,mtblEstPrfl
                                        ,arrSBeacon_XYZ[2],0.009,marrXRez,arrSBeacon_XYZ);
      pSolver2D_Angs = &Solver3D_Angs;
      break;

  default:
      return;
  }





  delete []parrMeas;

  //Графики
  double arrSKZ2[3] = {0.};
  double *parrBuffNeviaz = new double [QuantMeas * (pSolver2D_Angs->mDimY +1)];
  pSolver2D_Angs ->calcNeviazArray_(&marrXRez[5], parrBuffNeviaz,arrSKZ2);

  QuantMeas = pSolver2D_Angs->mVectBigMeasures.size();
  int iLenName = 30;
  wchar_t pwcharrColNames [4 * iLenName];

  wcscpy(pwcharrColNames, L"n");
  wcscpy(&pwcharrColNames[iLenName   ],  L"ant_bet");
  wcscpy(&pwcharrColNames[iLenName *2],  L"ant_eps");
  wcscpy(&pwcharrColNames[iLenName *3],  L"deltaq_from_deltaQ");

  double arrScale[100];
  for (int i = 0; i < 100;++i)
  {
   arrScale[i] = 180./M_PI;
  }





for (int i =1; i < (1 +pSolver2D_Angs ->mDimY); i++)
{
  TYrWriteShapeFile::WriteOneReport(wchOutPutFold// путь к папке
                                  ,parrBuffNeviaz // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,(1 + pSolver2D_Angs ->mDimY)  // - к-во переменных о корорых накоплена информация в буфере
                                  ,QuantMeas //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,0  //  номер переменной по оси X
                                  ,i  //  номер переменной по оси Y
                                  ,1.  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                   ) ;
}

// невязка по КУ в зависимости от КУ
double *parrBuffNeviaz1 = new double [QuantMeas * 2];

wchar_t wcharrColNames [400]= {0};

wcscpy(wcharrColNames, L"Q");
wcscpy(&wcharrColNames[iLenName   ],  L"nev_q");

QBigMeasure meas= pSolver2D_Angs->mVectBigMeasures.at(0);
double q0 = meas.marrMuZv[0];

parrBuffNeviaz1[0]  =q0 / M_PI * 180.;
parrBuffNeviaz1[ 1] = parrBuffNeviaz[1]/ M_PI * 180.;

for (int i =1; i < QuantMeas; ++i)

{
    QBigMeasure meas= pSolver2D_Angs->mVectBigMeasures.at(i);
    double q = fmod(meas.marrMuZv[0], 2. * M_PI);


    parrBuffNeviaz1[2 * i]  = q/ M_PI * 180.;
    parrBuffNeviaz1[2 * i + 1] = parrBuffNeviaz[(1 + pSolver2D_Angs ->mDimY) * i +1]/ M_PI * 180.;

}


TYrWriteShapeFile::WriteOneReport(wchOutPutFold// путь к папке
                                ,parrBuffNeviaz1 // массив с информацией - матрица nBuffRows x nBuffCols
                                ,2  // - к-во переменных о корорых накоплена информация в буфере
                                ,QuantMeas //  - к-во точек
                                ,wcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                ,iLenName // максимальная длина имени переменной
                                ,0 //  номер переменной по оси X
                                ,1  //  номер переменной по оси Y
                                ,1.  //  масштаб по оси X
                                ,50.// масштаб по оси Y
                                 ) ;

delete []parrBuffNeviaz;
delete []parrBuffNeviaz1;


wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntEstBeaconPos.shp");
TURPointXY pntEstBeaconPos(marrXRez[3], marrXRez[4]);
pntEstBeaconPos.WriteSetSHPFiles(wchAxesFileName0,&pntEstBeaconPos,1);



wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\Axes00.shp");
TYrWriteShapeFile::CreateShpAxes(wchAxesFileName0,-5000.,5000.
     ,-2000.,2000.);
  //!

int NumMeas = (int)(ui->doubleSpinBox_15->value()+0.1);
createPanoramePictXY(wchOutPutFold, marrBeaconPos, marrAntPosParams, pSolver2D_Angs ->mVectBigMeasures[NumMeas]);

// облако точек
readBigMeasuresOnly_acd(qstr,marrAntPosParams
                             , marrGPSPosParams, marrBeaconPos,mTObrabotki
                             , &vctBigMesure);
mQuantMeas0 = vctBigMesure.size();
TURPointXY *pntsXY = new TURPointXY[mQuantMeas0];
TURPointXY *pntsXZ= new TURPointXY[mQuantMeas0];
TURPointXY *pntsZY = new TURPointXY[mQuantMeas0];
int qPnts = 0;


double val_q = 0., val_e =0., val_t = 0.;
double arrFieldMean[3] = {0.};

double *parrSBeaconField = new double [3 *mQuantMeas0];
for (int i =0; i < QuantMeas; ++i)
{
    QBigMeasure meas = vctBigMesure.at(i);
    if (mTypeOfSolverTask == USBL_2D)
    {
        if(!transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, arrSBeacon_XYZ, meas.marrSVessZv, meas.marrMuZv
          ,arrAntPosParams,  &val_q, &val_e,  &val_t))
        {
            continue;
        }
        meas.mezv =val_e;
    }
    double valTZv = 0.;

    if(!transfMeasure_to_GSK(mtblEstPrfl,meas,arrAntPosParams,NULL
                              ,  &parrSBeaconField[3 * qPnts], &valTZv))
    {
        continue;
    }

    if (Norm3(&parrSBeaconField[3 * qPnts]) < 0.001)
    {
        int yyy=0;
    }
   MtrxMinusMatrx(&parrSBeaconField[3 * qPnts], arrSBeacon_XYZ,1, 3,&parrSBeaconField[3 * qPnts]) ;
    double *p =&parrSBeaconField[3 * qPnts];

    pntsXY[qPnts].X =  p[0];
    pntsXY[qPnts].Y =  p[1];

    pntsXZ[qPnts].X =  p[0];
    pntsXZ[qPnts].Y =  p[2];

    pntsZY[qPnts].X =  p[1];
    pntsZY[qPnts].Y =  p[2];

    qPnts++;

}

 for (int i =0 ;i < qPnts; ++i)
{
    MtrxSumMatrx(arrFieldMean, &parrSBeaconField[3 * i],1, 3, arrFieldMean) ;

}
MatrxMultScalar(arrFieldMean,1, 3, 1./ ((double)qPnts),arrFieldMean);
double arrDisp[3] = {0.};
for (int i =0 ;i < qPnts; ++i)
{
    double arrt[3]= {0.};
    MtrxMinusMatrx(&parrSBeaconField[3 * i], arrFieldMean,1, 3,arrt) ;
    arrDisp[0] += arrt[0] * arrt[0];
    arrDisp[1] += arrt[1] * arrt[1];
    arrDisp[2] += arrt[2] * arrt[2];

}
MatrxMultScalar(arrDisp, 1, 3, 1./ ((double)qPnts -1.),arrDisp);
arrDisp[2] +=0.01;

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntFieldXY.shp");
pntsXY[0].WriteSetSHPFiles(wchAxesFileName0,pntsXY,qPnts);

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntFieldXZ.shp");
pntsXY[0].WriteSetSHPFiles(wchAxesFileName0,pntsXZ,qPnts);

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntFieldZY.shp");
pntsXY[0].WriteSetSHPFiles(wchAxesFileName0,pntsZY,qPnts);


double arrElK_XY[4 ] ={0.};
arrElK_XY[0] = arrDisp[0];
arrElK_XY[3] = arrDisp[1];
TURPointXY pointCentre(0.,0.);
TURPolygon   ellXY =TURPolygon::fncCreateEllipse_( pointCentre
                                                   ,2. * sqrt(arrDisp[0]),2. * sqrt(arrDisp[1]), 501);
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\ElXY.shp");
ellXY.WriteSetSHPFiles(wchAxesFileName0,&ellXY,1);

double arrElK_XZ[4 ] ={0.};
arrElK_XZ[0] = arrDisp[0];
arrElK_XZ[3] = arrDisp[2];

TURPolygon   ellXZ =TURPolygon::fncCreateEllipse_( pointCentre
                                                  ,2. * sqrt(arrDisp[0]),2. * sqrt(arrDisp[2]), 501);
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\ElXZ.shp");
ellXY.WriteSetSHPFiles(wchAxesFileName0,&ellXZ,1);

double arrElK_ZY[4 ] ={0.};
arrElK_ZY[0] = arrDisp[1];
arrElK_ZY[3] = arrDisp[2];
TURPolygon   ellZY =TURPolygon::fncCreateEllipse_( pointCentre
                                                   ,2. * sqrt(arrDisp[1]),2. * sqrt(arrDisp[2]), 501);
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\ElZY.shp");
ellXY.WriteSetSHPFiles(wchAxesFileName0,&ellZY,1);

TURPointXY pntMeanXY(arrFieldMean[0], arrFieldMean[1]);
TURPointXY pntMeanXZ(arrFieldMean[0], arrFieldMean[2]);
TURPointXY pntMeanZY(arrFieldMean[1], arrFieldMean[2]);

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntMeanXY.shp");
pntMeanXY.WriteSetSHPFiles(wchAxesFileName0,&pntMeanXY,1);

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntMeanXZ.shp");
pntMeanXZ.WriteSetSHPFiles(wchAxesFileName0,&pntMeanXZ,1);

wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\pntMeanZY.shp");
pntMeanZY.WriteSetSHPFiles(wchAxesFileName0,&pntMeanZY,1);


delete []parrSBeaconField;

delete []pntsXY;
delete []pntsXZ;
delete []pntsZY;

}
//-------------------------------------------------------
void MainWindow::createPanoramePictXY(wchar_t *wchOutPutFold, double *arrBeaconPos, double *arrAntPosParams, const QBigMeasure &Meas)
{
wchar_t wchAxesFileName0[300] ={0};
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\AxesArr.shp");
TURPointXY pnt00(arrBeaconPos[0], arrBeaconPos[1]);
TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
,-10000., 10000.,30.,pnt00) ;

// вектор корабля на приеме
double length0= 100.;
double valQ = Meas.marrMuZv[0];
TURPointXY pntVess0(Meas.marrSVessZv[0],Meas.marrSVessZv[1]);
TURPointXY pntVess1(Meas.marrSVessZv[0] + length0 * sin(valQ),Meas.marrSVessZv[1]+ length0 * cos(valQ));
TURPolyLine pllVess( pntVess0, pntVess1) ;
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\plnVess.shp");
pllVess.WriteSetSHPFiles(wchAxesFileName0,&pllVess,1);
// !

// антенная сист кординат
TURPolyLine pln_aspk0 = TURPolyLine::fncCreateAxes(TURPointXY(0, 0.)
                                       ,TURPointXY(1000, 0.),TURPointXY(0, 0.),TURPointXY(0, 1000.)
                                      ,100.);
// !

//
TURPolyLine pln_psk = pln_aspk0.LinTransform(-valQ, TURPointXY(0.,0.)
                                               , pntVess0,1. );
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\plnPSK.shp");
pln_psk.WriteSetSHPFiles(wchAxesFileName0,&pln_psk,1);


double angRot = valQ + arrAntPosParams[3];//

TURPolyLine pln_aspk1 = pln_aspk0.LinTransform(-angRot, TURPointXY(0.,0.)
                                               , pntVess0,1. );
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\plnASPK.shp");
pln_aspk1.WriteSetSHPFiles(wchAxesFileName0,&pln_aspk1,1);

// !

double valBeaconDirZv = valQ + arrAntPosParams[3] +Meas.mqzv;
valBeaconDirZv =fmod(valBeaconDirZv, 2. * M_PI);
double length1= 1000.;
TURPointXY pntBeaconDir(Meas.marrSVessZv[0] + length1 * sin(valBeaconDirZv),Meas.marrSVessZv[1]+ length1 * cos(valBeaconDirZv));
TURPolyLine plnBeaconDirZv( pntVess0, pntBeaconDir) ;
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\plnBeaconDirZv.shp");
plnBeaconDirZv.WriteSetSHPFiles(wchAxesFileName0,&plnBeaconDirZv,1);
//

TURPolyLine plnBeaconDir( pntVess0, pnt00) ;
wcscpy(  wchAxesFileName0,  wchOutPutFold);
wcscat(wchAxesFileName0, L"\\plnBeaconDir.shp");
plnBeaconDir.WriteSetSHPFiles(wchAxesFileName0,&plnBeaconDir,1);

double rrealQ = atan2(pnt00.X-pntVess0.X, pnt00.Y-pntVess0.Y);
if (rrealQ < 0.)
{
    rrealQ += 2. * M_PI;
}

double nev = valBeaconDirZv - rrealQ;
double nev1 = fmod(nev, 2. * M_PI);

double nev2 = nev1 - 2. * M_PI;
if (fabs (nev2) < fabs(nev1))
{
   nev =nev2;
}
else
{
  nev =nev1;
}

int yt =0;

}
//---------------------------------------------------------
void MainWindow::read_tableWidget_8(double *arr2)
{
     for (int i=0; i < ui->tableWidget_8->rowCount(); ++i)
    {
          arr2 [i] =  ui->tableWidget_8->item(i,0)->text().toDouble();
    }
}
//-------------------------------------------------
void MainWindow::read_tableWidget_6(double *arrAprioriBeaconPos)
{
     for (int i=0; i < ui->tableWidget_6->columnCount(); ++i)
    {

          arrAprioriBeaconPos[i] =  ui->tableWidget_6->item(0,i)->text().toDouble();
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
             QString strNumber = QString::number(arr [ i * ncols + j], 'f', 2);
             ui->tableWidget_7->item(i,j)->setText(strNumber);
         }
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
//------------------------------------------------------

//---------------------------------------------------------------
void MainWindow::readHeaderDataFile(QString DataFileName,int &iN, double *arrP_Ant
                                    , double *arrP_Gps, double *arrS_Beacon, double &valTBeacon)
{
QFile TargFile;
TargFile.setFileName(DataFileName);
TargFile.open(QIODevice::ReadOnly | QIODevice::Text);
QString line0 = TargFile.readLine();
line0 = TargFile.readLine();

iN = line0.toInt();


line0 = TargFile.readLine();
line0 = TargFile.readLine();

QStringList w =line0.split(";");
for(int i =0; i < 3;++i )
{
    arrP_Ant[i] = w.at(i).toDouble();
}
for(int i =3; i < 6;++i )
{
    arrP_Ant[i] = w.at(i).toDouble() * M_PI/180.;
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

TargFile.close();
}
/*
//---------------------------------------
//Чтение файла с исх данными и формиорование массивов отсеков по типам воздействия
void MainWindow::readDataFile_(QString TargFileName,const int QuantMeas
                                              ,QBigMeasure *parrMeas)
{
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
   for (int i =0;i < QuantMeas; ++i)
   {
       line0 = TargFile.readLine();
       QStringList w =line0.split(";");
       parrMeas [i].mTzaprZv = w.at(0).toDouble();

       parrMeas [i].marrSVessWaveZv[0] = w.at(1).toDouble();
       parrMeas [i].marrSVessWaveZv[1] = w.at(2).toDouble();
       parrMeas [i].marrSVessWaveZv[2] = w.at(3).toDouble();

       parrMeas [i].marrMuWaveZv[0] = w.at(4).toDouble() * M_PI/ 180.;
       parrMeas [i].marrMuWaveZv[1] = w.at(5).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuWaveZv[2] = w.at(6).toDouble()* M_PI/ 180.;

       double arrGps_KGSK[3] = {0.};
       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( marrGPSPosParams, parrMeas [i].marrMuWaveZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(parrMeas [i].marrSVessWaveZv, arrGps_KGSK,1, 3, parrMeas [i].marrSVessWaveZv);


       parrMeas [i].mTotvZv = w.at(7).toDouble();

       parrMeas [i].marrSVessZv[0] = w.at(8).toDouble();
       parrMeas [i].marrSVessZv[1] = w.at(9).toDouble();
       parrMeas [i].marrSVessZv[2] = marrBeaconPos[2];//w.at(10).toDouble();

       parrMeas [i].marrMuZv[0] = w.at(11).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuZv[1] = w.at(12).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuZv[2] = w.at(13).toDouble()* M_PI/ 180.;

       QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( marrGPSPosParams, parrMeas [i].marrMuZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(parrMeas [i].marrSVessZv, arrGps_KGSK,1, 3, parrMeas [i].marrSVessZv);

       parrMeas [i].mSig_t = (ui->doubleSpinBox_8->value()) * 0.001;
       parrMeas [i].mTobr = mTObrabotki;
       ++iChetchik;
   }

    TargFile.close();
}
//---------------------------------------
//Чтение файла с исх данными и формиорование массивов отсеков по типам воздействия
void MainWindow::readDataFile(QString TargFileName,const int QuantMeas
                                              ,QBigMeasure *parrMeas)
{
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
   for (int i =0;i < QuantMeas; ++i)
   {
       line0 = TargFile.readLine();
       QStringList w =line0.split(";");
       parrMeas [i].mTzaprZv = w.at(0).toDouble();

       parrMeas [i].marrSVessWaveZv[0] = w.at(1).toDouble();
       parrMeas [i].marrSVessWaveZv[1] = w.at(2).toDouble();
       parrMeas [i].marrSVessWaveZv[2] = w.at(3).toDouble();

       parrMeas [i].marrMuWaveZv[0] = w.at(4).toDouble() * M_PI/ 180.;
       parrMeas [i].marrMuWaveZv[1] = w.at(5).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuWaveZv[2] = w.at(6).toDouble()* M_PI/ 180.;

       double arrGps_KGSK[3] = {0.};
     QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( marrGPSPosParams, parrMeas [i].marrMuWaveZv
                                             , NULL,arrGps_KGSK,3 );
      MtrxMinusMatrx(parrMeas [i].marrSVessWaveZv, arrGps_KGSK,1, 3, parrMeas [i].marrSVessWaveZv);


       parrMeas [i].mTotvZv = w.at(7).toDouble();

       parrMeas [i].marrSVessZv[0] = w.at(8).toDouble();
       parrMeas [i].marrSVessZv[1] = w.at(9).toDouble();
       parrMeas [i].marrSVessZv[2] = w.at(10).toDouble();

       parrMeas [i].marrMuZv[0] = w.at(11).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuZv[1] = w.at(12).toDouble()* M_PI/ 180.;
       parrMeas [i].marrMuZv[2] = w.at(13).toDouble()* M_PI/ 180.;

       parrMeas [i].mqzv = w.at(14).toDouble()* M_PI/ 180.;
       parrMeas [i].mezv = w.at(15).toDouble()* M_PI/ 180.;

      QPeaceVess::RecalcVect_PSK_CT_INTO_KGSK( marrGPSPosParams, parrMeas [i].marrMuZv
                                                , NULL,arrGps_KGSK,3 );
       MtrxMinusMatrx(parrMeas [i].marrSVessZv, arrGps_KGSK,1, 3, parrMeas [i].marrSVessZv);

       parrMeas [i].mSig_t = (ui->doubleSpinBox_8->value()) * 0.001;
       parrMeas [i].mTobr = mTObrabotki;
       ++iChetchik;
   }

    TargFile.close();
}
*/
//-------------------------------------------------
void MainWindow::fillFields()
{
    fill_tableWidget_4(marrAntPosParams);
    fill_tableWidget_6(marrBeaconPos);
    fill_tableWidget_8(marrGPSPosParams);
}
//----------------------------------------------------
void MainWindow::fill_tableWidget_4(double *arr)
{

    for (int j=0; j < 3; ++j)
    {
        ui->tableWidget_4->item(0,j)->setText(QString::number( arr [j]));

    }
    for (int j=3; j < 6; ++j)
    {
        ui->tableWidget_4->item(0,j)->setText(QString::number( arr [j]* 180./M_PI));

    }
}

//--------------------------------------------
void MainWindow::fill_tableWidget_6(double *arr)
{
    for (int j=0; j < 3; ++j)
    {
        QString strNumber = QString::number(arr [j], 'f', 2);
       // ui->tableWidget_6->item(0,j)->setText(QString::number( arr [j]));
                ui->tableWidget_6->item(0,j)->setText(strNumber);

    }
}

//--------------------------------------------
void MainWindow::fill_tableWidget_8(double *arr)
{
    for (int j=0; j < 3; ++j)
    {
        ui->tableWidget_8->item(0,j)->setText(QString::number( arr [j]));

    }
}



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
   ++mNumCkickStart;




  if(!mbOpenDialog_acd)
  {
      QMessageBox msgBox;
      msgBox.setWindowTitle("Внимание!");
      msgBox.setText("Укажите путь к <*.acd> файлу с данными калибровки!!!");
      msgBox.exec();
      ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ. Не указан путь к <*.acd> файлу \n"));
      ui->textEdit->append(QString("// ! \n"));
      return;
  }
  if (! mbOpenDialog_000)

    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Внимание!");
        msgBox.setText("Укажите путь к <*.000> файлу с верт. профилем скорости звука!!!");
        ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ. Не указан путь к <*.000> файлу \n"));
        ui->textEdit->append(QString("// ! \n"));
        msgBox.exec();
        return;
    }
/*  ui->doubleSpinBox_46->setVisible(false);
  ui->doubleSpinBox_47->setVisible(false);
  ui->labelTaskInfo1->setText("");
  ui->labelTaskInfo2->setText("");
  ui->labelTaskInfo3->setText("");
  ui->labelTaskInfo4->setText("");

    setTableZero(*(ui->tableWidget_5));
    setTableZero(*(ui->tableWidget_7));*/
  cleanInteface();

    // создаем поток и привязываем сигналы к слотам
    myThread = new MyThread;
    QObject::connect(myThread, SIGNAL(finished()), this, SLOT(onTaskFinished()));

    QObject::connect(myThread, SIGNAL(progress(int , QVector<double>, QVector<double>
                           , QVector<double> , int  ,bool , QVector<double>))
            , this, SLOT(onTaskProgress(int , QVector<double>, QVector<double>
                                        , QVector<double> , int  ,bool , QVector<double>)));

    QObject::connect(myThread, SIGNAL(kleaningEnded(int , QVector<int>, QVector<double>, QVector<double>))
            , this, SLOT(onKleaningEnded(int , QVector<int>, QVector<double>, QVector<double>)));

    QObject::connect(myThread, SIGNAL(calcfinished(int)), this, SLOT(onCalcFinished(int)));
    qDebug() << "MainWindow::on_pushButtonStartTask_clicked";

    // заполнение исходеых данных с формы
    inputData();

    QString qstrInp = createInputParamsString(mNumCkickStart);
    ui->textEdit->append(qstrInp);
    // !

    switch(mTypeOfSolverTask)
    {
    case LBL:
    ui->groupBox->setTitle(QString("РЕШЕНИЕ LBL ЗАДАЧИ"));
        break;
    case USBL_2D:
    ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-2D ЗАДАЧИ"));
        break;
    case USBL_3D:
    ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-3D ЗАДАЧИ"));
        break;
    default: break;
    }



    //чтение файла с исходными данными
   QString qstr = this->ui->lineEdit_5->text();

    QVector <QBigMeasure> vctMeas;


    readBigMeasuresOnly_acd(qstr,marrAntPosParams
                            , marrGPSPosParams, marrBeaconPos, mTObrabotki,&vctMeas);



    ui->doubleSpinBox_41->setValue(vctMeas.size());
    for (int i =0; i < vctMeas.size(); ++i)
    {
        QBigMeasure meas = vctMeas.at(i);
        meas.mSig_t =  0.001;
        meas.mSig_e = 0.02 * 0.001;
        meas.mSig_q = 0.02 * 0.001;
       // meas.marrMuZv[0] *= -1.;
       // meas.marrMuWaveZv[0] *= -1.;
        vctMeas.replace(i, meas);
    }

    // очистка интерфейса
    ui->pushButtonStartTask->setEnabled(false);
    ui->pushButtonStopTask->setEnabled(true);
    ui->progressBar->setValue(0);
    ui->progressBar->setRange(0, mMaxQuantIter_LBL);
    ui->labelTaskInfo1->setText("Этап 1: вычисление офсетов...");
    ui->labelTaskInfo2->setText("");
    ui->labelTaskInfo3->setText("");




    // задаём параметры выполнения задачи
    // здесь можно передать массив измерений, профиль скорости звука и пр. параметры

    marrXRez[0]= marrAntPosParams[0];
    marrXRez[1]= marrAntPosParams[1];
    marrXRez[2]= marrAntPosParams[2];

    marrXRez[3]= marrBeaconPos[0];
    marrXRez[4]= marrBeaconPos[1];


    marrXRez[5]= marrAntPosParams[3];
    marrXRez[6]= marrAntPosParams[4];
    marrXRez[7]= marrAntPosParams[5];


    QVector <double>vctX0(8);


    for (int i = 0; i < 8; ++i)
    {
      vctX0.replace(i,marrXRez[i]);

    }


// чистка замеров по геометрии
    double valTZv =-1.;
    int  nun_dist= 0 //удалены по превышению расстояния
            ,  nun_UM =0  // по ограничению УМ
            ,  nun_notCorrect =0;// по некорректности

    const int NUm0 = vctMeas.size();
    for (int i = (NUm0 - 1); i >=0; --i )
    {
        if (vctMeas.size() < 10)
        {
            QMessageBox msgBox;
            msgBox.setWindowTitle("Внимание!");
            msgBox.setText("Хороших замеров менее 10. Проверьте настройки !");
            // очистка интерфейса
            ui->pushButtonStartTask->setEnabled(true);
            ui->pushButtonStopTask->setEnabled(false);
            ui->progressBar->setValue(0);
            ui->labelTaskInfo1->setText("");
            ui->labelTaskInfo2->setText("");
            ui->labelTaskInfo3->setText("");
            ui->labelTaskInfo4->setText("");
            ui->labelTaskInfo6->setText("");
            ui->label_2->setText("");

            ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ ПО ОТСЕВУ ПО ГЕОМЕТРИИ: ОСТАЛОСЬ МЕНЕЕ 10 ИЗМЕРЕНИЙ. УДАЛЕНЫ: \n"));
            ui->textEdit->append(QString(" по ограничению  Д  %0 замеров\n").arg(nun_dist));
            ui->textEdit->append(QString(" по ограничению УМ  %0 замеров\n").arg(nun_UM));
            ui->textEdit->append(QString(" по по несовместимости  %0 замеров\n").arg(nun_notCorrect));
            ui->textEdit->append(QString("// !\n"));
            msgBox.exec();

            return;
        }


        QBigMeasure meas = vctMeas.at(i);
        double arr[2]={0.};
        MtrxMinusMatrx(meas.marrSVessWaveZv, marrBeaconPos,1, 2, arr);
        if(NormVect2(arr) > mWorkR)
        {
            vctMeas.remove(i);
            nun_dist++;
            continue;
        }
        if ( USBL_3D == mTypeOfSolverTask)
        {

            if((fabs(meas.mezv) < mTetta_min)||( fabs(meas.mezv) > mTetta_max))
            {
                vctMeas.remove(i);
                nun_UM++;
                continue;
            }


        }
        else
        {
            //
            double val_q = 0., val_e= 0.,  val_t = 0.;
            if (!transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, marrBeaconPos, meas.marrSVessWaveZv, meas.marrMuWaveZv
             ,marrAntPosParams,  &val_q, &val_e,  &val_t))

            {
               nun_notCorrect++;
               vctMeas.remove(i);
               continue;
            }

            if((fabs(val_e) < mTetta_min)|| (fabs(val_e) > mTetta_max))
            {

                vctMeas.remove(i);
                continue;
            }


            if (!transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, marrBeaconPos, meas.marrSVessZv, meas.marrMuZv
              ,marrAntPosParams,  &val_q, &val_e,  &val_t))
            {
                nun_notCorrect++;
                vctMeas.remove(i);
                continue;
            }
            if((fabs(val_e) < mTetta_min)|| (fabs(val_e) > mTetta_max))
            {
                nun_UM++;
                vctMeas.remove(i);
                continue;
            }
        }
    }

    const int numZam = vctMeas.size();
    ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ОТСЕВА ПО ГЕОМЕТРИИ:осталось %0 замеров,удалены %0 замеров \n").arg(numZam).arg(NUm0 - numZam));
    ui->doubleSpinBox_45->setValue((double)numZam );


    //чистка по максимальной невязке

    //по задержке времени

      // по времени
    QBigMeasure *parrMeas = (QBigMeasure *)malloc((vctMeas.size()) * sizeof(QBigMeasure) );
    for (int i = 0; i < vctMeas.size(); ++i )
    {
       parrMeas[i] =  vctMeas.at(i);
    }
   QLblSolver LblSolver(parrMeas, vctMeas.size()    ,mtblEstPrfl
                           ,marrBeaconPos[2], 0.015);

   double arrX0[8] ={0.},arrXRez [8] = {0.}, valNeviaz0 = -1.;
   for (int i = 0; i < 8; ++i)
   {
      arrX0[i] =  vctX0.at(i);
   }
   QPosSolver *pPosSolver = &LblSolver;


   QVector<int> vctGoodZamersNum;

     double valTreshold = ui->doubleSpinBox_44->value()* 1000.;
     pPosSolver->fncRoughCleaningMeasures(arrX0,valTreshold, 0, &vctGoodZamersNum);
     const int numZam1 = pPosSolver->mVectBigMeasures.size();
     if (numZam1 < 10)
     {
         QMessageBox msgBox;
         msgBox.setWindowTitle("Внимание!");
         msgBox.setText("Отсев по невязке времени. Хороших замеров менее 10.\n Проверьте настройки !");
         // очистка интерфейса
         ui->pushButtonStartTask->setEnabled(true);
         ui->pushButtonStopTask->setEnabled(false);
         ui->progressBar->setValue(0);
         ui->labelTaskInfo1->setText("");
         ui->labelTaskInfo2->setText("");
         ui->labelTaskInfo3->setText("");
         ui->labelTaskInfo4->setText("");
         ui->labelTaskInfo6->setText("");
         ui->label_2->setText("");
         ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ ПО ПРЕДВ. ОТСЕВУ ПО НЕВ. ВРЕМЕНИ: ОСТАЛОСЬ МЕНЕЕ 10 ИЗМЕРЕНИЙ. \n"));
         ui->textEdit->append(QString("// !\n"));
         msgBox.exec();
         return;
     }

    int ncur = numZam - numZam1;
     ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ПРЕДВ. ОТСЕВА ПО НЕВ. ВРЕМЕНИ:осталось %0 замеров,удалены %0 замеров \n").arg(numZam1).arg(ncur));
     ui->textEdit->append(QString("// !\n"));


     ui->doubleSpinBox_45->setValue((double)numZam );
     vctMeas = pPosSolver->mVectBigMeasures;


    ui->doubleSpinBox_46->setValue((double)vctMeas.size() );

    //
    QSolver3D_Angs Solver3D_Angs;
    QSolver2D_Angs Solver2D_Angs;
    QSolver2D_Angs *pSolver2D_Angs;

    if (mTypeOfSolverTask != LBL)
    {
      // по КУ
    for (int i = 0; i < vctMeas.size(); ++i )
    {
       parrMeas[i] =  vctMeas.at(i);
    }



    switch(mTypeOfSolverTask )
        {

        case USBL_2D:
            Solver2D_Angs = QSolver2D_Angs(parrMeas, vctMeas.size(),mtblEstPrfl
                 ,marrBeaconPos[2],0.0009,arrXRez,marrBeaconPos);

            pSolver2D_Angs = &Solver2D_Angs;
            pPosSolver = &Solver2D_Angs;
            ui->doubleSpinBox_47->setVisible(true);

            break;

        case USBL_3D:

          Solver3D_Angs = QSolver3D_Angs(parrMeas, vctMeas.size(),mtblEstPrfl
                                       ,marrBeaconPos[2],0.0009,arrXRez,marrBeaconPos);

        pSolver2D_Angs = &Solver3D_Angs;
        pPosSolver = &Solver3D_Angs;


        ui->doubleSpinBox_47->setVisible(true);
        ui->doubleSpinBox_48->setVisible(true);

            break;

        default:

            break;
        }


    // чистка массива замеров по курс углу
    valTreshold = ui->doubleSpinBox_42->value()* M_PI/ 180.;
    pPosSolver->fncRoughCleaningMeasures(&arrX0[5],valTreshold, 0, &vctGoodZamersNum);

    const int numZam2 = pPosSolver->mVectBigMeasures.size();
    if (numZam2 < 10)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Внимание!");
        msgBox.setText("Отсев по невязке КУ. Хороших замеров менее 10.\n Проверьте настройки !");
        // очистка интерфейса
        ui->pushButtonStartTask->setEnabled(true);
        ui->pushButtonStopTask->setEnabled(false);
        ui->progressBar->setValue(0);
        ui->labelTaskInfo1->setText("");
        ui->labelTaskInfo2->setText("");
        ui->labelTaskInfo3->setText("");
        ui->labelTaskInfo4->setText("");
        ui->labelTaskInfo6->setText("");
        ui->label_2->setText("");
        ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ ПО ОТСЕВУ ПО НЕВ. КУ: ОСТАЛОСЬ МЕНЕЕ 10 ИЗМЕРЕНИЙ. \n"));
        ui->textEdit->append(QString("// !\n"));
        msgBox.exec();
        return;
    }
    ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ОТСЕВА ПО НЕВ.КУ:осталось %0 замеров,удалены %0 замеров \n").arg(numZam2).arg(numZam1 - numZam2));
    ui->textEdit->append(QString("// !\n"));

    ui->doubleSpinBox_47->setValue((double)(pPosSolver->mVectBigMeasures.size()) );
    if (mTypeOfSolverTask == USBL_3D)
    {

     // чистка массива замеров по углу места
        valTreshold = ui->doubleSpinBox_43->value()* M_PI/ 180.;
        pPosSolver->fncRoughCleaningMeasures(&arrX0[5],valTreshold, 1,&vctGoodZamersNum);
       ui->doubleSpinBox_48->setValue((double)(pPosSolver->mVectBigMeasures.size()) );

     //
    }
    const int numZam3 = pPosSolver->mVectBigMeasures.size();
    if (numZam3 < 10)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Внимание!");
        msgBox.setText("Отсев по невязке УМ. Хороших замеров менее 10.\n Проверьте настройки !");
        // очистка интерфейса
        ui->pushButtonStartTask->setEnabled(true);
        ui->pushButtonStopTask->setEnabled(false);
        ui->progressBar->setValue(0);
        ui->labelTaskInfo1->setText("");
        ui->labelTaskInfo2->setText("");
        ui->labelTaskInfo3->setText("");
        ui->labelTaskInfo4->setText("");
        ui->labelTaskInfo6->setText("");
        ui->label_2->setText("");
        ui->textEdit->append(QString("АВАР. ЗАВЕРШЕНИЕ ПО ОТСЕВУ ПО НЕВ. УМ: ОСТАЛОСЬ МЕНЕЕ 10 ИЗМЕРЕНИЙ. \n"));
        ui->textEdit->append(QString("// !\n"));
        msgBox.exec();
        return;
    }
    ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ОТСЕВА ПО НЕВ.УМ:осталось %0 замеров,удалены %0 замеров \n").arg(numZam3).arg(numZam2 - numZam3));
    ui->textEdit->append(QString("// !\n"));
  }


 // построение облаков точек
      double arrMean[3] ={0.},arrDisp[3] ={0.};
     QVector<int> vctNumGood;
     QString qstrOutPutFold0;

     wchar_t wchFolg0[400] = {0};
     if (mqstrOutPutFold.length() >3)
     {
      qstrOutPutFold0 = mqstrOutPutFold + QString("//CloudAfterGeom%0").arg(0);
      qstrOutPutFold0.toWCharArray(wchFolg0);

      wchFolg0[qstrOutPutFold0.length()]=0;
      _wmkdir(wchFolg0);
     }


    QPntCloud :: createPictFiles(wchFolg0,pPosSolver->mVectBigMeasures,mtblEstPrfl
                                ,marrBeaconPos,marrAntPosParams, mTypeOfSolverTask
                                 ,arrMean, arrDisp,&vctNumGood,mTYPE_OF_OUTPUT_FILE
                               );

    ui->tableWidget->item(0,0)->setText(QString::number(sqrt(arrDisp[0]), 'f', 2));
    ui->tableWidget->item(1,0)->setText(QString::number(sqrt(arrDisp[1]), 'f', 2));
    ui->tableWidget->item(2,0)->setText(QString::number(sqrt(arrDisp[2]), 'f', 2));


    mvctMeasures_work = pPosSolver->mVectBigMeasures;
    // !
    myThread->setParams(mvctMeasures_work, mTypeOfSolverTask,vctX0, mtblEstPrfl
        ,marrGPSPosParams, marrBeaconPos[2]
        ,mMaxQuantIter_LBL, mMaxQuantIter_USBL
        ,mGape,mAngStep,  mALGOR_TYPE, mLBLcoeff,mUSBLcoeff
        ,mTYPE_OF_OUTPUT_FILE,mqstrOutPutFold);


   mNumOfTask = 1;

    // запускаем поток (при этом автоматически запустится функция MyThread::run())
    myThread->start();

}
//-----------------------------------------------------
// обработчик нажатия кнопки Stop (принудительная остановка задачи пользователем)
void MainWindow::on_pushButtonStopTask_clicked()
{
    qDebug() << "MainWindow::on_pushButtonStopTask_clicked";
    myThread->userBreak(); // устанавливаем в потоке флаг завершения
}
//-----------------------------------------------------

// обработчик сигнала о текущем прогрессе выполнения задачи
void MainWindow::onTaskProgress(int step, QVector<double>vctNevSquare, QVector<double> vctMean
              , QVector<double> vctDisp, int numPart ,bool bEndPart, QVector<double>vctX0)
{
   // qDebug() << "MainWindow::onTaskProgress" << step << error;
    // вывод в интерфейс

    for(int i =0;i < 8; ++i)
    {
      marrXRez[i] = vctX0.at(i);
    }

   if(numPart == 1)
   {      
       // выгрузка таблицы маяка
       for (int j =0; j < 2; ++j)
       {

        QString strNumber = QString::number(marrXRez[3 + j], 'f', 2);
        ui->tableWidget_7->item(0,j)->setText(strNumber);
       // ui->tableWidget_7->item(1,j)->setText(QString("???"));
       }

       int ia = ((int)(marrBeaconPos[2] * 1000.)) ;
       double a = ((double)ia)/ 1000.;
       ui->tableWidget_7->item(0,2)->setText(QString::number(a));


       // !

       // выгрузка таблицы антенны
       for (int j =0; j < 3; ++j)
       {
           int ia = ((int)(marrXRez[ j] * 1000.)) ;
           double a = ((double)ia)/ 1000.;
        ui->tableWidget_5->item(0,j)->setText(QString::number(a));


       }
   }
   else
   {
      for (int j =0; j < 3; ++j)
       {
           int ia = ((int)(marrXRez[5 + j] * 180./M_PI * 100.)) ;
           double a = ((double)ia)/ 100.;
        ui->tableWidget_5->item(0,3 +j)->setText(QString::number(a));

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

        QString qstr("LBL:");
        qstr += ui->labelTaskInfo2->text();
        qstr += "\n";
        ui->textEdit->append(qstr);

    }
    else
    {
        ui->labelTaskInfo3->setText(QString("q(мр): N= %1, СКО=  %2, Сред= %3, СКЗ= %4").arg(step).arg(sqrt(vctNevSquare.at(0)) *1000., 0, 'f', 3).arg(vctMean.at(0) *1000., 0, 'f', 3).arg(sqrt(vctDisp.at(0)) *1000., 0, 'f', 3));
        QString qstr("USBL:");
        qstr += ui->labelTaskInfo3->text();
        qstr += "\n";
        ui->textEdit->append(qstr);
        if(vctNevSquare.size() == 2)
        {

         ui->labelTaskInfo4->setText(QString("e(мр): N= %1, СКО=  %2, Сред= %3, СКЗ= %4").arg(step).arg(sqrt(vctNevSquare.at(1)) *1000., 0, 'f', 3).arg(vctMean.at(1) *1000., 0, 'f', 3).arg(sqrt(vctDisp.at(1)) *1000., 0, 'f', 3));
         QString qstr1("     ");
         qstr1 += ui->labelTaskInfo4->text();
         qstr1 += "\n";
         ui->textEdit->append(qstr1);
        }
    }

    QString qstr("        arrXRez[8] = ");
    for(int i = 0;i < 8; ++i)
    {
        double coef = 1.;
        if (i > 4)
        {
          coef = 180./M_PI;
        }
        QString strNumber = QString::number(coef *marrXRez[i], 'f', 2);
        qstr += strNumber;
        qstr += "; ";
    }
    qstr += "\n";
    ui->textEdit->append(qstr);
}
//---------------------------------------
// обработчик сигнала о текущем прогрессе выполнения задачи
void MainWindow::onKleaningEnded(int numPart, QVector<int>vctGoodZamersNum, QVector<double>vctMean, QVector<double>vctDisp)
{
    double arrDisp[3] = {0.};
    arrDisp[0] = vctDisp.at(0);
    arrDisp[1] = vctDisp.at(1);
    arrDisp[2] = vctDisp.at(2);



    ui->tableWidget->item(0,numPart)->setText(QString::number(sqrt(arrDisp[0]), 'f', 2));
    ui->tableWidget->item(1,numPart)->setText(QString::number(sqrt(arrDisp[1]), 'f', 2));
    ui->tableWidget->item(2,numPart)->setText(QString::number(sqrt(arrDisp[2]), 'f', 2));
switch(numPart)
{
case 1:
    ui->groupBox_11->setTitle(QString("К-во измерений на входе в алгоритм:"));
    ui->label_2->setText(QString("LBL: %1").arg(vctGoodZamersNum.size()));
    ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ОТСЕВА ПО НЕВЯЗКЕ t,  3 СИГМА:осталось %1 замеров \n").arg(vctGoodZamersNum.size()));
    if(nullptr !=mqstrOutPutFold)
    {
     ui->textEdit->append(QString("СТАТ. ХАРАКТЕРИСТИКИ ОБЛАКА ПОСЛЕ ОТСЕВА \n").arg(vctGoodZamersNum.size()));
     ui->textEdit->append(QString("  Средняя точка облака :X= %1; Y= %2; Z= %3 \n").arg(vctMean.at(0)).arg(vctMean.at(1)).arg(vctMean.at(2)));
     ui->textEdit->append(QString("  СКЗ облака :SigX= %1; SigY= %2; SigZ= %3 \n").arg(sqrt(arrDisp[0])).arg(sqrt(arrDisp[1])).arg(sqrt(arrDisp[2])));
    }
    ui->textEdit->append(QString("// !\n"));
    break;
case 2:
    if(nullptr !=mqstrOutPutFold)
    {
     ui->textEdit->append(QString("СТАТ. ХАРАКТЕРИСТИКИ ОБЛАКА ПОСЛЕ РЕШЕНИЯ ЗАДАЧИ LBL \n").arg(vctGoodZamersNum.size()));
     ui->textEdit->append(QString("  Средняя точка облака :X= %1; Y= %2; Z= %3 \n").arg(vctMean.at(0)).arg(vctMean.at(1)).arg(vctMean.at(2)));
     ui->textEdit->append(QString("  СКЗ облака :SigX= %1; SigY= %2; SigZ= %3 \n").arg(sqrt(arrDisp[0])).arg(sqrt(arrDisp[1])).arg(sqrt(arrDisp[2])));
     ui->textEdit->append(QString("// !\n"));
    }


    break;
case 3:
    ui->textEdit->append(QString("РЕЗУЛЬТАТЫ ОТСЕВА ПО НЕВЯЗКЕ УГЛОВ (ПЕРЕД РЕШЕНИЕМ UISBL),  3 СИГМА:осталось %0 замеров \n").arg(vctGoodZamersNum.size()));
    if(nullptr !=mqstrOutPutFold)
    {
     ui->textEdit->append(QString("СТАТ. ХАРАКТЕРИСТИКИ ОБЛАКА ПОСЛЕ ОТСЕВА \n").arg(vctGoodZamersNum.size()));
     ui->textEdit->append(QString("  Средняя точка облака :X= %1; Y= %2; Z= %3 \n").arg(vctMean.at(0)).arg(vctMean.at(1)).arg(vctMean.at(2)));
     ui->textEdit->append(QString("  СКЗ облака :SigX= %1; SigY= %2; SigZ= %3 \n").arg(sqrt(arrDisp[0])).arg(sqrt(arrDisp[1])).arg(sqrt(arrDisp[2])));
     ui->textEdit->append(QString("// !\n"));
    }

    ui->labelTaskInfo6->setText(QString("USBL: %1").arg(vctGoodZamersNum.size()));
    break;
case 4:
    if(nullptr !=mqstrOutPutFold)
    {
     ui->textEdit->append(QString("СТАТ. ХАРАКТЕРИСТИКИ ОБЛАКА ПОСЛЕ РЕШЕНИЯ ЗАДАЧИ USBL \n").arg(vctGoodZamersNum.size()));
     ui->textEdit->append(QString("  Средняя точка облака :X= %1; Y= %2; Z= %3 \n").arg(vctMean.at(0)).arg(vctMean.at(1)).arg(vctMean.at(2)));
     ui->textEdit->append(QString("  СКЗ облака :SigX= %1; SigY= %2; SigZ= %3 \n").arg(sqrt(arrDisp[0])).arg(sqrt(arrDisp[1])).arg(sqrt(arrDisp[2])));
     ui->textEdit->append(QString("// !\n"));
    }

    break;


default:

    break;
}

}
//---------------------------------------
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
        ui->labelTaskInfo1->setText("Вычисления завершились успешно");
        ui->progressBar->setValue(ui->progressBar->maximum());

        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1 ЗАВЕРШИЛИСЬ УСПЕШНО\n").arg(mNumCkickStart));
        break;
    case MYTHREAD_RESULT_MAX_STEP:
        ui->progressBar->setValue(ui->progressBar->maximum());
        ui->labelTaskInfo1->setText("Ошибка (достигнуто макс. число итер.)");

        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1 ПРОВАЛЕНЫ:\n").arg(mNumCkickStart));
        ui->textEdit->append(QString("Достигнуто максимальное число итераций\n"));
        break;
    case MYTHREAD_RESULT_USER_CANCEL:
        ui->labelTaskInfo1->setText(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1  ПРЕРВАНЫ ПОЛЬЗОВАТЕЛЕМ\n").arg(mNumCkickStart));

        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1 ПРЕРВАНЫ ПОЛЬЗОВАТЕЛЕМ:\n").arg(mNumCkickStart));

        break;
    case MYTHREAD_RESULT_ERROR:
        ui->labelTaskInfo1->setText("Ошибка в исходных данны\n");

        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1 ПРОВАЛЕНЫ:\n").arg(mNumCkickStart));
        ui->textEdit->append(QString("Ошибка в исходных данных\n"));
        break;
    case MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_LBL:
        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %1 ПРОВАЛЕНЫ:\n").arg(mNumCkickStart));
        ui->textEdit->append(QString("На входе LBL осталось менее 10 замеров\n"));
        break;
    case MYTHREAD_RESULT_ERROR_CLEAN_BEFOR_USBL:
        ui->textEdit->append(QString("ВЫЧИСЛЕНИЯ ПОПЫТКИ %0 ПРОВАЛЕНЫ:\n").arg(mNumCkickStart));
        ui->textEdit->append(QString("На входе USBL осталось менее 10 замеров\n"));
        break;
    default:
        ui->labelTaskInfo1->setText("Ошибка (неизвестная причина)\n");
        break;
    }

}



void MainWindow::on_tableWidget_6_cellChanged(int row, int column)
{
   if(!mbtableWidget_6Init)
   {
       return;
   }
  // ui->tableWidget_6->item(row, column)->setText(QString::number(marrBeaconPos[column]));


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

void MainWindow::on_tableWidget_4_cellChanged(int row, int column)
{
    if (!mbtableWidget_4Init)
    {
        return;
    }
    double coef = (column<3)?1.:180./M_PI;

 //  ui->tableWidget_4->item(row, column)->setText(QString::number(marrAntPosParams[column]));


}



void MainWindow::on_pushButton_clicked()
{
    QString str0 = QCoreApplication::applicationDirPath() ;
    mqstrDataFIle = QFileDialog::getOpenFileName(0,"Выбор .acd с исходными данными",str0,"*.acd");
    if(mqstrDataFIle == nullptr)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Внимание!");
        msgBox.setText("Укажите путь к <*.acd> файлу с данными калибровки!!!");
        msgBox.exec();
        return;
    }
    mbOpenDialog_acd = true;

    this->ui->lineEdit_5->setText(mqstrDataFIle);

    if(this->ui->lineEdit_5->text().length() >3)
    {
    readHeaderDataFile(this->ui->lineEdit_5->text(),mQuantMeas0, marrAntPosParams
                                        , marrGPSPosParams, marrBeaconPos, mTObrabotki);
    fillFields();
    ui->doubleSpinBox_41->setValue(mQuantMeas0);
    }

    QString str("PushButton ВЫБОР <*.ACD>:\n " +mqstrDataFIle) ;
    ui->textEdit->append(str);
    ui->textEdit->append("// !");
}

void MainWindow::on_pushButton_4_clicked()
{
    QString str0 = QCoreApplication::applicationDirPath() ;
    mqstrPrflFIle = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости",str0,"*.000");
    if(mqstrPrflFIle == nullptr)
    {
        QMessageBox msgBox;
        msgBox.setWindowTitle("Внимание!");
        msgBox.setText("Укажите путь к <*.000> файлу с верт. профилем скорости звука!!!");
        msgBox.exec();
        return;
    }

    mbOpenDialog_000 = true;


    this->ui->lineEdit_3->setText(mqstrPrflFIle);

    wchar_t wchPrflFIle [400] = {0};
    mqstrPrflFIle.toWCharArray(wchPrflFIle);
    wchPrflFIle[mqstrPrflFIle.length()] = 0;

    QSubWaterBeam::createProfileTbl(wchPrflFIle, &mtblEstPrfl,mType_of_000);



    mZonaR = calcXMCrit( 2., -marrBeaconPos[2], mtblEstPrfl);

    ui->doubleSpinBox_39-> setValue(mZonaR);
    ui->doubleSpinBox_40-> setValue(mZonaR *0.7);
    QString str("PushButton ВЫБОР <*.000>:\n " +mqstrPrflFIle) ;
    ui->textEdit->append(str);
    ui->textEdit->append("// !");


}

//--------------------------------------
// вычисление неизвестного вектора позиционирования методом наименьших квадратов
//arrX[0],arrX[1],arrX[2] - вектор параллакса антенны ПСК
//arrX[3],arrX[4],arrX[5] - вектор координат маяка в ГСК
//arrX[6],arrX[7],arrX[8] - вектор углов ориентации антенны в ПСК
// INPUT:
//arrX0[mDimX] - начальное приближение
// OUTPUT:
//arrXEst[mDimX] - решение
// valNeviaz0 - невязка
/*
bool MainWindow::estimateParams__( QPosSolver *pQPosSolver, double *arrX0,double *arrXRez, double &valNeviaz0)
{
    int iDimX = pQPosSolver->mDimX;
    //0.
   double *arrMtrx_dFGr_po_dX_Inv = new double[iDimX * iDimX];
   memset(arrMtrx_dFGr_po_dX_Inv, 0, iDimX * iDimX * sizeof(double));


    // 1. формирование вектора начальгного приближения
    memcpy(arrXRez, arrX0, iDimX * sizeof(double));

    ///

    // 2. вычисление начальной невязки
    int quantGoodZamers0 = 0;
    valNeviaz0 = pQPosSolver->calcSumNeviazka(arrXRez, quantGoodZamers0);
    ///

    // 3. итреац процесс
   double *arrXCur = new double[iDimX];
   memcpy(arrXCur, arrXRez, iDimX * sizeof(double));

   double valNeviaz1=0.;
   double *arrDel = new double[iDimX];

   memset(arrDel, 0, iDimX* sizeof(double));

   bool breturn = false;
    for (int i =0; i < 400; ++i)
    {
     if(!pQPosSolver->doOneIteration(arrXRez, arrDel,arrMtrx_dFGr_po_dX_Inv))
     {
         delete []arrMtrx_dFGr_po_dX_Inv;
         delete []arrXCur;
         delete []arrDel;
         return false;
     }

     int len = iDimX;
     if (iDimX == 8)
     {
       len = 5;
     }
     double del1 = NormVect(arrDel, len);
     double del2 = (iDimX == 8)?NormVect(&arrDel[5], 3):0.;
        double coef = 0.2;//((i % 3)==0)?0.99:0.5;
        MatrxMultScalar( arrDel, iDimX, 1, coef, arrDel);
        MtrxMinusMatrx(arrXRez, arrDel,iDimX, 1, arrXCur);

        memcpy(arrXRez, arrXCur, iDimX * sizeof(double));
      if ((del1 < 0.015)&&(del2 < 0.0005))
      {

          valNeviaz0 = valNeviaz1;
          breturn = true;
          break;
      }


        valNeviaz0 = valNeviaz1;
    }


    delete []arrMtrx_dFGr_po_dX_Inv;
    delete []arrXCur;
    delete []arrDel;
    return breturn;
}
*/
void MainWindow::on_comboBox_5_currentIndexChanged(int index)
{
    //******* для прогресс бара **************
       // очистка интерфейса
      /* ui->pushButtonStartTask->setEnabled(true);
       ui->pushButtonStopTask->setEnabled(false);
       ui->progressBar->setValue(0);
       ui->labelTaskInfo1->setText("");
       ui->labelTaskInfo2->setText("");
       ui->labelTaskInfo3->setText("");*/
    cleanInteface();
    QString str;
    switch(index)
    {
    case 0:
     ui->groupBox->setTitle(QString("РЕШЕНИЕ LBL ЗАДАЧИ"));
     mTetta_min = 5. / 180. * M_PI;
     mTetta_max = 90. / 180. * M_PI;
     ui->groupBox_tuneUSBL->hide();
     str = QString("LBL");
        break;

    case 1:
     ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-2D ЗАДАЧИ"));
     mTetta_min = 10. / 180. * M_PI;
     mTetta_max = 80. / 180. * M_PI;
     ui->groupBox_tuneUSBL->show();
     str = QString("USBL-2D");
        break;

    case 2:
     ui->groupBox->setTitle(QString("РЕШЕНИЕ USBL-3D ЗАДАЧИ"));
     mTetta_min = 10. / 180. * M_PI;
     mTetta_max = 80. / 180. * M_PI;
     ui->groupBox_tuneUSBL->show();
     str = QString("USBL-3D");
        break;

    default:

        break;

    }

    ui->doubleSpinBox_35->setValue(mTetta_min * 180./ M_PI);
    ui->doubleSpinBox_36->setValue(mTetta_max * 180./ M_PI);

    // в журнал сессии посылка
    ui->textEdit->append("PushButton Антенна: " + str);
    ui->textEdit->append("// !");

}

void MainWindow::setTable(QTableWidget& tbl, double *arr)
{
    for (int i =0; i < tbl.rowCount(); ++i)
        for (int j =0; j < tbl.columnCount() ; ++j)
    {
     int ia = ((int)(arr[i *tbl.columnCount() + j] * 1000.)) ;
     double a = ((double)ia)/ 1000.;
     tbl.item(i,j)->setText(QString::number(a));

    }
}

void MainWindow::setTableZero(QTableWidget& tbl)
{
    for (int i =0; i < tbl.rowCount(); ++i)
        for (int j =0; j < tbl.columnCount() ; ++j)
    {

     tbl.item(i,j)->setText(QString::number(0));

    }
}

void MainWindow::on_comboBox_5_activated(const QString &arg1)
{

}

void MainWindow::on_tableWidget_6_cellClicked(int row, int column)
{

}



void MainWindow::on_comboBox_9_currentIndexChanged(int index)
{
    tuneGroupBoxUSBL(index);
}
//---------------------------
void MainWindow::on_pushButton_5_clicked()
{
    double arr[50]={0.};
    read_tableWidget_7(arr);
            fill_tableWidget_6(arr);


            read_tableWidget_5(arr);
            arr[3] *= M_PI/ 180.;
            arr[4] *= M_PI/ 180.;
            arr[5] *= M_PI/ 180.;
            fill_tableWidget_4(arr);
            ui->textEdit->append("PushButton ВЗЯТЬ ЗА ИД ");
            ui->textEdit->append("// !");

}
//-------------------------------------------------
void MainWindow::read_tableWidget_6()
{
     for (int i=0; i < ui->tableWidget_6->columnCount(); ++i)
    {
          marrBeaconPos[i] =  ui->tableWidget_6->item(0,i)->text().toDouble();

    }
}

//-------------------------

// очистка интерфейса
void MainWindow::cleanInteface()
{
ui->pushButtonStartTask->setEnabled(true);
ui->pushButtonStopTask->setEnabled(false);
ui->progressBar->setValue(0);
ui->labelTaskInfo1->setText("");
ui->labelTaskInfo2->setText("");
ui->labelTaskInfo3->setText("");
ui->labelTaskInfo4->setText("");
ui->labelTaskInfo6->setText("");
ui->label_2->setText("");

switch (ui->comboBox_5->currentIndex())
{
case 0:
    ui->groupBox_ku->setEnabled(false);
    ui->groupBox_um->setEnabled(false);
    break;

case 1:
    ui->groupBox_ku->setEnabled(true);
    ui->groupBox_um->setEnabled(false);
    break;

case 2:
    ui->groupBox_ku->setEnabled(true);
    ui->groupBox_um->setEnabled(true);
    break;

default: break;
}

setTableZero(*(ui->tableWidget_5));
setTableZero(*(ui->tableWidget_7));
setTableZero(*(ui->tableWidget));
}
//-----------------------------------------
QString MainWindow::createInputParamsConstructorString()
{
    QString strOut;
   // QString strOut1(mpwchDataFIle);
    strOut = " 1. Путь к файлу .ACD: \n";
    strOut += mqstrDataFIle +"\n";
    strOut +=" 2. Путь к файлу .000: \n";
    strOut +=mqstrPrflFIle+"\n";
    strOut +=" 3. Путь к директории с графикой: \n";
    strOut +=mqstrOutPutFold+"\n";

    QString str0;
    switch(mTypeOfSolverTask)
    {
    case LBL:
        mTypeOfSolverTask = LBL;
     str0 = "LBL";
        break;

    case USBL_2D:
        str0 = "USBL_2D";
        break;

    case USBL_3D:
       str0 = "USBL_3D";
        break;

    default:
        break;
    }
    strOut +=QString(" 4. Тип антенны: ") + str0 + "\n";

    switch(ui->comboBox_6->currentIndex())
    {
    case 0:
        str0 = "Реальный";
        break;
    case 1:
        str0 = "Средний";
        break;
    default:
        break;

    }
    strOut +=QString(" 5. Тип профиля звука: ") + str0 + "\n";



    return strOut;
}
//----------------------------------------
QString MainWindow::createInputParamsString(int num)
{
    QString strOut = QString("НАЧАЛО ВЫЧИСЛИТЕЛЬНОГО ЦИКЛА №%0. ВХОД ИНФОРМАЦИЯ.\n ").arg(num);

    strOut += " 1. Путь к файлу .ACD: \n";
    strOut += mqstrDataFIle +"\n";
    strOut += "\\ !\n";
    strOut +=" 2. Путь к файлу .000: \n";
    strOut +=mqstrPrflFIle+"\n";
    strOut += "\\ !\n";
    strOut +=" 3. Путь к директории с графикой: \n";
    strOut +=mqstrOutPutFold+"\n";
    strOut += "\\ !\n";

    // тип антенны
    QString str0;
    switch(mTypeOfSolverTask)
    {
    case LBL:
        mTypeOfSolverTask = LBL;
     str0 = "LBL";
        break;

    case USBL_2D:
        str0 = "USBL_2D";
        break;

    case USBL_3D:
       str0 = "USBL_3D";
        break;

    default:
        break;
    }
    strOut +=QString(" 4. Тип антенны: ") + str0 + "\n";
    strOut += "\\ !\n";

    // тип используемого профиля
    switch(ui->comboBox_6->currentIndex())
    {
    case 0:
        str0 = "Реальный";
        break;
    case 1:
        str0 = "Средний";
        break;
    default:
        break;

    }
    strOut +=QString(" 5. Тип профиля звука: ") + str0 + "\n";
    strOut += "\\ !\n";

   // путь к папке с графикой
    strOut +=QString(" 6. Путь к папке с графикой: ") + mqstrOutPutFold + "\n";
    strOut += "\\ !\n";

    // 7. расширение графических файлов
     QString qstr0 = (ui->comboBox->currentIndex()==0)? QString("<*.SHP>"):QString("<*.CSV>");
     strOut +=QString(" 7. расширение графических файлов: ") + qstr0  + "\n";
     strOut += "\\ !\n";


     // 8. макс дальность прямой видимости
      strOut +=QString(" 8. макс дальность прямой видимости = %1 м\n").arg(mZonaR) ;
      strOut += "\\ !\n";

      // 9. исходное количество замеров
      int numZam0= (int)(ui->doubleSpinBox_41->value());
       strOut +=QString(" 9. исходное количество замеров = %0 \n").arg(numZam0) ;
       strOut += "\\ !\n";

       // 10. априорная позиция маяка:
        strOut +=QString("10. априорная позиция маяка: X= %2 м; Y= %2 м; Z= %2 м \n").arg(marrBeaconPos[0]).arg(marrBeaconPos[1]).arg(marrBeaconPos[2]) ;
        strOut += "\\ !\n";

        // 11. позиция GPS:
         strOut +=QString("11. позиция GPS: X= %2 м; Y= %2 м; Z= %2 м \n").arg(marrGPSPosParams[0]).arg(marrGPSPosParams[1]).arg(marrGPSPosParams[2]) ;
         strOut += "\\ !\n";

         // 12. априорный офсет антенны:
          strOut +=QString("12. априорный офсет антенны: X= %2 м; Y= %2 м; Z= %2 м \n").arg(marrAntPosParams[0]).arg(marrAntPosParams[1]).arg(marrAntPosParams[2]) ;
          strOut += "\\ !\n";

          // 13. априорные углы Эйлера антенны:
           strOut +=QString("13. априорные углы Эйлера антенны(град): bet= %2; eps= %2; alf= %2\n").arg(marrAntPosParams[3] * 180./ M_PI).arg(marrAntPosParams[4]* 180./ M_PI).arg(marrAntPosParams[5]* 180./ M_PI) ;
           strOut += "\\ !\n";


           // 14. параметры предварительного отсева по геометрии:
            strOut +=QString("14. параметры предварительного отсева по геометрии:\n");
            strOut +=QString("   Дгор<= %1 м; %3 град <= УМ <=%3 град\n").arg(mWorkR).arg(mTetta_min* 180./ M_PI).arg(mTetta_max* 180./ M_PI);
            strOut += "\\ !\n";


            // 15. параметры предварительного отсева по невязке:
            double valTreshold_t = ui->doubleSpinBox_44->value();
             double valTreshold_q = ui->doubleSpinBox_42->value();
              double valTreshold_e = ui->doubleSpinBox_43->value();
             strOut +=QString("15. параметры предварительного отсева по невязке:\n");
             strOut +=QString("   dt<= %1 мc; dКУ <=%1 град; dУМ <=%1 град\n").arg(valTreshold_t).arg(valTreshold_q).arg(valTreshold_e);
             strOut += "\\ !\n";

             // 16. параметры настройки алгоритма LBL:
              strOut +=QString("16. параметры настройки алгоритма LBL: \n");
              strOut +=QString("   Коеф= %4 ; Макс. ит= %0 град\n").arg(mLBLcoeff).arg(mMaxQuantIter_LBL);
              strOut += "\\ !\n";

              // 17. параметры настройки алгоритма USBL:
              strOut +=QString("17. параметры настройки алгоритма USBL:\n");
              switch(ui->comboBox_9->currentIndex())
              {
              case 0:
                  strOut +=QString("   Ручной режим (перебор): Gap= %2 град; Step= %3 град\n").arg(mGape).arg(mAngStep);
                  break;

              case 1:
                  strOut +=QString("   Метод Ньютона: Коеф= %4 ; Макс. ит= %0 \n").arg(mUSBLcoeff).arg(mMaxQuantIter_USBL);
                  break;

              default: break;
              }
            strOut += "\\ !\n";
               strOut += "'\\' /*КОНЕЦ БЛОКА ВХОДНАЯ ИНФОРМАЦИЯ*/!\n";
               strOut += "\\************************************************************************* !\n";
               strOut += "\\***** НАЧАЛО ВЫЧИСЛЕНИЙ ************************************************ !\n";


    return strOut;
}




void MainWindow::on_comboBox_6_activated(const QString &arg1)
{

}

void MainWindow::on_comboBox_6_currentIndexChanged(int index)
{

}

void MainWindow::on_comboBox_6_currentIndexChanged(const QString &arg1)
{
     ui->textEdit->append("ComboBox Profile: " + arg1);
     ui->textEdit->append("// !");
}

void MainWindow::on_comboBox_currentIndexChanged(const QString &arg1)
{
    ui->textEdit->append("ComboBox ТИП РАСШИРЕНИЯ ФАЙЛА ГРАФИКИ: " + arg1);
    ui->textEdit->append("// !");
}

void MainWindow::on_comboBox_9_currentTextChanged(const QString &arg1)
{
    ui->textEdit->append("ComboBox НАСТРОЙКА USBL: " +arg1);
    ui->textEdit->append("// !");
}

void MainWindow::on_comboBox_activated(const QString &arg1)
{

}

void MainWindow::on_pushButton_6_clicked()
{
    QString str0 = QCoreApplication::applicationDirPath() ;
    if (mqstrOutPutFold != nullptr)
    {
      str0 =mqstrOutPutFold;
    }
    str0 += QString("\\report.txt");
    QFile fileOut(str0); // Связываем объект с файлом str0

    if(fileOut.open(QFile::WriteOnly |
                    QFile::Text))
    { // Если файл успешно открыт для записи в текстовом режиме

        QTextStream writeStream(&fileOut); // Создаем объект класса QTextStream
// и передаем ему адрес объекта fileOut
      //  writeStream << "Text, text, text."; // Посылаем строку в поток для записи

        QString ss = ui->textEdit->toPlainText();
        writeStream<<ss;

        fileOut.flush();
        fileOut.close(); // Закрываем файл
    }

}

void MainWindow::on_pushButton_NEW_clicked()
{
    // README!    ВАЖНО!
    /* Предварительно надо выбрать файл .acd
     * информация из хидера вставляется в таблицы интерфейса
     */

    // 1. чтение пути к файлу профиля звука
    mqstrPrflFIle = this->ui->lineEdit_3->text();
    wchar_t wchPrflFIle [400] = {0};
    mqstrPrflFIle.toWCharArray(wchPrflFIle);
    wchPrflFIle[mqstrPrflFIle.length()] = 0;

    // !1

    // 2. директория с результатом
    mqstrOutPutFold = this->ui->lineEdit->text();
    wchar_t wchOutPutFold[400];
    mqstrOutPutFold.toWCharArray(wchOutPutFold);
    wchOutPutFold[mqstrOutPutFold.size()] = 0;
    // !2

    // 3.чтение пути к файлу с исходными данными .acd
    mqstrDataFIle = this->ui->lineEdit_5->text();
    wchar_t wchDataFIle [400] = {0};
    mqstrDataFIle.toWCharArray(wchDataFIle);
    wchDataFIle[mqstrDataFIle.length()] = 0;

    // !3


   // 4. создание профиля скорости звука
      //wchPrflFIle - путь к файлу wchar_t*
    // mtblEstPrfl - таблица скорости звука
    // mType_of_000 - тип фала .000, равен VAR0
    //
   QSubWaterBeam::createProfileTbl(wchPrflFIle, &mtblEstPrfl,mType_of_000);
  // 4!



    // присвоение относительных координат антенны (6 мерный массив - вектиор офсета, 3 угла)
  for (int i = 0; i < 6; ++i)
  {
      if (i <3)
      {
          marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();

      }
      else
      {
          marrAntPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;

      }
  }

  // присвоение координат жпс
  read_tableWidget_8(marrGPSPosParams);

  // присвоение координат маяка
  read_tableWidget_6(marrBeaconPos);
  ///






    //чтение файла с исходными данными
   QString qstr = this->ui->lineEdit_5->text();

    QVector <QBigMeasure> vctMeas;

    // чтение массива с измерениями из файла .acd в QVector vctMeas
    readBigMeasuresOnly_acd(qstr,marrAntPosParams
                            , marrGPSPosParams, marrBeaconPos, mTObrabotki,&vctMeas);
    //буфер вывода
    int qrows=31;
    double *parrBuff = new double [qrows * vctMeas.size()];
    memset (parrBuff, 0, sizeof(double) *qrows * vctMeas.size());



    // основной цикл
    for (int i =0;i <vctMeas.size(); ++i )
    {
        // 1. работаем с одним измерением. вытаскиваем его и обрабатываем.
        QBigMeasure meas = vctMeas.at(i);
        // !1

        // 2.  вычисление вектора положения антенны в момент приема в ГСК - arrSAnt_GSK[3]  3 ++

        double arrV1[3] = {0.}, arrSAnt_KGSK[3] = {0.}, arrSAnt_GSK[3] = {0.};
        QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,meas.marrMuZv
                                                                ,marrAntPosParams, arrSAnt_KGSK);
        MtrxSumMatrx(meas.marrSVessZv, arrSAnt_KGSK,1, 3, arrSAnt_GSK) ;
        // 2!

        // 3. вычисление вектора положения антенны в момент запроса в ГСК - arrSAntWave_GSK[3]   3 ++
        double  arrSAntWave_KGSK[3] = {0.}, arrSAntWave_GSK[3] = {0.};
        QPeaceVess::recalcPositionFromGdgSphericalSK_to_KGSK(arrV1,meas.marrMuWaveZv
                                                                ,marrAntPosParams, arrSAntWave_KGSK);
        MtrxSumMatrx(meas.marrSVessWaveZv, arrSAntWave_KGSK,1, 3, arrSAntWave_GSK) ;
        // 3!

        // 4. расстояние от антенны до маяка в момент запроса - val_Lp1       1  +
        double arrt[3] = {0.};
        MtrxMinusMatrx(arrSAntWave_GSK, marrBeaconPos,1, 3, arrt);
        double val_Lp1 = Norm3( arrt) ;
        // 4!

        // 5. расстояние от антенны до маяка в момент приема - val_Lp2        1  ++
        MtrxMinusMatrx(arrSAnt_GSK, marrBeaconPos,1, 3, arrt);
        double val_Lp2 = Norm3( arrt) ;
        // 5!

        // 6.  valLs1 - длина кривой, по которой должен пройти сигнал запроса от антенны к маяку    1

            // 6.1 вычисление горизонтального расстояния от антенны до маяка на момент запроса
        double valBeaconHorDist = sqrt((arrSAntWave_GSK[0] - marrBeaconPos[0]) * (arrSAntWave_GSK[0] - marrBeaconPos[0])
                + (arrSAntWave_GSK[1] - marrBeaconPos[1]) * (arrSAntWave_GSK[1] - marrBeaconPos[1]));
             // ! 6.1

              // 6.2 вычисление угла скольжения луча на момент запроса - valTetta_request
        double valTetta_request = 0.;
        bool brez = calcTetta( -arrSAntWave_GSK[2],-marrBeaconPos[2],valBeaconHorDist
                       , mtblEstPrfl, valTetta_request);
              // 6.3 вычисление длины дуги луча
        double valLs1 = calc_CurveLength( -arrSAntWave_GSK[2],-marrBeaconPos[2]
                                           ,valTetta_request, mtblEstPrfl);
             // 6.3!


        // 7.  valLs2 - длина кривой, по которой должен пройти сигнал ответа от антенны к маяку      1

            // 7.1 вычисление горизонтального расстояния от антенны до маяка на момент запроса
        double valBeaconHorDist1 = sqrt((arrSAnt_GSK[0] - marrBeaconPos[0]) * (arrSAnt_GSK[0] - marrBeaconPos[0])
                + (arrSAnt_GSK[1] - marrBeaconPos[1]) * (arrSAnt_GSK[1] - marrBeaconPos[1]));
             // ! 7.1

              // 7.2 вычисление угла скольжения луча на момент запроса - valTetta_request
        double valTetta_response = 0.;
        brez = calcTetta( -arrSAnt_GSK[2], -marrBeaconPos[2],valBeaconHorDist
                       , mtblEstPrfl, valTetta_response);
              // 7.3 вычисление длины дуги луча
        double valLs2 = calc_CurveLength( -arrSAnt_GSK[2],-marrBeaconPos[2]
                                           ,valTetta_response, mtblEstPrfl);
             // 7.3!

        // 8.t, КУ, УМ - измерения, которые в теории должна была измерить антенна       3

        // Преобразование координат точки из ГСК в Антенную сферическую систему координат
        //INPUT:
        //tblPrfl - профиль скорости звука
        //arrGSK_XYZ[3] - положение точки в ГСК
        // arrSVessGSK[3] - положение корабля в ГСК
        //arrMu[3] - вектор палубных углов
        //arrAntPosParams[6] - вектор параметров позиционирования антенны
        //OUTPUT:
        //pval_q - курсовой угол сигнала в АСфСК
        //pval_e - угол места сигнала в АСфСК
        //pval_t - время прохождения сигнала

        double val_q = 0., val_e = 0., val_t = 0.;
        bool brez1 =  transf_GSK_XYZ_to_USBL3D(mtblEstPrfl, marrBeaconPos, meas.marrSVessZv, meas.marrMuZv
          ,marrAntPosParams,  &val_q, &val_e,  &val_t);

        //


        // заполняем строку буфера вывода
        meas.fill_info_row(marrGPSPosParams,&parrBuff[qrows * i]);
      double *p = &parrBuff[qrows * i + 16];
      memcpy(p, arrSAnt_GSK, 3 * sizeof(double));

      memcpy(&p[3], arrSAntWave_GSK, 3 * sizeof(double));

      p[6] = val_Lp1;

      p[7] = val_Lp2;

      p[8] = valLs1;

      p[9] = valLs2;

      p[10] = val_q *180./M_PI ;
      p[11] = val_e *180./M_PI ;
      p[12] = val_t ;

      // угол скольжения расчетный
      p[13] =  valTetta_response*180./M_PI ;

      // вычисление угла скольжения измеренного
      double val_KGSK_q = 0., val_KGSK_e = 0.;
      calcSphericalAnglesKGSK(mtblEstPrfl, meas.mqzv, meas.mezv, meas.marrMuZv
             ,marrAntPosParams, &val_KGSK_q, &val_KGSK_e);

       p[14] =  -val_KGSK_e*180./M_PI ;
    }
 // вывод в CSV файл
    wchar_t wchBuff[500] = {0};
    QString  qstrBuff = mqstrOutPutFold + QString("//Buff.csv").arg(0);
    qstrBuff.toWCharArray(wchBuff);
    wchBuff[qstrBuff.length()]=0;


    int iNumCols = qrows;
    int iLenName = 30;


    wchar_t *pwcharrColNames = new wchar_t[qrows* 30];
    memset(pwcharrColNames, 0,qrows* 30  * sizeof(wchar_t));
    wcscpy( pwcharrColNames,  L"TZapr");

    wcscpy( &pwcharrColNames[iLenName],  L"SzaprX");
    wcscpy( &pwcharrColNames[iLenName *2],  L"SzaprY");
    wcscpy( &pwcharrColNames[iLenName *3],  L"SzaprZ");

    wcscpy( &pwcharrColNames[iLenName *4],  L"ZaprCurse_0");
    wcscpy( &pwcharrColNames[iLenName *5],  L"ZaprPitch");
    wcscpy( &pwcharrColNames[iLenName *6],  L"ZaprRoll");

    wcscpy( &pwcharrColNames[iLenName *7],  L"TOtv");

    wcscpy( &pwcharrColNames[iLenName *8],  L"SotvX");
    wcscpy( &pwcharrColNames[iLenName *9],  L"SotvY");
    wcscpy( &pwcharrColNames[iLenName *10],  L"SotvZ");

    wcscpy( &pwcharrColNames[iLenName *11],  L"OtvCurse");
    wcscpy( &pwcharrColNames[iLenName *12],  L"OtvPitch");
    wcscpy( &pwcharrColNames[iLenName *13],  L"OtvRoll");

    wcscpy( &pwcharrColNames[iLenName *14],  L"qZv");
    wcscpy( &pwcharrColNames[iLenName *15],  L"eZv");

    wcscpy( &pwcharrColNames[iLenName *16],  L"PriemAntX");
    wcscpy( &pwcharrColNames[iLenName *17],  L"PriemAntY");
    wcscpy( &pwcharrColNames[iLenName *18],  L"PriemAntZ");

    wcscpy( &pwcharrColNames[iLenName *19],  L"ZaprAntX");
    wcscpy( &pwcharrColNames[iLenName *20],  L"ZaprAntY");
    wcscpy( &pwcharrColNames[iLenName *21],  L"ZaprAntZ");

    wcscpy( &pwcharrColNames[iLenName *22],  L"Lp1");
    wcscpy( &pwcharrColNames[iLenName *23],  L"Lp2");

    wcscpy( &pwcharrColNames[iLenName *24],  L"Ls1");
    wcscpy( &pwcharrColNames[iLenName *25],  L"Ls2");

    wcscpy( &pwcharrColNames[iLenName *26],  L"val_q");
    wcscpy( &pwcharrColNames[iLenName *27],  L"val_e");
    wcscpy( &pwcharrColNames[iLenName *28],  L"val_t");
    wcscpy( &pwcharrColNames[iLenName *29],  L"Tetta");
    wcscpy( &pwcharrColNames[iLenName *30],  L"TettaZv");

    TYrWrite::WriteMassiveInFIleSCV(wchBuff,parrBuff, vctMeas.size(), iNumCols
                                 ,NULL,pwcharrColNames, iLenName);

    delete []parrBuff;
    delete []pwcharrColNames;
}
