#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <float.h>
#include <wchar.h>
#include <math.h>
#include <string.h>
#include <QFileDialog>
#include <QMessageBox>

#include <stdlib.h>

#include "YrRead.h"
#include "Plane.h"



#include "YrWriteShapeFile.h"

#include "dir.h"
#include "Comp.h"
#include "Equations.h"

#include "MatrixProccess.h"
#include "YrWrite.h"
#include "LinDiffEq.h"

#include "LinOptCtrlSyst_dim21.h"
#include "LinOptCtrlSyst.h"
#include "UrPointXY.h"
#include "URPolyLine.h"
#include "Brls.h"
#include "Segment.h"
#include "Lnsgm.h"

#include "Ctrl.h"
#include "StatSolutionParams.h"

#include "DriverShipImit.h"
#include "TargTrackingOperation.h"
#include "Device.h"

#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsReport0;

extern const double CONST_ZP ;
extern const double GEARS[8];

extern const double VAL_BRLS_LENGTH ;
extern const double VAL_BRLS_HEIGHT ;
extern const double VAL_BRLS_MASS ;
extern const double VAL_BRLS_Cv ;
extern const double VAL_BRLS_dOm_po_dt;
extern const double VAL_BRLS_Cx;

extern const double VAL_TEST_Cx;
extern const double VAL_TEST_J;
extern const double VAL_TEST_Cv ;
extern const double VAL_TEST_dOm_po_dt;
extern const double VAL_TEST_Cx;

//extern const bool BEZ_SHUMOV = true;
extern  bool BEZ_SHUMOV = true;
//шаг фильрации

extern  double constDbl_mh = 1./20000.;

double VAlIntegrStep = 1. / 20000.;

extern const double CONST_J_ROTOR;

extern  bool bIDEAL_FILT = false;



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

     ui->setupUi(this);

     mbtableWidget_2Init = false;
     mbtableWidget_Init = false;
     mbtableWidget_4Init = false;
    // ui->comboBox->setCurrentIndex(0);
     // создание таблицы параметров
        ui->tableWidget->setColumnCount(8);
        ui->tableWidget->setRowCount(2 );


        QHeaderView *pHorHeader = ui->tableWidget->horizontalHeader();
        pHorHeader->setVisible(true);
        QStringList lst;
        lst<< "L, мГ"<<"R, Ом"<<"Psi_f"<<"J_load"<<"C_om"<<"Cv"<<"AmpMResid"<<"PhMResid";


        ui->tableWidget->setHorizontalHeaderLabels(lst);
        QHeaderView *pVertHeader = ui->tableWidget->verticalHeader();
        pVertHeader->setVisible(true);


        for (int i=0; i<  ui->tableWidget->rowCount() ; i++)
                     for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                     {
                         QTableWidgetItem* ptwi0 = nullptr;
                         ptwi0 = new QTableWidgetItem(QString::number(0.));
                         ui->tableWidget->setItem(i,j,ptwi0);
                     }
        ui->tableWidget->item(0,0)->setText(QString::number(0.0006* 1000.));
        ui->tableWidget->item(0,1)->setText(QString::number(0.2));
        ui->tableWidget->item(0,2)->setText(QString::number(0.066));


        ui->tableWidget->item(0,3)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(0,4)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(0,5)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(0,6)->setText(QString::number(0.));
        ui->tableWidget->item(0,7)->setText(QString::number(0.));



        ui->tableWidget->item(1,0)->setText(QString::number(0.0006* 1000.));
        ui->tableWidget->item(1,1)->setText(QString::number(0.2));
        ui->tableWidget->item(1,2)->setText(QString::number(0.066));

        ui->tableWidget->item(1,3)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(1,4)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(1,5)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(1,6)->setText(QString::number(0.9 * 0.8));
        ui->tableWidget->item(1,7)->setText(QString::number(0.));

        mbtableWidget_Init  = true;
///
      // double arr[8] = {-11.,-12.,-14.,0.
       //                 ,0.,0.,0.,0.};
   //init_tableWidget_2(4, arr);

       ui->tableWidget_2->setColumnCount(4);
       ui->tableWidget_2->setRowCount(2 );
       QHeaderView *pHorHeader2 = ui->tableWidget_2->horizontalHeader();
       pHorHeader2->setVisible(true);
   QHeaderView *pVertHeader2 = ui->tableWidget_2->verticalHeader();
       pVertHeader2->setVisible(true);
    QStringList lst1;
        lst1<< "Lamb1"<<"Lamb2"<<"Lamb3"<<"Lamb4";
       ui->tableWidget_2->setHorizontalHeaderLabels(lst1);
    QStringList lst2;
        lst2<< "Re"<<"Im";
       ui->tableWidget_2->setVerticalHeaderLabels(lst2);



      QTableWidgetItem* ptwi0 = nullptr;
       QTableWidgetItem* ptwi01 = nullptr;
       QTableWidgetItem* ptwi02 = nullptr;
       QTableWidgetItem* ptwi03 = nullptr;
       QTableWidgetItem* ptwi13 = nullptr;
       QString qstr0 = QString::number(0.);
       ptwi0 = new QTableWidgetItem(qstr0);
       ptwi01 = new QTableWidgetItem(qstr0);
       ptwi02 = new QTableWidgetItem(qstr0);
       ptwi03 = new QTableWidgetItem(qstr0);
       ptwi13 = new QTableWidgetItem(qstr0);
       ui->tableWidget_2->setItem(1,0,ptwi0);
       ui->tableWidget_2->setItem(1,1,ptwi01);
       ui->tableWidget_2->setItem(1,2,ptwi02);
       ui->tableWidget_2->setItem(1,3,ptwi03);

       QTableWidgetItem* ptwi1 = nullptr;
       QString qstr1 = QString::number(-100.);
       ptwi1 = new QTableWidgetItem(qstr1);
       ui->tableWidget_2->setItem(0,0,ptwi1);

       QTableWidgetItem* ptwi2 = nullptr;
       QString qstr2 = QString::number(-101.);
       ptwi2 = new QTableWidgetItem(qstr2);
       ui->tableWidget_2->setItem(0,1,ptwi2);

       QTableWidgetItem* ptwi3 = nullptr;
       QString qstr3 = QString::number(-102.);
       ptwi3 = new QTableWidgetItem(qstr3);
       ui->tableWidget_2->setItem(0,2,ptwi3);

       QTableWidgetItem* ptwi30 = nullptr;
       QString qstr30 = QString::number(-103.);
       ptwi30 = new QTableWidgetItem(qstr30);
       ui->tableWidget_2->setItem(0,3,ptwi30);



       ///



       ui->tableWidget_3->setColumnCount(1);
       ui->tableWidget_3->setRowCount(2 );
       ui->tableWidget_3->horizontalHeader()->setVisible(true);
       //pHorHeader2->setVisible(true);
       ui->tableWidget_3->verticalHeader()->setVisible(true);

       QStringList lst3;
       lst3<< "Нач. условия";
       ui->tableWidget_3->setHorizontalHeaderLabels(lst3);
       QStringList lst4;
       lst4<< "Tetta"<<"Omega";
       ui->tableWidget_3->setVerticalHeaderLabels(lst4);


       QTableWidgetItem* ptwi = nullptr;
       QTableWidgetItem* ptwi00 = nullptr;

       QString qstr = QString::number(0.);



    ptwi = new QTableWidgetItem(qstr);
    ptwi00 = new QTableWidgetItem(qstr);
     ui->tableWidget_3->setItem(0,0,ptwi);
     ui->tableWidget_3->setItem(1,0,ptwi00);
     //---------------------------------------------

  ui->comboBox_3->setCurrentIndex(0);


  ui->tableWidget_6->setColumnCount(1);
  ui->tableWidget_6->setRowCount(3);


  QTableWidgetItem* ptwi20 = nullptr;
  QString qstr20 = QString::number(0.);
  ptwi20 = new QTableWidgetItem(qstr20);
  ui->tableWidget_6->setItem(0,0,ptwi20);

  QTableWidgetItem* ptwi21 = nullptr;
  QString qstr21 = QString::number(180.);
  ptwi21 = new QTableWidgetItem(qstr21);
  ui->tableWidget_6->setItem(1,0,ptwi21);

  QTableWidgetItem* ptwi31 = nullptr;
  QString qstr31 = QString::number(1000.);
  ptwi31 = new QTableWidgetItem(qstr31);
  ui->tableWidget_6->setItem(2,0,ptwi31);

  ui->tableWidget_6->item(0,0)->setText(QString::number(0.));
  ui->tableWidget_6->item(1,0)->setText(QString::number(90.));
  ui->tableWidget_6->item(2,0)->setText(QString::number(1000.));

  mDriveMeasImitator = QMeasurmentImitator ();

  mFiltr = QFiltr();

  ui->comboBox_4->setCurrentIndex(1);
  switch(ui->comboBox_4->currentIndex())
  {
  case 0:
      ui->frame->setVisible(false);
      ui->tableWidget_6->setVisible(true);
      break;
  case 1:
      ui->frame->setVisible(true);
      ui->tableWidget_6->setVisible(false);
      ui->doubleSpinBox_15->setVisible(false);
      break;
  case 2:
      ui->frame->setVisible(true);
      ui->tableWidget_6->setVisible(false);
      ui->doubleSpinBox_15->setVisible(true);
  default:

      break;

  }

  // установка параметров позиционирования
  memset(marrOecPosParams, 0, 6 * sizeof(double));
  memset(marrOecPosParams_Zv, 0, 6 * sizeof(double));

  marrOecPosParams[0] = marrOecPosParams_Zv[0] = 15.;
  marrOecPosParams[1] = marrOecPosParams_Zv[1] = 15.;
  marrOecPosParams[2] = marrOecPosParams_Zv[2] = 15.;

  double arr[18] = {0.};
  for (int j =0; j < 6; ++j)
  {
      arr[j] = marrOecPosParams[j];
      arr[12 + j] = marrOecPosParams_Zv[j];
  }


  mbtableWidget_4Init = false;
  for (int i=0; i<  ui->tableWidget_4->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = nullptr;
                     ptwi0 = new QTableWidgetItem(QString::number(0.));
                     ui->tableWidget_4->setItem(i,j,ptwi0);
                 }
  //write_tableWidget_4(arr);


  for (int j =0; j < ui->tableWidget_4->columnCount(); ++j)
  {
      int ia = 10000.* marrOecPosParams[j];
      double temp = ia/10000.;
   ui->tableWidget_4->item(0,j)->setText(QString::number(temp));

   ia = 10000.* marrOecPosParams_Zv[j];
   temp = ia/10000.;
   ui->tableWidget_4->item(2,j)->setText(QString::number(temp));
  }

  mbtableWidget_4Init = true;
  ///

  ui->tableWidget_5->setColumnCount(3);
  ui->tableWidget_5->setRowCount(2 );
  ui->tableWidget_5->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_5->verticalHeader()->setVisible(true);

  QStringList lst5;
  lst5<<"Сист. ош."<<"СКЗ"<<"Сумм. СКО";
  ui->tableWidget_5->setHorizontalHeaderLabels(lst5);
  QStringList lst6;
  lst6<< "Тетта"<<"Омега";
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

  ui->tableWidget_7->setColumnCount(3);
  ui->tableWidget_7->setRowCount(2 );
  ui->tableWidget_7->horizontalHeader()->setVisible(true);
  //pHorHeader2->setVisible(true);
  ui->tableWidget_7->verticalHeader()->setVisible(true);

  QStringList lst7;
  lst5<<"Сист. ош."<<"СКЗ"<<"Сумм. СКО";
  ui->tableWidget_7->setHorizontalHeaderLabels(lst7);
  QStringList lst8;
  lst6<< "Тетта"<<"Омега";
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
 constDbl_mh = VAlIntegrStep = 1./20000;
 ui->doubleSpinBox_27->setValue(1. /constDbl_mh);

 if(ui->checkBox_2->isChecked())
 {
     bIDEAL_FILT = true;
 }
 else
 {
     bIDEAL_FILT = false;
 }
}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------


void MainWindow:: inputData()
{
    // начальное время
    mT0 = 0.;
    ///
    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;



    // момент сухого трения
      mMomDryFriction = ui->doubleSpinBox_28->value();
    // МАТЕМАТИЧЕСКАЯ МОДЕЛЬ
      mModelInductL =   (ui->tableWidget->item(0,0)->text().toDouble())/ 1000.;

      mModelResist = ui->tableWidget->item(0,1)->text().toDouble();

      mModelPsi_f = ui->tableWidget->item(0,2)->text().toDouble();

      //мом инерции ротора
        mModelJ0= CONST_J_ROTOR;

      // момент сопротивления эл двигателя при обесточенных обмотках, не более
       mModelMaxAmpMomResidual = 0.;// ui->tableWidget->item(0,7)->text().toDouble();
       mModelPhMomResidual = 0.;// ui->tableWidget->item(0,8)->text().toDouble()/ 180. * M_PI;
       //mModelMomOut = ui->tableWidget->item(0,3)->text().toDouble();

  // НАГРУЗКА
      //мом инерции нагрузки
        mModelJPayLoad= ui->tableWidget->item(0,3)->text().toDouble();
      // коэффициент сопротивления нагрузки
        mModelCx= ui->tableWidget->item(0,4)->text().toDouble();
        mModelCv= ui->tableWidget->item(0,5)->text().toDouble();

        // создание нагрузки
        int iii = ui->comboBox_3->currentIndex();
        if(ui->comboBox_3->currentIndex()==0) // нагрузка тестовая
        {
            mHorLoadModel = QLoad ( mModelJPayLoad , mModelCx,mModelCv, VAL_TEST_dOm_po_dt);
            mpHorLoadModel = &mHorLoadModel;

            mVertLoadModel = QLoad ( mModelJPayLoad , mModelCx,mModelCv, VAL_TEST_dOm_po_dt);
            mpVertLoadModel = &mVertLoadModel;
        }

        if(ui->comboBox_3->currentIndex()==1)
        {
         QSegment SgmBarrier;
          TPlane Plane;
          mHorBrlsModel =  QBrls ( mModelCv, mModelCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                              ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane);
          mpHorLoadModel = &mHorBrlsModel;

          mVertBrlsModel =  QBrls ( mModelCv, mModelCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                              ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane);
          mpVertLoadModel = &mVertBrlsModel;
        }
       //


        marrSpreadParams[0] = mModelInductL * 0.02;
        marrSpreadParams[1] = mModelResist * 0.005;
        marrSpreadParams[2] = mModelPsi_f * 0.03;
        marrSpreadParams[3] = mpHorLoadModel->mJPayLoad * 0.05;
        marrSpreadParams[4] = mpHorLoadModel->mCx *0.05;
        marrSpreadParams[5] = mpHorLoadModel->mCv *0.05;

        double valModelMomDryFriction = 0.;
        mModelOecHorElectMotor  = QElectMotor (mModelInductL , mModelResist, mModelPsi_f, mModelJ0
                ,mModelMaxAmpMomResidual, mModelPhMomResidual
                       ,0.,0.,valModelMomDryFriction);
        mModelOecHorDriver  = QElectDriver(mModelOecHorElectMotor  , mpHorLoadModel, marrSpreadParams );

        mModelOecVertElectMotor  = QElectMotor (mModelInductL , mModelResist, mModelPsi_f, mModelJ0
                ,mModelMaxAmpMomResidual, mModelPhMomResidual
                       ,0.,0.,valModelMomDryFriction);
        mModelOecVertDriver  = QElectDriver(mModelOecVertElectMotor  , mpVertLoadModel, marrSpreadParams);


        // РЕАЛЬНЫЙ ПРИВОД
 // QString str =ui->tableWidget-> item(1,0)->text();
        mRealInductL =   (ui->tableWidget-> item(1,0)->text().toDouble())/ 1000.;

        mRealResist = ui->tableWidget-> item(1,1)->text().toDouble();

        mRealPsi_f = ui->tableWidget-> item(1,2)->text().toDouble();

        //мом инерции ротора
          mRealJ0= CONST_J_ROTOR;

        // момент сопротивления эл двигателя при обесточенных обмотках, не более
          mRealMaxAmpMomResidual =  ui->tableWidget->item(1,6)->text().toDouble();
          mRealPhMomResidual =  ui->tableWidget->item(1,7)->text().toDouble()/ 180. * M_PI;

        // mRealMomOut = ui->tableWidget-> item(1,3)->text().toDouble();

    // НАГРУЗКА
        //мом инерции нагрузки
          mRealJPayLoad= ui->tableWidget-> item(1,3)->text().toDouble();
        // коэффициент сопротивления нагрузки
          mRealCx= ui->tableWidget-> item(1,4)->text().toDouble();
           mRealCv= ui->tableWidget-> item(1,5)->text().toDouble();
           if(ui->comboBox_3->currentIndex()==0)
           {
               mHorLoadReal = QLoad ( mRealJPayLoad , mRealCx,mRealCv, VAL_TEST_dOm_po_dt);
               mpHorLoadReal =  &mHorLoadReal;

               mVertLoadReal = QLoad ( mRealJPayLoad , mRealCx,mRealCv, VAL_TEST_dOm_po_dt);
               mpVertLoadReal =  &mVertLoadReal;
           }

           if(ui->comboBox_3->currentIndex()==1)
           {
            QSegment SgmBarrier;
            TPlane Plane;
            mHorBrlsReal =  QBrls ( mRealCv, mRealCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                                ,VAL_BRLS_dOm_po_dt, SgmBarrier, Plane);
            mpHorLoadReal = &mHorBrlsReal;

            mVertBrlsReal =  QBrls ( mRealCv, mRealCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                                ,VAL_BRLS_dOm_po_dt, SgmBarrier, Plane);
            mpVertLoadReal = &mVertBrlsReal;
           }


           mOecHorElectMotor  = QElectMotor (mRealInductL , mRealResist ,
           mRealPsi_f , mRealJ0 , mRealMaxAmpMomResidual, mRealPhMomResidual, QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual)
             ,QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual) ,mMomDryFriction  );
            mOecHorDriver  = QElectDriver(mOecHorElectMotor  , mpHorLoadReal, marrSpreadParams);

            mOecVertElectMotor  = QElectMotor (mRealInductL , mRealResist ,
            mRealPsi_f , mRealJ0 , mRealMaxAmpMomResidual, mRealPhMomResidual, QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual)
              ,QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual) ,mMomDryFriction  );
             mOecVertDriver  = QElectDriver(mOecVertElectMotor  , mpVertLoadReal, marrSpreadParams);


          for (int i =0; i < 4; ++i)
          {
              mCmpArrLamb[i].m_Re = ui->tableWidget_2->item(0,i)->text().toDouble();
              mCmpArrLamb[i].m_Im = ui->tableWidget_2->item(1,i)->text().toDouble();
          }

    //  начальное угловая скорость
    mOmegaBegin = (ui->tableWidget_3->item(1,0)->text().toDouble())/ 180. * M_PI;

    //  начальное угловое пложение
    mTettaBegin = (ui->tableWidget_3->item(0,0)->text().toDouble())/ 180. * M_PI;

   marrObjective[0] = (ui->tableWidget_6->item(0,0)->text().toDouble())/ 180. * M_PI;
   marrObjective[1] = (ui->tableWidget_6->item(1,0)->text().toDouble())/ 180. * M_PI;
   mTargR = (ui->tableWidget_6->item(2,0)->text().toDouble());
   // без шумов
   if(ui->checkBox->isChecked()==true)
   {
      BEZ_SHUMOV = true;
   }
   else
   {
    BEZ_SHUMOV = false;
   }

   // интервал времени решения задачи управления

  VAlIntegrStep =constDbl_mh = 1./ ui->doubleSpinBox_27->value();

  //

   //const double Wind_V = ui->doubleSpinBox_30->value();
  // const double Wind_Alf = (ui->doubleSpinBox_29->value()) * M_PI/ 180.;
  // mEnvironment =TEnvironment ( Wind_V,  Wind_Alf, 0.);

   double valBearing = (ui->doubleSpinBox_2->value()) * M_PI/ 180.;
           double valTargCourse = (ui->doubleSpinBox_3->value()) * M_PI/ 180.;
           double valTargZenitAng  = (ui->doubleSpinBox_14->value()) * M_PI/ 180.;
           double valV = (ui->doubleSpinBox_4->value());
           double valR = (ui->doubleSpinBox->value());
           double valH = (ui->doubleSpinBox_13->value());
           //double valT = 0.;

   QInitTargData InitTargData (valBearing, valTargCourse, valTargZenitAng ,valV
       ,valR,valH, mT0);



   switch(ui->comboBox_4->currentIndex())
   {
   case 0:
       mPrimitiveCirleUniformTraj = QPrimitiveCirleUniformTraj(mT0, mT0, marrObjective,mTargR);
       mpTargTraj = &mPrimitiveCirleUniformTraj;
       ui->doubleSpinBox_15->setVisible(false);
       break;
   case 1:
       mUniformRectLinMotion = QUniformRectLinMotion(mT0, ui->doubleSpinBox_5->value(), InitTargData );
       mpTargTraj = &mUniformRectLinMotion;
       ui->doubleSpinBox_15->setVisible(false);
       break;
   case 2:
       mEquallyAcceleratedMotion = QEquallyAcceleratedMotion(mT0, ui->doubleSpinBox_5->value()
                     , InitTargData,ui->doubleSpinBox_15->value() );
       mpTargTraj = &mEquallyAcceleratedMotion;
       ui->doubleSpinBox_15->setVisible(true);
       break;
   default:
       break;
   }


   for (int i = 0; i < 6; ++i)
   {
       if (i <3)
       {
           marrOecPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble();
           marrOecPosParams_Zv[i] =  ui->tableWidget_4->item(2,i)->text().toDouble();
       }
       else
       {
           marrOecPosParams[i] =  ui->tableWidget_4->item(0,i)->text().toDouble() * M_PI/ 180.;
           marrOecPosParams_Zv[i] =  ui->tableWidget_4->item(2,i)->text().toDouble()* M_PI/ 180.;
       }
   }
   ///

   // создание корабля
       // атмосфера
       double valWind_V = ui->doubleSpinBox_30->value();

       double valWind_VertV = 0.;

       double valWind_Alf =  ui->doubleSpinBox_29->value();
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
        mTimeTempSins = 1./ ui->doubleSpinBox_12->value();
       mSins = QWarSins (mEnvironment,valMaxSig_Q, valMaxSig_Psi, valMaxSig_Tet
                            ,valMaxSig_dQdt, valMaxSig_dPsidt,valMaxSig_dTetdt
                            ,valMaxSig_H,valMaxSig_VH,valK1
                            ,valSigV, mT0,mTimeTempSins );
       ///
            // 3.4  параметры корабля
          // константы
       double valVesselWidth = 40; // ширина(м)
       double valVesselLength = 200.; // длина(м)

       // парамеитры движения
       double valQ0 = ui->doubleSpinBox_8->value() * M_PI / 180. ; // генеральный курс

       double valVVess = ui->doubleSpinBox_7->value();
       mvalMaxQ =    3./180.*M_PI; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
       double valT_Q = 18.; // период рыскания
       mvalMaxPsi =      3./180.*M_PI;// максимальный угол килевой качки(амплитуда)
       double valT_Psi =  12; // период килевой качки
       mvalMaxTet =      12./180.*M_PI; //максимальный угол боротовой качки(амплитуда)
       double valT_Tet =  6.; // период бортовой качки
       double valMaxVert =     1. ;


       double valMaxAmp_AftFlexure  = 1. * M_PI/180.;
       // период колебаний кормового изгиба
       double valT_AftFlexure = 4.;
       //максимальная амплитуда бортового изгиба корабля в рад на 100 м
       double valMaxAmp_BoardFlexure =  1. * M_PI/180.;
       // период колебаний бортового изгиба
       double valT_BoardFlexure = 2.;

       // создание платформы ОЭК

         // создание устройства
       const double SigV = ui->doubleSpinBox_31->value()/ 1000.;
       const double SigU = SigV;
       const double SigR = ui->doubleSpinBox_33->value();
       const double DiagWidth =  1. * M_PI/3000. ;
       const double TimeTemp = 1./ ui->doubleSpinBox_32->value();
       mOEChannel = QOEChannel(SigV, SigU, SigR
                               , DiagWidth, TimeTemp,mT0 -TimeTemp);
       mpDevice = &mOEChannel;

       mDriveMeasImitator = QMeasurmentImitator (10., 10., 0.002 * 0.001/3.);

       const double VAlHorMomOut = 0.;
       const double VAlVertMomOut = 0.;
       mMobilePlatfOec = QMobilePlatf(mpDevice,  marrOecPosParams
                        ,mOecHorDriver, VAlHorMomOut  , mOecVertDriver, VAlVertMomOut
                        ,mTettaBegin, mOmegaBegin
                        ,mTettaBegin, mOmegaBegin, mDriveMeasImitator, mT0);
       mpPltfOec = &mMobilePlatfOec;


 mWarShip =  QWarShip (mEnvironment,valVesselWidth,valVesselLength
             ,mvalMaxQ ,valT_Q
             ,mvalMaxPsi,valT_Psi ,mvalMaxTet
             ,valT_Tet,valMaxVert, valQ0,valVVess
             ,valMaxAmp_AftFlexure,valT_AftFlexure,valMaxAmp_BoardFlexure
             ,valT_BoardFlexure, mSins,  mpPltfOec, mT0);

  mMovingT = ui->doubleSpinBox_6->value();

  if(ui->checkBox_2->isChecked())
  {
      bIDEAL_FILT = true;
  }
  else
  {
      bIDEAL_FILT = false;
  }



}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\АМЕТИСТ_2019\\Электропривод\\ОТЧЕТ\\РАБ_МАТЕРИАЛЫ");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}
//--------------------------------------------
void MainWindow::on_pushButton_3_clicked()
{

/*
#include "CoordSystTrsf.h"
    // ОТЛАДКА
    double arrKGSK0[6] = {-1000., 5000., 1000.,100., -300., 200.};
    double arrEilerCntrKP0[3]= {-0.2, 0.2,-0.3};// когда = {-M_PI/4., M_PI/4.,0.}; работает!!!
    double arrVectOmegaPSK0[3]= {0.5,-1, 1.};
    double arrGdgPosParams[6] = {100., 300.,-500.,0.,0.,0.};
   // 1.
    double arrPSK0[6] = {0.}, arrPSK1[6] = {0.};
    QWarShip::RecalcVect_KGSK_INTO_GdgPSK(arrKGSK0,arrEilerCntrKP0,arrVectOmegaPSK0, arrGdgPosParams
                            ,arrPSK0,6 );


    //////
    /// \brief dt
    ///
    double arrTargV0[3], arrTargOm0[3];


    double dt = 0.001;
    double arrKGSK1[6] = {0.};
    memcpy(arrKGSK1, arrKGSK0, 6 * sizeof(double));
    arrKGSK1[0] = arrKGSK0[0] +  arrKGSK0[3] * dt;
    arrKGSK1[1] = arrKGSK0[1] +  arrKGSK0[4] * dt;
    arrKGSK1[2] = arrKGSK0[2] +  arrKGSK0[5] * dt;

    double arrEilerCntrKP1[3] = {0.};

    arrEilerCntrKP1[0] = arrEilerCntrKP0[0] + arrVectOmegaPSK0[2] * dt;//Q
    arrEilerCntrKP1[1] = arrEilerCntrKP0[1] + arrVectOmegaPSK0[0] * dt;//Psi
    arrEilerCntrKP1[2] = arrEilerCntrKP0[2] + arrVectOmegaPSK0[1] * dt;// Tet







    QWarShip::RecalcVect_KGSK_INTO_GdgPSK(arrKGSK1,arrEilerCntrKP1,arrVectOmegaPSK0, arrGdgPosParams
                            ,arrPSK1,6 );
    double arrVPSK[3] = {0.};
    for (int i =0; i < 3;++i)
    {
      arrVPSK[i] = (arrPSK1[i] - arrPSK0[i]) /dt;
    }

    double arrPSK2[6] = {0.};
    QWarShip::RecalcVect_KGSK_INTO_GdgPSK_differentiation( arrKGSK0  , arrEilerCntrKP0, arrVectOmegaPSK0, arrGdgPosParams
                                                ,arrPSK2,6 );


    int uu=0;
   return;

    */
  inputData();
    wchar_t wchOutPutFold0[400] = {0};

    wcscpy(wchOutPutFold0, mwchOutPutFold);

  // начальная коррел матрица фильтра ускорения
    double arrFiltK_Begin[9] = {0.};

   // дисперсия корости ускорения для фильтра углвого ускорнеия привода
    double valDispAngAccel = mModelMaxAmpMomResidual * mModelMaxAmpMomResidual/2.
            +mWarShip.calcMaxDispAngAcceleration(mvalMaxQ,mvalMaxPsi,mvalMaxTet) ;

//valDispAngAccel = 9.;
    /// выбор алгоритма управления
    switch(ui->comboBox->currentIndex())
    {
    case 0:  // адаптация по W
        arrFiltK_Begin[0] = mDriveMeasImitator.mSgmTetta * mDriveMeasImitator.mSgmTetta;
        arrFiltK_Begin[4] = 0.01* 0.01;
        arrFiltK_Begin[8] = 0.01* 0.01;
        mCtrlHorFollowAdapt1 = QCtrlFollowAdapt1(mModelOecHorDriver.mElectMotor,mModelOecHorDriver.mpLoad
                        ,marrSpreadParams, 0., mTettaBegin, mOmegaBegin, mT0, mT0,constDbl_mh
                        ,mCmpArrLamb, arrFiltK_Begin,valDispAngAccel);
        mpHorCtrl = &mCtrlHorFollowAdapt1;


        mCtrlVertFollowAdapt1 = QCtrlFollowAdapt1(mModelOecVertDriver.mElectMotor,mModelOecVertDriver.mpLoad
                        ,marrSpreadParams, 0., mTettaBegin, mOmegaBegin, mT0, mT0,constDbl_mh
                        ,mCmpArrLamb, arrFiltK_Begin,valDispAngAccel);
        mpVertCtrl = &mCtrlVertFollowAdapt1;
        break;
    case 1:  // без адаптации (простое слежение)
        mCtrlHorFollow = QCtrlFollow(mModelOecHorDriver.mElectMotor,mModelOecHorDriver.mpLoad
                                     ,marrSpreadParams, 0., mTettaBegin, mOmegaBegin, mT0, mT0,constDbl_mh
                                     ,mCmpArrLamb);
        mpHorCtrl = &mCtrlHorFollow;

        mCtrlVertFollow = QCtrlFollow(mModelOecVertDriver.mElectMotor,mModelOecVertDriver.mpLoad
                                     ,marrSpreadParams, 0., mTettaBegin, mOmegaBegin, mT0, mT0,constDbl_mh
                                     ,mCmpArrLamb);
        mpVertCtrl = &mCtrlVertFollow;
        break;

    default:
        break;
    }

    // вычисление начального курсового угла цели в
    // собственной системе координат платформы ОЭК
    double arrEilerCntrKP[3] = {0.}, arrVGdg[3] = {0.};
    mWarShip.mSins.getEstArrEilers(mT0, arrEilerCntrKP);
    QWarShip::recalcPositionFromKGSK_to_GdgSphericalSK(arrEilerCntrKP
             , mpTargTraj->marrVectSostGSK,marrOecPosParams_Zv,arrVGdg);
    ///
    mOecMobilePlatfCtrl = QMobilePlatfCtrl(marrOecPosParams_Zv
            , mpHorCtrl, mpVertCtrl, mT0, arrVGdg[0] );

    QPlatfCtrl *pOecMobilePlatfCtrl = &mOecMobilePlatfCtrl;
    mOperCtrl = QOperativeCtrl (mT0,mT0,pOecMobilePlatfCtrl,mpTargTraj->marrVectSostGSK);



   const int QUantRows = int(mMovingT/VAlIntegrStep ) + 1;

   double *arrBuff = new double[QUantColsReport0 * QUantRows];
   memset(arrBuff, 0, sizeof(double) * QUantColsReport0 *(QUantRows ));
   int quantDoneSteps = -1;

   // вычисление истинного вектора углового положения
   // и скорости цели в собственной системе координат ОЭК
   double arrTargTruObjHor[2] = {0.}, arrTargTruObjVert[2] = {0.};
   double arrMu[3] = {0.}, arrShipOm[3] = {0.};
   mWarShip.integrateEilersVect (arrMu);
   mWarShip.integrate_arr_dEilers_po_dt(arrShipOm);

   QWarShip::recalcVectKGSK_to_ObjVects(mpTargTraj->marrVectSostGSK
       ,arrMu, arrShipOm, mWarShip.mpPltfOec->marrPosParams
       , arrTargTruObjHor,arrTargTruObjVert);


   double valQEst = arrTargTruObjHor[0];
   QWarShip::adjustQourse( valQEst, arrTargTruObjHor[0]);

   QTargTrackingOperation TargTrackingOperation(mOperCtrl,mWarShip
      ,mEnvironment, mpTargTraj,arrTargTruObjHor, arrTargTruObjVert);

   //
   wchar_t wchInpDataFile[400] = {0};
       wcscpy(  wchInpDataFile,  mwchOutPutFold);

      wcscat(wchInpDataFile, L"\\report.txt");
      _wremove(wchInpDataFile);
       createInputDataReport(wchInpDataFile,true,TargTrackingOperation) ;
       ///


       TargTrackingOperation.move(mT0, mMovingT, VAlIntegrStep
                               ,arrBuff, &quantDoneSteps);



   // вывод графиков
   int iLenName = 30;
   wchar_t *pwcharrColNames = new wchar_t[QUantColsReport0 * iLenName];
   memset(pwcharrColNames, 0, QUantColsReport0 * iLenName*sizeof(wchar_t));

   wcscpy(pwcharrColNames, L"t");
   wcscpy(&pwcharrColNames[iLenName   ],  L"HorDrivTetta");
   wcscpy(&pwcharrColNames[iLenName *2],  L"HorDrivOm");
   wcscpy(&pwcharrColNames[iLenName *3],  L"HorMomOutReal");

   wcscpy(&pwcharrColNames[iLenName *4],  L"VertDrivTetta");
   wcscpy(&pwcharrColNames[iLenName *5],  L"VertDrivOm");
   wcscpy(&pwcharrColNames[iLenName *6],  L"VertMomOutReal");

   wcscpy(&pwcharrColNames[iLenName *7],  L"HorTargRealTetta");
   wcscpy(&pwcharrColNames[iLenName *8],  L"HorTargRealOm");

   wcscpy(&pwcharrColNames[iLenName *9],  L"VertTargRealTetta");
   wcscpy(&pwcharrColNames[iLenName *10],  L"VertTargRealOm");

   wcscpy(&pwcharrColNames[iLenName *11],  L"HorEstMomOut");
   wcscpy(&pwcharrColNames[iLenName *12],  L"VertEstMomOut");

   wcscpy(&pwcharrColNames[iLenName *13],  L"HorTettaEst");
   wcscpy(&pwcharrColNames[iLenName *14],  L"HorOmEst");

   wcscpy(&pwcharrColNames[iLenName *15],  L"VertTettaEst");
   wcscpy(&pwcharrColNames[iLenName *16],  L"VertOmEst");




   // оси  координат
   wchar_t wchAxesFileName0[300] ={0};
   wcscpy(  wchAxesFileName0,  wchOutPutFold0);
   wcscat(wchAxesFileName0, L"\\AxesArr.shp");
   TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
   ,-10000., 10000.,30.) ;
   // график
   double arrScale[100];
   for (int i = 0; i < 100;++i)
   {
    arrScale[i] = 1.;
   }

   arrScale[1] =arrScale[2] =arrScale[4] =arrScale[5] =arrScale[7]=arrScale[8] = 180./ M_PI;
   arrScale[9]=arrScale[10] = 180./ M_PI;
   arrScale[13]=arrScale[14] =arrScale[15]=arrScale[16] = 180./ M_PI;


 for (int i =1; i < QUantColsReport0; i++)
 {
   TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                   ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                   ,QUantColsReport0  // - к-во переменных о корорых накоплена информация в буфере
                                   ,quantDoneSteps //  - к-во точек
                                   ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                   ,iLenName // максимальная длина имени переменной
                                   ,0  //  номер переменной по оси X
                                   ,i  //  номер переменной по оси Y
                                   ,100.  //  масштаб по оси X
                                   ,arrScale[i]  // масштаб по оси Y
                                    ) ;
 }



 // построение графиков ошибок приводов по тетта и омега
    // формирование разностей
 for (int i = 0; i < quantDoneSteps; ++i)
 {
     double *p = &arrBuff[i * QUantColsReport0];
     p[1] -= p[7];
     p[2] -= p[8];
     p[3] = p[4] - p[9];
     p[4] = p[5] -p[10];
 }
 wcscpy(&pwcharrColNames[iLenName   ],  L"HorDelTetta");
 wcscpy(&pwcharrColNames[iLenName *2],  L"HorDelOm");
 wcscpy(&pwcharrColNames[iLenName *3],  L"VertDelTetta");
 wcscpy(&pwcharrColNames[iLenName *4],  L"VertDelOm");
 arrScale[1] = arrScale[2] = arrScale[3] = arrScale[4] = 1000.;
 for (int i =1; i < 5; i++)
 {
   TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                   ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                   ,QUantColsReport0  // - к-во переменных о корорых накоплена информация в буфере
                                   ,quantDoneSteps //  - к-во точек
                                   ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                   ,iLenName // максимальная длина имени переменной
                                   ,0  //  номер переменной по оси X
                                   ,i  //  номер переменной по оси Y
                                   ,100.  //  масштаб по оси X
                                   ,arrScale[i]  // масштаб по оси Y
                                    ) ;
 }
 ///




 // анализ разбросов траектории на последнем интервале времени на имитационной модели
  // [mMovingT - VAl_t; mMovingT]
 const double VAl_t = ((mMovingT - 3.) > 0.)?mMovingT - 3.:mMovingT-0.001;

 double arrHorSyst[2] = {0.},arrHorDisp[2] = {0.},arrVertSyst[2] = {0.},arrVertDisp[2] = {0.};

 QTargTrackingOperation::processingStatistic(arrBuff, quantDoneSteps, VAl_t
           , arrHorSyst, arrHorDisp,  arrVertSyst, arrVertDisp);

 ui->tableWidget_5->item(0,0)->setText(QString::number((arrHorSyst[0] * 1000.)));
 ui->tableWidget_5->item(0,1)->setText(QString::number((sqrt(arrHorDisp[0]) * 1000.)));
 ui->tableWidget_5->item(0,2)->setText(QString::number((sqrt(arrHorDisp[0] + arrHorSyst[0]* arrHorSyst[0]) * 1000.)));

 ui->tableWidget_5->item(1,0)->setText(QString::number((arrHorSyst[1] * 1000.)));
 ui->tableWidget_5->item(1,1)->setText(QString::number((sqrt(arrHorDisp[1]) * 1000.)));
 ui->tableWidget_5->item(1,2)->setText(QString::number((sqrt(arrHorDisp[1] + arrHorSyst[1]* arrHorSyst[1]) * 1000.)));

 ui->tableWidget_7->item(0,0)->setText(QString::number((arrVertSyst[0] * 1000.)));
 ui->tableWidget_7->item(0,1)->setText(QString::number((sqrt(arrVertDisp[0]) * 1000.)));
 ui->tableWidget_7->item(0,2)->setText(QString::number((sqrt(arrVertDisp[0] + arrVertSyst[0]* arrVertSyst[0]) * 1000.)));

 ui->tableWidget_7->item(1,0)->setText(QString::number((arrVertSyst[1] * 1000.)));
 ui->tableWidget_7->item(1,1)->setText(QString::number((sqrt(arrVertDisp[1]) * 1000.)));
 ui->tableWidget_7->item(1,2)->setText(QString::number((sqrt(arrVertDisp[1] + arrVertSyst[1]* arrVertSyst[1]) * 1000.)));


 createOutputDataReport(wchInpDataFile, false
                ,arrHorSyst,arrHorDisp,  arrVertSyst,arrVertDisp);

 delete [] arrBuff;

}








void MainWindow::on_comboBox_3_currentIndexChanged(int index)
{
    QLoad Load;
    QSegment SgmBarrier;
    TPlane Plane;
    switch(index)
    {
    case 0: // тестовая нагрузка
        ui->tableWidget->item(0,3)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(0,4)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(0,5)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(1,3)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(1,4)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(1,5)->setText(QString::number(VAL_TEST_Cv));
        break;

    case 1:

        mBrls  = QBrls (VAL_BRLS_Cv,VAL_BRLS_Cx,VAL_BRLS_MASS,VAL_BRLS_LENGTH,VAL_BRLS_HEIGHT,VAL_BRLS_dOm_po_dt
                        , SgmBarrier,Plane);
        ui->tableWidget->item(0,3)->setText(QString::number(mBrls.mJPayLoad));
        ui->tableWidget->item(0,4)->setText(QString::number(mBrls.mCx));
        ui->tableWidget->item(0,5)->setText(QString::number(mBrls.mCv));
        ui->tableWidget->item(1,3)->setText(QString::number(mBrls.mJPayLoad));
        ui->tableWidget->item(1,4)->setText(QString::number(mBrls.mCx));
        ui->tableWidget->item(1,5)->setText(QString::number(mBrls.mCv));

        break;

    default:
        break;

    };

}



void MainWindow::on_tableWidget_2_cellChanged(int row, int column)
{
     int jj =  ui->tableWidget_2->columnCount();
     int jv =  ui->tableWidget_2->rowCount();
    if (mbtableWidget_2Init)
    {


        if ((row == 1) && ((column == 0) || (column == 1)))
        {
         ui->tableWidget_2->item(row,column)->setText(QString::number(0.));
        }

        double arr[8] ={0.};
        for (int i =0; i < 2; ++i)
          for (int j =0; j < 4 ; ++j)
          {
             arr[i * 4 + j] =  ui->tableWidget_2->item(i,j)->text().toDouble();
          }

        if (fabs(arr[6]) <= 0.0000001)
        {
        arr[6] = arr[7] = 0.;
        }
        else
        {
          arr[3] = arr[2] ;
          arr[7] = -arr[6] ;
        }

        for (int i =0; i < 2; ++i)
          for (int j =0; j < 4 ; ++j)
          {
             ui->tableWidget_2->item(i,j)->setText(QString::number( arr [ i * 4 + j]));

          }


    }    

}

void MainWindow::read_tableWidget_2(double *arr)
{
  for (int i=0; i < 2; ++i)
      for (int j=0; j < 4; ++j)
      {
       arr [ i * 4 + j] =  ui->tableWidget_2->item(i,j)->text().toDouble();
      }

}

void MainWindow::init_tableWidget_2(const int ncols, double *arr)
{
  //ui->tableWidget_2->setColumnCount(ncols) ;
/*
 for (int i=0; i < ui->tableWidget_2->rowCount(); ++i)
      for (int j=0; j < ncols; ++j)
      {
           QString qstr1 = QString::number(arr[i *ncols + j]);
           ptwi[i * 4 + j] = new QTableWidgetItem(qstr1);
           ui->tableWidget_2->setItem(0,j,ptwi[i * 4 + j]);
      }*/

}

void MainWindow::on_comboBox_4_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:
        ui->frame->setVisible(false);
        ui->tableWidget_6->setVisible(true);
        break;
    case 2:
        ui->frame->setVisible(true);
        ui->tableWidget_6->setVisible(false);
        ui->doubleSpinBox_15->setVisible(true);
        break;
    case 1:
        ui->frame->setVisible(true);
        ui->tableWidget_6->setVisible(false);
        ui->doubleSpinBox_15->setVisible(false);
    default:
        ui->frame->setVisible(true);
        ui->tableWidget_6->setVisible(false);
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


void MainWindow::on_tableWidget_cellChanged(int row, int column)
{
    if (!mbtableWidget_Init)
    {
        return;
    }
    if((row == 0)&((column == 6)||(column == 7)))
    {
      ui->tableWidget->item(row,column)->setText(QString::number( 0.));
    }
}
//--------------------------------------------
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

if(bIDEAL_FILT )
{
  fprintf(fw,"  идеальная фильтрация параметров движения цели - ДА\n");
}
else
{
  fprintf(fw,"  идеальная фильтрация параметров движения цели - НЕТ\n");
}
if(BEZ_SHUMOV )
{
  fprintf(fw,"  без шумов - ДА\n");
}
else
{
  fprintf(fw,"  без шумов - НЕТ\n");
}

fprintf(fw,"  шаг интегрирования (С) = %8.7f\n",VAlIntegrStep);
fclose(fw);
TargTrackingOperation.createInputDataReport(FileName, false);

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

