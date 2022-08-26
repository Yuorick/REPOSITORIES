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
//#include "URPolygon.h"
//#include "UrPointXY.h"
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
#include "CtrlVelo.h"
#include "Ctrl.h"
#include "StatSolutionParams.h"
#include  "CtrlVeloQuiq.h"
#include "CtrlPos.h"
//#include "ParamsID.h"







#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsCSVReport;
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
//extern const double constDbl_mh = 1./20000.;
extern  double constDbl_mh = 1./20000.;

const double VAlIntegrStep = 1. / 20000.;

extern const double CONST_J_ROTOR;

extern const double COeff_L ;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

     ui->setupUi(this);

     mbtableWidget_2Init = false;
    // ui->comboBox->setCurrentIndex(0);
     // создание таблицы параметров
        ui->tableWidget->setColumnCount(9);
        ui->tableWidget->setRowCount(2 );


        QHeaderView *pHorHeader = ui->tableWidget->horizontalHeader();
        pHorHeader->setVisible(true);
        QStringList lst;
        lst<< "L, мГ"<<"R, Ом"<<"Psi_f"<<"MomOut"<<"J_load"<<"C_om"<<"Cv"<<"AmpMResid"<<"PhMResid";


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

        ui->tableWidget->item(0,3)->setText(QString::number(0.));
        ui->tableWidget->item(0,4)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(0,5)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(0,6)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(0,7)->setText(QString::number(0.9 * 0.8));
        ui->tableWidget->item(0,8)->setText(QString::number(0.));



        ui->tableWidget->item(1,0)->setText(QString::number(0.0006* 1000.));
        ui->tableWidget->item(1,1)->setText(QString::number(0.2));
        ui->tableWidget->item(1,2)->setText(QString::number(0.066));
        //ui->tableWidget->item(1,3)->setText(QString::number(0.015));

        ui->tableWidget->item(1,3)->setText(QString::number(0.));
        ui->tableWidget->item(1,4)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(1,5)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(1,6)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(1,7)->setText(QString::number(0.9 * 0.8));
        ui->tableWidget->item(1,8)->setText(QString::number(0.));
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
       QString qstr1 = QString::number(-50.);
       ptwi1 = new QTableWidgetItem(qstr1);
       ui->tableWidget_2->setItem(0,0,ptwi1);

       QTableWidgetItem* ptwi2 = nullptr;
       QString qstr2 = QString::number(-51.);
       ptwi2 = new QTableWidgetItem(qstr2);
       ui->tableWidget_2->setItem(0,1,ptwi2);

       QTableWidgetItem* ptwi3 = nullptr;
       QString qstr3 = QString::number(-52.);
       ptwi3 = new QTableWidgetItem(qstr3);
       ui->tableWidget_2->setItem(0,2,ptwi3);

       QTableWidgetItem* ptwi30 = nullptr;
       QString qstr30 = QString::number(-53.);
       ptwi30 = new QTableWidgetItem(qstr30);
       ui->tableWidget_2->setItem(0,3,ptwi30);


       mCtrlVelo = true;
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


          ui->tableWidget_4->setColumnCount(10);
         ui->tableWidget_4->setRowCount(3 );
         ui->tableWidget_4->horizontalHeader()->setVisible(true);
         //pHorHeader2->setVisible(true);
         ui->tableWidget_4->verticalHeader()->setVisible(true);

         QStringList lst5;
         lst5<< "Tetta"<<"Omega"<<"I_d"<<"I_qu"<<"U_d"<<"U_qu"<<"TetSyst, мрад"<<"Tet SKZ, мрад"<<"OmSyst, мрад/с"<<"OM SKZ, мрад/с";
         ui->tableWidget_4->setHorizontalHeaderLabels(lst5);
         QStringList lst6;
         lst6<< "Мат. модель"<<"Диф.ур-я"<<"Реал. модель";
         ui->tableWidget_4->setVerticalHeaderLabels(lst6);
         QString qstr00 = QString::number(0.);

         for (int i = 0; i < (ui->tableWidget_4->rowCount()); ++i)
             for (int j = 0; j < (ui->tableWidget_4->columnCount()); ++j)
             {
                QTableWidgetItem* ptwi = nullptr;
                ptwi = new QTableWidgetItem(qstr00);
                ui->tableWidget_4->setItem(i,j,ptwi);
             }

         //---------------------------------------------

             ui->tableWidget_5->setColumnCount(4);
             ui->tableWidget_5->setRowCount(2 );
             ui->tableWidget_5->horizontalHeader()->setVisible(true);
             //pHorHeader2->setVisible(true);
             ui->tableWidget_5->verticalHeader()->setVisible(true);

             QStringList lst7;
             lst7<< "Lamb1"<<"Lamb2"<<"Lamb3"<<"Lamb4";
             ui->tableWidget_5->setHorizontalHeaderLabels(lst7);

             QStringList lst8;
             lst8<< "Re"<<"Im";
             ui->tableWidget_5->setVerticalHeaderLabels(lst8);

             for (int i = 0; i < (ui->tableWidget_5->rowCount()); ++i)
                 for (int j = 0; j < (ui->tableWidget_5->columnCount()); ++j)
                 {
                    QTableWidgetItem* ptwi = nullptr;
                    ptwi = new QTableWidgetItem(qstr00);
                    ui->tableWidget_5->setItem(i,j,ptwi);
                 }

         // ui->tableWidget_5->setColumnCount(3);


   //ui->comboBox->setCurrentIndex(3);
   ui->doubleSpinBox_23->setVisible(false);
   ui->doubleSpinBox_24->setVisible(false);

  ui->comboBox_3->setCurrentIndex(0);


  ui->tableWidget_6->setColumnCount(1);
  ui->tableWidget_6->setRowCount(2);


  QTableWidgetItem* ptwi20 = nullptr;
  QString qstr20 = QString::number(0.);
  ptwi20 = new QTableWidgetItem(qstr20);
  ui->tableWidget_6->setItem(0,0,ptwi20);

  QTableWidgetItem* ptwi21 = nullptr;
  QString qstr21 = QString::number(180.);
  ptwi21 = new QTableWidgetItem(qstr21);
  ui->tableWidget_6->setItem(1,0,ptwi21);

  ui->tableWidget_6->item(0,0)->setText(QString::number(0.));
  ui->tableWidget_6->item(1,0)->setText(QString::number(90.));

  mMeasurmentImitator = QMeasurmentImitator ();

  mFiltr = QFiltr();


}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------


void MainWindow:: inputData()
{

    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

    mQuiqMaxU0 =   ui->doubleSpinBox_23->value();

    // момент сухого трения
    mMomDryFriction = ui->doubleSpinBox_28->value();
    // МАТЕМАТИЧЕСКАЯ МОДЕЛЬ
      mInductL =   (ui->tableWidget->item(0,0)->text().toDouble())/ 1000.;

      mResist = ui->tableWidget->item(0,1)->text().toDouble();

      mPsi_f = ui->tableWidget->item(0,2)->text().toDouble();

      //мом инерции ротора
        mJ0= CONST_J_ROTOR;

      // момент сопротивления эл двигателя при обесточенных обмотках, не более
       mMaxAmpMomResidual =  ui->tableWidget->item(0,7)->text().toDouble();
       mPhMomResidual =  ui->tableWidget->item(0,8)->text().toDouble()/ 180. * M_PI;
       mModelMomOut = ui->tableWidget->item(0,3)->text().toDouble();

  // НАГРУЗКА
      //мом инерции нагрузки
        mJPayLoad= ui->tableWidget->item(0,4)->text().toDouble();
      // коэффициент сопротивления нагрузки
        mCx= ui->tableWidget->item(0,5)->text().toDouble();
        mCv= ui->tableWidget->item(0,6)->text().toDouble();


        int iii = ui->comboBox_3->currentIndex();
        if(ui->comboBox_3->currentIndex()==0) // нагрузка тестовая
        {
         mLoadModel = QLoad ( mJPayLoad , mCx,mCv, VAL_TEST_dOm_po_dt);
         mpLoadModel = &mLoadModel;
        }

        if(ui->comboBox_3->currentIndex()==1)
        {
            QSegment SgmBarrier;
            TPlane Plane;
         mBrls =  QBrls ( mCv, mCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                             ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane);
         mpLoadModel = &mBrls;
        }
       //


        marrSpreadParams[0] = mInductL * 0.02;
        marrSpreadParams[1] = mResist * 0.005;
        marrSpreadParams[2] = mPsi_f * 0.03;
        marrSpreadParams[3] = mpLoadModel->mJPayLoad * 0.05;
        marrSpreadParams[4] = mpLoadModel->mCx *0.05;
        marrSpreadParams[5] = mpLoadModel->mCv *0.05;
        mElectMotorModel  = QElectMotor (mInductL , mResist, mPsi_f, mJ0, mMaxAmpMomResidual, mPhMomResidual
                       ,QElectMotor::calc_MaxDispMomRezid (mMaxAmpMomResidual),0.,mMomDryFriction);
        mDriverModel  = QElectDriver(mElectMotorModel  , mpLoadModel, marrSpreadParams );


        // РЕАЛЬНЫЙ ПРИВОД
  QString str =ui->tableWidget-> item(1,0)->text();
        mRealInductL =   (ui->tableWidget-> item(1,0)->text().toDouble())/ 1000.;

        mRealResist = ui->tableWidget-> item(1,1)->text().toDouble();

        mRealPsi_f = ui->tableWidget-> item(1,2)->text().toDouble();

        //мом инерции ротора
          mRealJ0= CONST_J_ROTOR;

        // момент сопротивления эл двигателя при обесточенных обмотках, не более
          mRealMaxAmpMomResidual =  ui->tableWidget->item(1,7)->text().toDouble();
          mRealPhMomResidual =  ui->tableWidget->item(1,8)->text().toDouble()/ 180. * M_PI;

         mRealMomOut = ui->tableWidget-> item(1,3)->text().toDouble();

    // НАГРУЗКА
        //мом инерции нагрузки
          mRealJPayLoad= ui->tableWidget-> item(1,4)->text().toDouble();
        // коэффициент сопротивления нагрузки
          mRealCx= ui->tableWidget-> item(1,5)->text().toDouble();
           mRealCv= ui->tableWidget-> item(1,6)->text().toDouble();
           if(ui->comboBox_3->currentIndex()==0)
           {
            mLoadReal = QLoad ( mRealJPayLoad , mRealCx,mRealCv, VAL_TEST_dOm_po_dt);
            mpLoadReal =  &mLoadReal;
           }

           if(ui->comboBox_3->currentIndex()==1)
           {
            QSegment SgmBarrier;
            TPlane Plane;
            mBrlsReal =  QBrls ( mRealCv, mRealCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                                ,VAL_BRLS_dOm_po_dt, SgmBarrier, Plane);
            mpLoadReal = &mBrlsReal;
           }


          mElectMotorReal  = QElectMotor (mRealInductL , mRealResist ,
         mRealPsi_f , mRealJ0 , mRealMaxAmpMomResidual, mRealPhMomResidual, QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual)
           ,QElectMotor:: calc_MaxDispMomRezid (mRealMaxAmpMomResidual) ,mMomDryFriction  );
          mDriverReal  = QElectDriver(mElectMotorReal  , mpLoadReal, marrSpreadParams);


    mMovingT = ui->doubleSpinBox_22->value();

     mTPeriodParamsRenew = ui->doubleSpinBox_25->value();

     mTimeIntegrator = ui->doubleSpinBox_26->value();

    switch(ui->comboBox->currentIndex())
    {
    case 0:
    case 1:
    mCmpArrLamb[0].m_Re = ui->tableWidget_2->item(0,0)->text().toDouble();
    mCmpArrLamb[0].m_Im = 0.;
    mCmpArrLamb[1].m_Re  = ui->tableWidget_2->item(0,1)->text().toDouble();
    mCmpArrLamb[1].m_Im = ui->tableWidget_2->item(1,1)->text().toDouble();
    mCmpArrLamb[2].m_Re = ui->tableWidget_2->item(0,2)->text().toDouble();
    mCmpArrLamb[2].m_Im = -mCmpArrLamb[1].m_Im;
        break;
    case 2:
    case 3:
    case 4:
    case 5:
        for (int i =0; i < 4; ++i)
        {
            mCmpArrLamb[i].m_Re = ui->tableWidget_2->item(0,i)->text().toDouble();
            mCmpArrLamb[i].m_Im = ui->tableWidget_2->item(1,i)->text().toDouble();
        }

        break;
    default:
        break;
    };
 //
   mMeasurmentImitator = QMeasurmentImitator (10., 10., 0.002 * 0.001/3.);
   //  mMeasurmentImitator = QMeasurmentImitator (0.000001, 0.000001, 0.000001);  // отладка   !!!


    //  начальное угловая скорость
    mOmegaBegin = (ui->tableWidget_3->item(1,0)->text().toDouble())/ 180. * M_PI;

    //  начальное угловое пложение
    mTettaBegin = (ui->tableWidget_3->item(0,0)->text().toDouble())/ 180. * M_PI;

   marrObjective[0] = (ui->tableWidget_6->item(0,0)->text().toDouble())/ 180. * M_PI;
   marrObjective[1] = (ui->tableWidget_6->item(1,0)->text().toDouble())/ 180. * M_PI;

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
  constDbl_mh = 1./ ui->doubleSpinBox_27->value();

  //

   const double Wind_V = ui->doubleSpinBox_30->value();
   const double Wind_Alf = (ui->doubleSpinBox_29->value()) * M_PI/ 180.;
   mEnvironment =TEnvironment ( Wind_V,  Wind_Alf, 0.);


    mPrimitiveCirleUniformTraj = QPrimitiveCirleUniformTraj(0.,0., marrObjective,1000.  );
    mpTargTraj = &mPrimitiveCirleUniformTraj;
}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\АМЕТИСТ_2019\\Электропривод\\ОТЧЕТ\\РАБ_МАТЕРИАЛЫ");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}



void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{

}



void MainWindow::on_doubleSpinBoxMass_valueChanged(const QString &arg1)
{

}

void MainWindow::on_lineEdit_cursorPositionChanged(int arg1, int arg2)
{

}

void MainWindow::on_pushButton_3_clicked()
{

    // отладка
   /* double parrInp0[6] = { 1., 2.
                         ,4.,4.
                         ,3.,3.
                         };
    double parrInp2[6] = { 5., 7.,8.
                           ,14.,13.,31.
                         };
    double parrInp1[4] = { 5., -2.
                         ,-4.,24.
                         };
    double parrOut[9] = {0.}, parrOut1[9] = {0.};

MtrxMultMatrx_MultMatrxTransp(parrInp0,parrInp1,parrInp2,3,2, parrOut);


double arrT[6] = {0.};
MtrxMultMatrx(parrInp0,3, 2, parrInp1,2, arrT) ;
MtrxMultMatrxTransp(arrT,3, 2, parrInp2,3, parrOut1);
  int iir = 1;

  return;
  */




    ui->doubleSpinBox_24->setValue(-1.);
  inputData();
    wchar_t wchOutPutFold0[400] = {0};

    wcscpy(wchOutPutFold0, mwchOutPutFold);

    ///
    switch(ui->comboBox->currentIndex())
    {
    case 0: // управление по скорости вращения без быстродействия
        mModelCtrlVelo = QCtrlVelo (mDriverModel.mElectMotor,mDriverModel.mpLoad, marrSpreadParams, mModelMomOut
                                    ,mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpModelCtrl = &mModelCtrlVelo;

        mRealCtrlVelo = QCtrlVelo (mDriverReal.mElectMotor,mDriverReal.mpLoad
                                   ,marrSpreadParams, mRealMomOut, mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlVelo;
        break;
    case 1:// управление по скорости вращения + быстродействие
        mModelCtrlVeloQuiq = QCtrlVeloQuiq(mDriverModel.mElectMotor,mDriverModel.mpLoad
                              ,marrSpreadParams, mModelMomOut,mTettaBegin, mOmegaBegin, 0., 0.
                                          ,marrObjective,mQuiqMaxU0,constDbl_mh,mCmpArrLamb);
        mpModelCtrl = &mModelCtrlVeloQuiq;

        mRealCtrlVeloQuiq = QCtrlVeloQuiq (mDriverReal.mElectMotor,mDriverReal.mpLoad
                      ,marrSpreadParams, mRealMomOut, mTettaBegin, mOmegaBegin, 0., 0.
                                           ,marrObjective,mQuiqMaxU0,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlVeloQuiq;
        break;

    case 2:
        mModelCtrlPos= QCtrlPos(mDriverModel.mElectMotor,mDriverModel.mpLoad
                     ,marrSpreadParams, mModelMomOut, mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpModelCtrl = &mModelCtrlPos;

        mRealCtrlPos= QCtrlPos(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlPos;
        break;
    case 3:
        mModelCtrlFollow = QCtrlFollow(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mModelMomOut, mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpModelCtrl = &mModelCtrlFollow;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlFollow;
        break;
    case 4:

       mModelCtrlDnmkCorrect_M_R = QCtrlDnmkCorrect_M_R(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mModelMomOut, mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb
                         , mTPeriodParamsRenew, mTPeriodParamsRenew, 0.);
        mpModelCtrl = &mModelCtrlDnmkCorrect_M_R;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlFollow;

        break;
    case 5:
         mModelCtrlDnmkCorrect_W = QCtrlDnmkCorrect_W(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mModelMomOut, mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb
                        , mTPeriodParamsRenew,mTPeriodParamsRenew, mTPeriodParamsRenew);

        mpModelCtrl = &mModelCtrlDnmkCorrect_W;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0., 0.,constDbl_mh,mCmpArrLamb);
        mpRealCtrl = &mRealCtrlFollow;

        break;

     case 6:

        mRealCtrlTrad = QCtrlTrad ( mDriverReal.mElectMotor,mDriverReal.mpLoad
                                , marrSpreadParams,  mRealMomOut,mTettaBegin, mOmegaBegin
                                , 0., 0.,constDbl_mh,mCmpArrLamb
                                , mQuiqMaxU0 , mMovingT, mMovingT);
        mpRealCtrl = &mRealCtrlTrad;

        mModelCtrlTrad = QCtrlTrad ( mDriverReal.mElectMotor,mDriverReal.mpLoad
                                , marrSpreadParams,  mRealMomOut,mTettaBegin, mOmegaBegin
                                , 0., 0.,constDbl_mh,mCmpArrLamb
                                , mQuiqMaxU0 , mMovingT, mMovingT);
        mpModelCtrl = &mModelCtrlTrad;
        break;
    default:
        break;
    };


    QStatSolutionParams StatSolutionParams;
    mpModelCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams,mModelMomOut);
    //mpModelCtrl->mStatSolutionParams = StatSolutionParams;
   // (*mpRealCtrl).mStatSolutionParams = StatSolutionParams;


    TComp cmparrRealEigenValues[4];

    mpRealCtrl->CalcEigenValues(marrObjective[0], marrObjective[1],  StatSolutionParams.marrGears,mRealMomOut,cmparrRealEigenValues );

    switch(ui->comboBox->currentIndex())
    {
    case 0:
    case 1:
        ui->tableWidget_5->item(0,0)->setText(QString::number(cmparrRealEigenValues[0].m_Re));
        ui->tableWidget_5->item(1,0)->setText(QString::number(cmparrRealEigenValues[0].m_Im));

        ui->tableWidget_5->item(0,1)->setText(QString::number(cmparrRealEigenValues[1].m_Re));
        ui->tableWidget_5->item(1,1)->setText(QString::number(cmparrRealEigenValues[1].m_Im));

        ui->tableWidget_5->item(0,2)->setText(QString::number(cmparrRealEigenValues[2].m_Re));
        ui->tableWidget_5->item(1,2)->setText(QString::number(cmparrRealEigenValues[2].m_Im));
         break;
    case 2:
    case 3:
        ui->tableWidget_5->item(0,0)->setText(QString::number(cmparrRealEigenValues[0].m_Re));
        ui->tableWidget_5->item(1,0)->setText(QString::number(cmparrRealEigenValues[0].m_Im));

        ui->tableWidget_5->item(0,1)->setText(QString::number(cmparrRealEigenValues[1].m_Re));
        ui->tableWidget_5->item(1,1)->setText(QString::number(cmparrRealEigenValues[1].m_Im));

        ui->tableWidget_5->item(0,2)->setText(QString::number(cmparrRealEigenValues[2].m_Re));
        ui->tableWidget_5->item(1,2)->setText(QString::number(cmparrRealEigenValues[2].m_Im));

        ui->tableWidget_5->item(0,3)->setText(QString::number(cmparrRealEigenValues[3].m_Re));
        ui->tableWidget_5->item(1,3)->setText(QString::number(cmparrRealEigenValues[3].m_Im));

        break;
    default:
        break;
    };



   const int QUantRows = int(mMovingT/VAlIntegrStep ) + 1;

   double *arrBuff = new double[QUantColsCSVReport * QUantRows];
   memset(arrBuff, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
   int quantDoneSteps = -1;

   double valt0 = 0., valt1 = 0.;
   int isignum0 = -2;   

   mDriveMoveImit = QDriveMoveImit (mpModelCtrl,mDriverReal,mMeasurmentImitator
                                    , mEnvironment,mpTargTraj,mRealMomOut, mModelMomOut);

   //const double VAlIntegratorTime = 2;

   mDriveMoveImit.move(marrObjective,0., mMovingT, VAlIntegrStep
                    ,arrBuff, &quantDoneSteps);

   ui->tableWidget->item(0,1)->setText(QString::number(mDriveMoveImit.mpCtrl->mElectMotor.mResist));
   ui->tableWidget->item(0,3)->setText(QString::number(mModelMomOut));


   // ОТЛАДКА ПОТОМ УБРАТЬ !!!
 //  const double valU1 = mpModelCtrl->mStatSolutionParams.marrStatU[0];
 //  const double valU2 = mpModelCtrl->mStatSolutionParams.marrStatU[1];
  // double arrStatPhVect[QVARS] = {0.};
  // mpModelCtrl->calcStableSolution( valU1,  valU2, arrStatPhVect);
  // int yyu = 0;
   /// КОНЕЦ ОТЛАДКИ

   // анализ разбросов траектории на последнем интервале времени на имитационной модели
    // [mMovingT - TPeriod; mMovingT]

   // вектор систематической ошибки
   double arrSyst[2] = {0.};
   // веутор диспесий M(x-xsyst)*(x-xsyst)
   double arrDisp[2] = {0.};
   // ПЕРИОД гармонической составляющей остат момента
   const double VAl_TPeriod = 2. * M_PI /(mDriverReal.mElectMotor.mZp * marrObjective[1])*2. ;
           //(mDriverReal.mElectMotor.mZp * marrObjective[1] / 2./M_PI) * 3.;
   mDriveMoveImit.StatisticProcess(arrBuff, quantDoneSteps,VAl_TPeriod, arrSyst, arrDisp);


   ui->tableWidget_4->item(1,6)->setText(QString::number((arrSyst[0] * 1000.)));
   ui->tableWidget_4->item(1,7)->setText(QString::number((sqrt(arrDisp[0]) * 1000.)));

   ui->tableWidget_4->item(1,8)->setText(QString::number((arrSyst[1] * 1000.)));
   ui->tableWidget_4->item(1,9)->setText(QString::number((sqrt(arrDisp[1]) * 1000.)));

   ///
  // анализ разбросов траектории в стационарном режиме аналитическим методом
   double arrW[2] = {0.};
   double arr_dx_po_dAlf[QVARS * QPARAMS] = {0.}, arrAHarmSquare[QVARS ] = {0.}, arrSystErr[QVARS ] = {0.}
         , arrAHarm[QVARS ] = {0.}, arrSumSig[QVARS ] = {0.}, arrRandomSig[QVARS ] = {0.};

   mDriveMoveImit.calcAnaliticSumTrajScatters(marrObjective,arrW, arr_dx_po_dAlf, arrAHarm, arrSystErr
                                              , arrRandomSig, arrSumSig);
   ui->tableWidget_4->item(0,6)->setText(QString::number((arrSystErr[3] * 1000.)));
   ui->tableWidget_4->item(0,7)->setText(QString::number(arrRandomSig[3]* 1000.));
   ui->tableWidget_4->item(0,8)->setText(QString::number((arrSystErr[0] * 1000.)));
   ui->tableWidget_4->item(0,9)->setText(QString::number(arrRandomSig[0] * 1000.));
   // вывод графиков
   int iLenName = 30;
   wchar_t *pwcharrColNames = new wchar_t[QUantColsCSVReport * iLenName];
   memset(pwcharrColNames, 0, QUantColsCSVReport * iLenName*sizeof(wchar_t));

   wcscpy(pwcharrColNames, L"t");
   wcscpy(&pwcharrColNames[iLenName   ],  L"Omega");
   wcscpy(&pwcharrColNames[iLenName *2],  L"I_d");
   wcscpy(&pwcharrColNames[iLenName *3],  L"I_qu");
   wcscpy(&pwcharrColNames[iLenName *4],  L"Tetta");
   wcscpy(&pwcharrColNames[iLenName *5],  L"U_d");
   wcscpy(&pwcharrColNames[iLenName *6],  L"U_qu");
   wcscpy(&pwcharrColNames[iLenName *7],  L"dOmega");
   wcscpy(&pwcharrColNames[iLenName *8],  L"Est_Omega");
   wcscpy(&pwcharrColNames[iLenName *9],  L"Est_I_d");
   wcscpy(&pwcharrColNames[iLenName *10],  L"Est_I_qu");
   wcscpy(&pwcharrColNames[iLenName *11],  L"Est_Tetta");

   wcscpy(&pwcharrColNames[iLenName *12],  L"TargTetta");
   wcscpy(&pwcharrColNames[iLenName *13],  L"TargOmega");



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

   arrScale[1] =arrScale[4] =arrScale[8] =arrScale[11] =arrScale[12] =arrScale[13] = 180./ M_PI;
   arrScale[2] =arrScale[3] =arrScale[9] =arrScale[10] = 1000.;

 for (int i =1; i < QUantColsCSVReport; i++)
 {
   TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                   ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                   ,QUantColsCSVReport  // - к-во переменных о корорых накоплена информация в буфере
                                   ,quantDoneSteps //  - к-во точек
                                   ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                   ,iLenName // максимальная длина имени переменной
                                   ,0  //  номер переменной по оси X
                                   ,i  //  номер переменной по оси Y
                                   ,100.  //  масштаб по оси X
                                   ,arrScale[i]  // масштаб по оси Y
                                    ) ;
 }


// построение графиков стационарных решений
 // создается массив arrBuff1 [26]
 //
 //туда собираются значения стационарных решений
 int quantcols = 13;
 double arrBuff1[13 *2] ={0.};
 arrBuff1[0] = 0.;
 arrBuff1[13] = mMovingT ;
 double valCurrentIq = 0., valUd = 0., valUq = 0., valDelFi = 0.;
 QStatSolutionParams StatSolutionParams0;
 mpModelCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams0,mModelMomOut);

 // переприсваивания , от "лени"
 double valOmegaStat = StatSolutionParams0.marrStatPhVect[0];
 valCurrentIq =  StatSolutionParams0.marrStatPhVect[2];
 valUd = StatSolutionParams0.marrStatU[0];
 valUq = StatSolutionParams0.marrStatU[1];

 arrBuff1[1] = arrBuff1[13 +1 ] = marrObjective[1];
 arrBuff1[3] = arrBuff1[13 +3 ] = valCurrentIq ;
 arrBuff1[5] = arrBuff1[13 +5 ] = valUd;
 arrBuff1[6] = arrBuff1[13 +6 ] = valUq;
 arrBuff1[4] = marrObjective[0]; // Tetta0
 arrBuff1[13 + 4] = marrObjective[0] + marrObjective[1] * mMovingT; // Tetta0 + Omega0 * t

  double valRealCurrentIq = 0., valRealUd = 0., valRealUq = 0., valRealDelFi = 0.;
  QStatSolutionParams StatSolutionParams1;
 mpRealCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams1,mRealMomOut);
 // переприсваивания , от "лени"

 valRealCurrentIq =  StatSolutionParams1.marrStatPhVect[2];
 valRealUd  = StatSolutionParams1.marrStatU[0];
 valRealUq = StatSolutionParams1.marrStatU[1];

 arrBuff1[6 + 1] = arrBuff1[13 + 6 +1 ] =marrObjective[1];
 arrBuff1[6 + 3] = arrBuff1[13 + 6 +3 ] =valRealCurrentIq ;
 arrBuff1[6 + 5] = arrBuff1[13 + 6 +5 ] =valRealUd;
 arrBuff1[6 + 6] = arrBuff1[13 + 6 +6 ] =valRealUq;


 wchar_t *pwcharrColNames1 = new wchar_t[quantcols * iLenName];
 memset(pwcharrColNames1, 0, quantcols * iLenName*sizeof(wchar_t));

 wcscpy(pwcharrColNames1, L"t");
 wcscpy(&pwcharrColNames1[iLenName   ],  L"ModelStat_Omega");
 wcscpy(&pwcharrColNames1[iLenName *2],  L"ModelStat_I_d");
 wcscpy(&pwcharrColNames1[iLenName *3],  L"ModelStat_I_qu");
 wcscpy(&pwcharrColNames1[iLenName *4],  L"ModelStat_Tetta");
 wcscpy(&pwcharrColNames1[iLenName *5],  L"ModelStat_U_d");
 wcscpy(&pwcharrColNames1[iLenName *6],  L"ModelStat_U_qu");

 wcscpy(&pwcharrColNames1[iLenName *7],  L"RealStat_Omega");
 wcscpy(&pwcharrColNames1[iLenName *8],  L"RealStat_I_d");
 wcscpy(&pwcharrColNames1[iLenName *9],  L"RealStat_I_qu");
 wcscpy(&pwcharrColNames1[iLenName *10],  L"RealStat_Tetta");
 wcscpy(&pwcharrColNames1[iLenName *11],  L"RealStat_U_d");
 wcscpy(&pwcharrColNames1[iLenName *12],  L"RealStat_U_qu");

 for(int j=0 ;j < quantcols; ++j)
 {
   arrScale[j] = 1.;
 }
 arrScale[1] =arrScale[4]=arrScale[7]=arrScale[10] = arrScale[4] = arrScale[13 + 4] =180./ M_PI;
 arrScale[2] = arrScale[3]= arrScale[8]= arrScale[9] = 1000.;
 for (int i =1; i < quantcols; i++)
 {
   TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                   ,arrBuff1 // массив с информацией - матрица nBuffRows x nBuffCols
                                   ,quantcols  // - к-во переменных о корорых накоплена информация в буфере
                                   ,2 //  - к-во точек
                                   ,pwcharrColNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
                                   ,iLenName // максимальная длина имени переменной
                                   ,0  //  номер переменной по оси X
                                   ,i  //  номер переменной по оси Y
                                   ,100.  //  масштаб по оси X
                                   ,arrScale[i] // масштаб по оси Y
                                    ) ;
 }

 // вывод в таблицу
 double *pp = &(arrBuff[(quantDoneSteps -1) *QUantColsCSVReport]);//pp[0] -  время
 // расчетные с заложенной в мат модель

 ui->tableWidget_4->item(0,0)->setText(QString::number((marrObjective[0] + marrObjective[1] *pp[0])/ M_PI * 180.));
 ui->tableWidget_4->item(0,1)->setText(QString::number(marrObjective[1]/ M_PI * 180.));
 ui->tableWidget_4->item(0,2)->setText(QString::number(0.));
 ui->tableWidget_4->item(0,3)->setText(QString::number(valCurrentIq * 1000));
 ui->tableWidget_4->item(0,4)->setText(QString::number(valUd));
 ui->tableWidget_4->item(0,5)->setText(QString::number(valUq));

  // расчетные с реальной модели
 ui->tableWidget_4->item(2,0)->setText(QString::number((marrObjective[0] + marrObjective[1] *pp[0])/ M_PI * 180.));
 ui->tableWidget_4->item(2,1)->setText(QString::number(marrObjective[1]/ M_PI * 180.));
 ui->tableWidget_4->item(2,2)->setText(QString::number(0.));
 ui->tableWidget_4->item(2,3)->setText(QString::number(valRealCurrentIq * 1000));
 ui->tableWidget_4->item(2,4)->setText(QString::number(valRealUd));
 ui->tableWidget_4->item(2,5)->setText(QString::number(valRealUq));

  // с системы диф уравнений


 ui->tableWidget_4->item(1,0)->setText(QString::number(pp[4]/ M_PI * 180.));
 ui->tableWidget_4->item(1,1)->setText(QString::number(pp[1] / M_PI * 180.));
 ui->tableWidget_4->item(1,2)->setText(QString::number(pp[2] * 1000));
 ui->tableWidget_4->item(1,3)->setText(QString::number(pp[3] * 1000));
 ui->tableWidget_4->item(1,4)->setText(QString::number(pp[5]));
 ui->tableWidget_4->item(1,5)->setText(QString::number(pp[6]));

 //double del = fmod(marrObjective[0] + marrObjective[1] *pp[0] -pp[4], 2. * M_PI);
 // построение графика ошибки по Omega
 double *parrx = new double[quantDoneSteps];
 double *parry = new double[quantDoneSteps];


 for (int i = 0; i < quantDoneSteps; ++i)
 {
     double *p = &(arrBuff[i * QUantColsCSVReport ]);
     parrx[i]= p[0] * 100.;
     parry[i] = (p[1]  - p[13] /*marrObjective[1]*/) * 1000000.;

 }
 TURPolyLine plnDelOmega(parrx,parry,quantDoneSteps);
 wchar_t wchplnDelOmega[400] ={0};
 wcscpy(  wchplnDelOmega,  wchOutPutFold0);
 wcscat(wchplnDelOmega, L"\\DelOmega.shp");
 plnDelOmega.WriteSetSHPFiles(wchplnDelOmega,&plnDelOmega, 1) ;
 ///
// построение графика ошибки по Tetta

 double valTStart = arrBuff[0 ];
 for (int i = 0; i < quantDoneSteps; ++i)
 {
     double *p = &(arrBuff[i * QUantColsCSVReport ]);
     parrx[i]= p[0] * 100.;
     parry[i] = (p[4]  - p[12]/*(marrObjective[0] + marrObjective[1] *(p[0] - valTStart))*/) * 1000000.;

 }
 TURPolyLine plnDelTetta(parrx,parry,quantDoneSteps);
 wchar_t wchplnDelTetta[400] ={0};
 wcscpy(  wchplnDelTetta,  wchOutPutFold0);
 wcscat(wchplnDelTetta, L"\\DelTetta.shp");
 plnDelTetta.WriteSetSHPFiles(wchplnDelTetta,&plnDelTetta, 1) ;
 delete []parrx;
 delete []parry;
 ///


 ui->doubleSpinBox_24->setValue(valt1 * 1000.);

 delete[] pwcharrColNames;
 delete []pwcharrColNames1;
 delete [] arrBuff;

 // построение графика прямой - ограничения по угловой скорости
 const TURPointXY  pnt1(0., mLoadModel.mMax_dOm_po_dt)
         , pnt2(mMovingT * 100., mLoadModel.mMax_dOm_po_dt);
 TURPolyLine pln_MaxdOm(  pnt1,   pnt2) ;
 wcscpy(  wchAxesFileName0,  wchOutPutFold0);
 wcscat(wchAxesFileName0, L"\\Max_dOm_po_dt.shp");
 pln_MaxdOm.WriteSetSHPFiles(wchAxesFileName0,&pln_MaxdOm, 1) ;


}

void MainWindow::on_pushButton_10_clicked()
{
    /*
    mqstrCandidatesSCVFIle = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с комплектами перед. чисел"
                         , "D:\\REPOSITORIES\\aircraft-model\\OUT_DRIVER");
   //////this->ui->lineEdit_2->setText(mqstrCandidatesSCVFIle);

     wchar_t array [400] = {0};
     mqstrCandidatesSCVFIle.toWCharArray(array);
     array[mqstrCandidatesSCVFIle.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
   int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) -2;

   // ui->doubleSpinBox_18  ->setValue((double)numCandidates );




 //   TYrRead::ReadHdrFromGearFileForDriver(mpwchPickedOutCandidates
 //       , &mbVelo, &mInductL, &mResist,  &mPsi_f, &mJPayLoad
  //                            , &mCx_om , &mOmegaStat);
    if (mbVelo)
    {
   // ui->lineEdit_3->setText("VELOCITY");
    }
    else
    {
      //  ui->lineEdit_3->setText("POSITION");
    }
    ui->doubleSpinBox_5->setValue(mInductL * 1000.);
    ui->doubleSpinBox_11->setValue(mResist);
    ui->doubleSpinBox_8->setValue(mPsi_f);
    ui->doubleSpinBox_10->setValue(mJPayLoad);
    ui->doubleSpinBox_7->setValue(mCx_om);
    ui->doubleSpinBox_25->setValue(mOmegaStat);

    ui->doubleSpinBox_15->setValue(mInductL * 1000.);
    ui->doubleSpinBox_21->setValue(mResist);
    ui->doubleSpinBox_19->setValue(mPsi_f);
    ui->doubleSpinBox_24->setValue(mJPayLoad);
    ui->doubleSpinBox_23->setValue(mCx_om);
    ui->doubleSpinBox_26->setValue(mOmegaStat);

    mRealInductL = mInductL;
    mRealResist = mResist;
    mRealPsi_f = mPsi_f;
    mRealJPayLoad = mJPayLoad;
    mRealCx_om = mCx_om;
    mRealOmegaStat = mOmegaStat;
    */
}




void MainWindow::on_pushButton_3_pressed()
{

}

void MainWindow::on_comboBox_currentIndexChanged(int index)
{


    ui->doubleSpinBox_23->setVisible(false);
    ui->doubleSpinBox_24->setVisible(false);
    switch(index)
    {
    case 6:
        ui->doubleSpinBox_23->setVisible(true);
        break;

    case 1:
        ui->doubleSpinBox_23->setVisible(true);
        ui->doubleSpinBox_24->setVisible(true);
        break;
    default:
        break;
    };

    int jjj = ui->tableWidget_2->rowCount();
    int jdd = ui->tableWidget_2->columnCount();
     // ui->tableWidget_2->setRowCount(2);
    ///
    QStringList lst1, lst2,lst3,lst4;
    if( index < 2)
    {
        if (!mCtrlVelo)
        {
            ui->tableWidget_2->item(1,0)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,1)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,2)->setText(QString::number(0.));

            ui->tableWidget_2->item(0,0)->setText(QString::number(-8.));
            ui->tableWidget_2->item(0,1)->setText(QString::number(-9.));
            ui->tableWidget_2->item(0,2)->setText(QString::number(-10.));
            ui->tableWidget_2->item(0,3)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,3)->setText(QString::number(0.));
            mbtableWidget_2Init = true;
            mCtrlVelo = true;

            ui->tableWidget_6->item(0,0)->setText(QString::number(0.));
            ui->tableWidget_6->item(1,0)->setText(QString::number(5.));
        }

      ui->tableWidget_5->setColumnCount(3);
      lst4<<"Lamb1"<<"Lamb2"<<"Lamb3";
    }


    if( index >= 2)
    {
        if (mCtrlVelo)
        {
            ui->tableWidget_2->item(1,0)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,1)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,2)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,3)->setText(QString::number(0.));

            ui->tableWidget_2->item(0,0)->setText(QString::number(-30.));
            ui->tableWidget_2->item(0,1)->setText(QString::number(-31.));
            ui->tableWidget_2->item(0,2)->setText(QString::number(-33.));
            ui->tableWidget_2->item(0,3)->setText(QString::number(-34.));

             mbtableWidget_2Init = true;
             mCtrlVelo = false;
             ui->tableWidget_6->item(0,0)->setText(QString::number(180.));
             ui->tableWidget_6->item(1,0)->setText(QString::number(0.));
             if (index == 3)
             {
               ui->tableWidget_6->item(1,0)->setText(QString::number(10.));
             }
        }

        ui->tableWidget_5->setColumnCount(4);
        lst4<<"Lamb1"<<"Lamb2"<<"Lamb3"<<"Lamb4";

    }


   ui->tableWidget_2->setHorizontalHeaderLabels(lst1);
   ui->tableWidget_5->setHorizontalHeaderLabels(lst4);
   for (int i=0; i<  ui->tableWidget_5->rowCount() ; i++)
                for (int j = 0; j < ui->tableWidget_5->columnCount(); j++)
                {
                    QTableWidgetItem* ptwi0 = nullptr;
                    ptwi0 = new QTableWidgetItem(QString::number(0.));
                    ui->tableWidget_5->setItem(i,j,ptwi0);
                }
///



}

void MainWindow::on_comboBox_activated(const QString &arg1)
{

}

void MainWindow::on_comboBox_editTextChanged(const QString &arg1)
{

}



void MainWindow::on_comboBox_3_currentIndexChanged(int index)
{
    QLoad Load;
    QSegment SgmBarrier;
    TPlane Plane;
    switch(index)
    {
    case 0: // тестовая нагрузка
        ui->tableWidget->item(0,4)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(0,5)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(0,6)->setText(QString::number(VAL_TEST_Cv));
        ui->tableWidget->item(1,4)->setText(QString::number(VAL_TEST_J));
        ui->tableWidget->item(1,5)->setText(QString::number(VAL_TEST_Cx));
        ui->tableWidget->item(1,6)->setText(QString::number(VAL_TEST_Cv));
        break;

    case 1:

        mBrls  = QBrls (VAL_BRLS_Cv,VAL_BRLS_Cx,VAL_BRLS_MASS,VAL_BRLS_LENGTH,VAL_BRLS_HEIGHT,VAL_BRLS_dOm_po_dt
                        , SgmBarrier,Plane);
        ui->tableWidget->item(0,4)->setText(QString::number(mBrls.mJPayLoad));
        ui->tableWidget->item(0,5)->setText(QString::number(mBrls.mCx));
        ui->tableWidget->item(0,6)->setText(QString::number(mBrls.mCv));
        ui->tableWidget->item(1,4)->setText(QString::number(mBrls.mJPayLoad));
        ui->tableWidget->item(1,5)->setText(QString::number(mBrls.mCx));
        ui->tableWidget->item(1,6)->setText(QString::number(mBrls.mCv));

        break;

    default:
        break;

    };

}

void MainWindow::on_comboBox_3_currentIndexChanged(const QString &arg1)
{
}

void MainWindow::on_tableWidget_2_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{

}

void MainWindow::on_tableWidget_2_cellChanged(int row, int column)
{
     int jj =  ui->tableWidget_2->columnCount();
     int jv =  ui->tableWidget_2->rowCount();
    if (mbtableWidget_2Init)
    {
        if (ui->comboBox->currentIndex() <= 1)
        {
            if ((row == 1) && (column == 0))
            {
             ui->tableWidget_2->item(row,column)->setText(QString::number(0.));
            }

            double arr[6] ={0.};
            for (int i =0; i < 2; ++i)
              for (int j =0; j < 3 ; ++j)
              {
                 arr[i * 3 + j] =  ui->tableWidget_2->item(i,j)->text().toDouble();
              }

            if (fabs(arr[4]) <= 0.0000001)
            {
            arr[4] = arr[5] = 0.;
            }
            else
            {
              arr[1] = arr[2] ;
              arr[5] = -arr[4] ;
            }

            for (int i =0; i < 2; ++i)
              for (int j =0; j < 3 ; ++j)
              {
                 ui->tableWidget_2->item(i,j)->setText(QString::number( arr [ i * 3 + j]));

              }
            ui->tableWidget_2->item(0,3)->setText(QString::number(0.));
            ui->tableWidget_2->item(1,3)->setText(QString::number(0.));

        }
        ///

        if (ui->comboBox->currentIndex() > 1)
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



void MainWindow::on_comboBox_currentIndexChanged(const QString &arg1)
{

}
