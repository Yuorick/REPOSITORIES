#include <Qt>
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
#include "ParamsID.h"

#include <QCheckBox>






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

extern const bool BEZ_SHUMOV = true;
//шаг фильрации
extern const double constDbl_mh = 1./20000.;//1./5000.;//

const double VAlIntegrStep = 1. / 20000.;//1. / 5000.;//

extern const double CONST_GARMON_COEFF;

extern const double CONST_J_ROTOR;

extern const double COeff_L ;

class QRezPointTraj;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

     ui->setupUi(this);

     mbtableWidget_2Init = false;

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
       QString qstr1 = QString::number(-40.);
       ptwi1 = new QTableWidgetItem(qstr1);
       ui->tableWidget_2->setItem(0,0,ptwi1);

       QTableWidgetItem* ptwi2 = nullptr;
       QString qstr2 = QString::number(-41.);
       ptwi2 = new QTableWidgetItem(qstr2);
       ui->tableWidget_2->setItem(0,1,ptwi2);

       QTableWidgetItem* ptwi3 = nullptr;
       QString qstr3 = QString::number(-42.);
       ptwi3 = new QTableWidgetItem(qstr3);
       ui->tableWidget_2->setItem(0,2,ptwi3);

       QTableWidgetItem* ptwi4 = nullptr;
       QString qstr4 = QString::number(-43.);
       ptwi3 = new QTableWidgetItem(qstr3);
       ui->tableWidget_2->setItem(0,3,ptwi3);


       mbtableWidget_2Init = true;

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


          ui->tableWidget_4->setColumnCount(6);
         ui->tableWidget_4->setRowCount(3 );
         ui->tableWidget_4->horizontalHeader()->setVisible(true);
         //pHorHeader2->setVisible(true);
         ui->tableWidget_4->verticalHeader()->setVisible(true);

         QStringList lst5;
         lst5<< "Tetta"<<"Omega"<<"I_d"<<"I_qu"<<"U_d"<<"U_qu";
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







  ui->comboBox_3->setCurrentIndex(0);


  ui->tableWidget_6->setColumnCount(1);
  ui->tableWidget_6->setRowCount(2);


  QTableWidgetItem* ptwi20 = nullptr;
  QString qstr20 = QString::number(0.);
  ptwi20 = new QTableWidgetItem(qstr20);
  ui->tableWidget_6->setItem(0,0,ptwi20);

  QTableWidgetItem* ptwi21 = nullptr;
  QString qstr21 = QString::number(5.);
  ptwi21 = new QTableWidgetItem(qstr21);
  ui->tableWidget_6->setItem(1,0,ptwi21);

  ui->tableWidget_6->item(0,0)->setText(QString::number(0.));
  ui->tableWidget_6->item(1,0)->setText(QString::number(5.));

  mMeasurmentImitator = QMeasurmentImitator ();

  mFiltr = QFiltr();

  ui->comboBox_2->setCurrentIndex(0);
  setIdentedParams(ui->comboBox_2->currentIndex());
}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------







void MainWindow:: inputData()
{

    create_IarrTarg();

    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

    mQuiqMaxU0 =   ui->doubleSpinBox_23->value();

   mMomDryFriction = ui->doubleSpinBox_29->value();

  // МАТЕМАТИЧЕСКАЯ МОДЕЛЬ
    mInductL =   (ui->tableWidget->item(0,0)->text().toDouble())/ 1000.;

    mResist = ui->tableWidget->item(0,1)->text().toDouble();

    mPsi_f = ui->tableWidget->item(0,2)->text().toDouble();

    //мом инерции ротора
      mJ0= CONST_J_ROTOR;

    // момент сопротивления эл двигателя при обесточенных обмотках, не более
     mMaxAmpMomResidual =  ui->tableWidget->item(0,7)->text().toDouble();
     mPhMomResidual =  ui->tableWidget->item(0,8)->text().toDouble()/ 180. * M_PI;
     mMomOut = ui->tableWidget->item(0,3)->text().toDouble();

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
      mDriverModel  = QElectDriver(mElectMotorModel  , mpLoadModel, marrSpreadParams, mMomOut );


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
        mDriverReal  = QElectDriver(mElectMotorReal  , mpLoadReal, marrSpreadParams, mRealMomOut);

        //
    mMovingT = ui->doubleSpinBox_25->value();
    mTimeU   = ui->doubleSpinBox_22->value();
    mModU    = ui->doubleSpinBox_23->value();
    mQuantIsp = ui->doubleSpinBox_26->value();
    mTrgOm = ui->doubleSpinBox_28->value() * M_PI / 180.;
        // Начальное угловое положение ротора серии испытаний
    mTetRotor0 = ui->doubleSpinBox_27->value()* M_PI / 180.;
    ///

    for (int i =0; i < 4; ++i)
    {
        mCmpArrLamb[i].m_Re = ui->tableWidget_2->item(0,i)->text().toDouble();
        mCmpArrLamb[i].m_Im = ui->tableWidget_2->item(1,i)->text().toDouble();
    }
    mCmpArrLamb[0].m_Im = 0.;


 //
   mMeasurmentImitator = QMeasurmentImitator (10., 10., 0.002 * 0.001/3.);
   //  mMeasurmentImitator = QMeasurmentImitator (0.000001, 0.000001, 0.000001);  // отладка   !!!


    //  начальное угловая скорость
    mOmegaBegin = (ui->tableWidget_3->item(1,0)->text().toDouble())/ 180. * M_PI;

    //  начальное угловое пложение
    mTettaBegin = (ui->tableWidget_3->item(0,0)->text().toDouble())/ 180. * M_PI;

   marrObjective[0] = (ui->tableWidget_6->item(0,0)->text().toDouble())/ 180. * M_PI;
   marrObjective[1] = (ui->tableWidget_6->item(1,0)->text().toDouble())/ 180. * M_PI;

   mTypeOfID = ui->comboBox_2->currentIndex();
}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\REPOSITORIES\\aircraft-model\\OUT_DRIVER");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}









void MainWindow::on_pushButton_3_clicked()
{


  inputData();
    wchar_t wchOutPutFold0[400] = {0};

    wcscpy(wchOutPutFold0, mwchOutPutFold);

    mModelCtrlFollow = QCtrlFollow(mDriverModel.mElectMotor,mDriverModel.mpLoad
                    ,marrSpreadParams, mMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
    mpModelCtrl = &mModelCtrlFollow;

    mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                 ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
    mpRealCtrl = &mRealCtrlFollow;

    ///



   // mCtrlModel=  QCtrl (mDriverModel, mTettaBegin, mOmegaBegin,arrC0);
    QStatSolutionParams StatSolutionParams;
    mpModelCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams);
    mpModelCtrl->mStatSolutionParams = StatSolutionParams;
    mpRealCtrl->mStatSolutionParams = StatSolutionParams;






    TComp cmparrRealEigenValues[4];
    mpRealCtrl->CalcEigenValues(marrObjective[0], marrObjective[1],  StatSolutionParams.marrGears,cmparrRealEigenValues );



   const int QUantRows = int(mMovingT/VAlIntegrStep ) + 1;

   double *arrBuff = new double[QUantColsCSVReport * QUantRows];
   memset(arrBuff, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
   int quantDoneSteps = -1;

  // double valt0 = 0.;//, valt1 = 0.;
   int isignum0 = -2;

   mDriveMoveImit = QDriveMoveImit (mpModelCtrl,mDriverReal,mMeasurmentImitator );
   const double VAlIntegratorTime  = 1.;
   mDriveMoveImit.move(0.,mMovingT, VAlIntegrStep,
                     VAlIntegratorTime ,arrBuff, &quantDoneSteps);


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

   arrScale[1] =arrScale[4] =arrScale[8] =arrScale[11] = 180./ M_PI;

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
 mpModelCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams0);

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
 mpRealCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams1);
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
 // расчетные с заложенной в мат модель

 ui->tableWidget_4->item(0,0)->setText(QString::number((marrObjective[0] + marrObjective[1] *mMovingT)/ M_PI * 180.));
 ui->tableWidget_4->item(0,1)->setText(QString::number(marrObjective[1]/ M_PI * 180.));
 ui->tableWidget_4->item(0,2)->setText(QString::number(0.));
 ui->tableWidget_4->item(0,3)->setText(QString::number(valCurrentIq * 1000));
 ui->tableWidget_4->item(0,4)->setText(QString::number(valUd));
 ui->tableWidget_4->item(0,5)->setText(QString::number(valUq));

  // расчетные с реальной модели
 ui->tableWidget_4->item(2,0)->setText(QString::number((marrObjective[0] + marrObjective[1] *mMovingT)/ M_PI * 180.));
 ui->tableWidget_4->item(2,1)->setText(QString::number(marrObjective[1]/ M_PI * 180.));
 ui->tableWidget_4->item(2,2)->setText(QString::number(0.));
 ui->tableWidget_4->item(2,3)->setText(QString::number(valRealCurrentIq * 1000));
 ui->tableWidget_4->item(2,4)->setText(QString::number(valRealUd));
 ui->tableWidget_4->item(2,5)->setText(QString::number(valRealUq));

  // с системы диф уравнений
 double *pp = &(arrBuff[(quantDoneSteps -1) *QUantColsCSVReport]);

 ui->tableWidget_4->item(1,0)->setText(QString::number(pp[4]/ M_PI * 180.));
 ui->tableWidget_4->item(1,1)->setText(QString::number(pp[1] / M_PI * 180.));
 ui->tableWidget_4->item(1,2)->setText(QString::number(pp[2] * 1000));
 ui->tableWidget_4->item(1,3)->setText(QString::number(pp[3] * 1000));
 ui->tableWidget_4->item(1,4)->setText(QString::number(pp[5]));
 ui->tableWidget_4->item(1,5)->setText(QString::number(pp[6]));




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








void MainWindow::on_comboBox_2_currentIndexChanged(int index)
{
 setIdentedParams(index);
}

void MainWindow::setIdentedParams(int index)
{
    switch(index)
    {
    case 0:
        ui->checkBox_8->setChecked(false);
        ui->checkBox_5->setChecked(false);
        ui->checkBox_7->setChecked(false);

        ui->checkBox_8->setVisible(false);
        ui->checkBox_5->setVisible(false);
        ui->checkBox_7->setVisible(false);

        ui->checkBox_8->setCheckable(false);
        ui->checkBox_5->setCheckable(false);
        ui->checkBox_7->setCheckable(false);

        break;
    case 1:
        ui->checkBox_8->setChecked(true);
        ui->checkBox_5->setChecked(true);
        ui->checkBox_7->setChecked(true);

        ui->checkBox_8->setVisible(true);
        ui->checkBox_5->setVisible(true);
        ui->checkBox_7->setVisible(true);

        ui->checkBox_8->setCheckable(true);
        ui->checkBox_5->setCheckable(true);
        ui->checkBox_7->setCheckable(true);

        break;
    default:
        break;
    };
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


void MainWindow::on_tableWidget_2_cellChanged(int row, int column)
{

    if (!mbtableWidget_2Init)
    {
        return;
    }
     int jj =  ui->tableWidget_2->columnCount();
     int jv =  ui->tableWidget_2->rowCount();



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

void MainWindow::read_tableWidget_2(double *arr)
{
  for (int i=0; i < 2; ++i)
      for (int j=0; j < 4; ++j)
      {
       arr [ i * 4 + j] =  ui->tableWidget_2->item(i,j)->text().toDouble();
      }

}




void MainWindow::on_pushButton_clicked()
{
    inputData();


      if(0 == mTypeOfID)
      {

          // делаем 2 набора стационарных решений с корнями -40, -41, -42,-43 и -50, -51, -52, -53
        // стационарные решения ищем для слежения с шагом по скорости deltaOmega от 0 град/с до 720 град/сек
            // и шагом по углу deltaTetta от 0 град до 360 град
            double deltaTetta = 18. / 180. * M_PI;
            int quantTetta = 1;
            double valOmega0 = 5./ 180. * M_PI;
            double deltaOmega = 72. / 180. * M_PI;
            int quantOm =  6;



            TComp cmpArrLamb[8];
            cmpArrLamb[0] = TComp(-20., 0.);
            cmpArrLamb[1] = TComp(-21., 0.);
            cmpArrLamb[2] = TComp(-22., 0.);
            cmpArrLamb[3] = TComp(-23., 0.);

            cmpArrLamb[4] = TComp(-50., 0.);
            cmpArrLamb[5] = TComp(-51., 0.);
            cmpArrLamb[6] = TComp(-52., 0.);
            cmpArrLamb[7] = TComp(-53., 0.);

            const int QUantLambSets = 1;


          //  mDriverModel = QParamsID::static_identify_by_StatMeth(mDriverModel, mDriverReal, mMeasurmentImitator, marrSpreadParams
                //                                           , deltaTetta,  quantTetta, valOmega0,deltaOmega,  quantOm
                //                                           ,mMovingT,cmpArrLamb,  QUantLambSets,miarrTargNums
                  //                                         , mlenTargNums,VAlIntegrStep );

            QStatSolutionParams *parrStatSolutionParams = new QStatSolutionParams[QUantLambSets * quantTetta * quantOm];
            double *arrDynamicDelX = new double[QUantLambSets * quantTetta * quantOm * QVARS];

            QParamsID ParamsID( mDriverModel, mDriverReal
                                 , mMeasurmentImitator,marrSpreadParams, miarrTargNums, mlenTargNums);

            ParamsID.imitateStatExperimentsRez(deltaTetta, quantTetta, valOmega0,deltaOmega, quantOm
                                               , mMovingT,cmpArrLamb, QUantLambSets, parrStatSolutionParams
                                               , arrDynamicDelX, VAlIntegrStep );

            int quantExprs = QUantLambSets * quantTetta * quantOm;
            //ParamsID.calcDelAlf_StatMeth(quantExprs, parrStatSolutionParams, arrDynamicDelX,arrDElAlf);
            QElectMotor MotorRez;
            QLoad LoadRez;
            double valMomOut = -1.;


            double valNeviaz = ParamsID.identify_by_StatMeth(quantExprs, parrStatSolutionParams
                     , arrDynamicDelX, &MotorRez, &LoadRez, &valMomOut);

            delete []parrStatSolutionParams;
            delete []arrDynamicDelX;
            //  присвоения в таблицу
            QElectDriver DriverRez(MotorRez,  &LoadRez, mDriverModel.marrSpreadParams, valMomOut);
           // ui->tableWidget->item(0,0)->setText(QString::number(1./DriverRez.getParam(0)* 1000./ COeff_L));
            ui->tableWidget->item(0,0)->setText(QString::number(1./DriverRez.getParam(0)* 1000.));
            ui->tableWidget->item(0,1)->setText(QString::number(DriverRez.getParam(1)));
            ui->tableWidget->item(0,2)->setText(QString::number(DriverRez.getParam(2)));
           // ui->tableWidget->item(0,3)->setText(QString::number(0.015));
           // ui->tableWidget->item(0,4)->setText(QString::number(0.));
            ui->tableWidget->item(0,3)->setText(QString::number(DriverRez.getParam(6)));
            ui->tableWidget->item(0,4)->setText(QString::number(DriverRez.getJPayLoad()));
            ui->tableWidget->item(0,5)->setText(QString::number(DriverRez.getParam(4)));
            ui->tableWidget->item(0,6)->setText(QString::number(DriverRez.getParam(5)));
            ui->tableWidget->item(0,7)->setText(QString::number(DriverRez.getParam(7)));
            ui->tableWidget->item(0,8)->setText(QString::number(DriverRez.getParam(8)/ M_PI * 180.));

            ui->doubleSpinBox_3->setValue(valNeviaz);
            return;

      }


      if(1 == mTypeOfID)
      {
          wchar_t wchOutPutFold0[400] = {0};

          wcscpy(wchOutPutFold0, mwchOutPutFold);



          // для динамического метода
          const int QUantRows = int(mMovingT/VAlIntegrStep ) + 1;

          double *arrBuff = new double[QUantColsCSVReport * QUantRows];
          memset(arrBuff, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
          int quantDoneSteps = -1;

          QParamsID ParamsID( mDriverModel, mDriverReal
                               , mMeasurmentImitator,marrSpreadParams,miarrTargNums
                              , mlenTargNums);

          ParamsID.imitateDynamicExperimentsRez(mModU,  mTimeU,mMovingT,mTrgOm, mTetRotor0, mQuantIsp
                                                 , VAlIntegrStep ,arrBuff, &quantDoneSteps);
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

          arrScale[1] =arrScale[4] =arrScale[8] =arrScale[11] = 180./ M_PI;

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

        delete []pwcharrColNames;
return;
        // формирование массива QRezPointTraj результатов экспериментов

        QRezPointTraj *parrRezPointTraj = new QRezPointTraj[quantDoneSteps];
        for (int i =0; i < quantDoneSteps; ++i)
        {
          parrRezPointTraj[i] =  QRezPointTraj(arrBuff[i * QUantColsCSVReport]
                  , &(arrBuff[i * QUantColsCSVReport + 1]), &(arrBuff[i * QUantColsCSVReport + 5]) );
          int uuu=0;
        }
        delete []arrBuff;
        QElectMotor MotorRez;
        QLoad LoadRez;
        double valMomOut = -1.;

        double valNeviaz = ParamsID.identify_by_DynamicMeth_2(VAlIntegrStep,parrRezPointTraj, quantDoneSteps
                                                          , &MotorRez, &LoadRez, &valMomOut);
        delete []parrRezPointTraj ;


        //  присвоения в таблицу
        QElectDriver DriverRez(MotorRez,  &LoadRez, mDriverModel.marrSpreadParams, valMomOut);
       // ui->tableWidget->item(0,0)->setText(QString::number(1./DriverRez.getParam(0)* 1000./ COeff_L));
        ui->tableWidget->item(0,0)->setText(QString::number(1./DriverRez.getParam(0)* 1000.));
        ui->tableWidget->item(0,1)->setText(QString::number(DriverRez.getParam(1)));
        ui->tableWidget->item(0,2)->setText(QString::number(DriverRez.getParam(2)));
       // ui->tableWidget->item(0,3)->setText(QString::number(0.015));
       // ui->tableWidget->item(0,4)->setText(QString::number(0.));
        ui->tableWidget->item(0,3)->setText(QString::number(DriverRez.getParam(6)));
        ui->tableWidget->item(0,4)->setText(QString::number(DriverRez.getJPayLoad()));
        ui->tableWidget->item(0,5)->setText(QString::number(DriverRez.getParam(4)));
        ui->tableWidget->item(0,6)->setText(QString::number(DriverRez.getParam(5)));
        ui->tableWidget->item(0,7)->setText(QString::number(DriverRez.getParam(7)));
        ui->tableWidget->item(0,8)->setText(QString::number(DriverRez.getParam(8)/ M_PI * 180.));

        ui->doubleSpinBox_3->setValue(valNeviaz);
        return;
}

      ///


}

//------------------------------------------
void MainWindow:: create_IarrTarg()
{
   mlenTargNums = 0;
   for (int i =0; i < (QPARAMS +2); ++i)
   {
      miarrTargNums[i] = -1;
   }
   QCheckBox *pCheckBox[QPARAMS +2];
   pCheckBox[0] = ui->checkBox;
   pCheckBox[1] = ui->checkBox_2;
   pCheckBox[2] = ui->checkBox_10;
   pCheckBox[3] = ui->checkBox_8;
   pCheckBox[4] = ui->checkBox_9;
   pCheckBox[5] = ui->checkBox_4;
   pCheckBox[6] = ui->checkBox_3;
   pCheckBox[7] = ui->checkBox_7;
   pCheckBox[8] = ui->checkBox_5;

   for (int i =0; i < (QPARAMS +2); ++i)
   {
     if(pCheckBox[i]->isChecked())
    {
       miarrTargNums[mlenTargNums] = i;
       mlenTargNums++;

    }
   }



}

void MainWindow::on_checkBox_stateChanged(int arg1)
{
inputData();
}

void MainWindow::on_checkBox_2_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_10_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_3_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_8_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_9_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_4_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_7_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_checkBox_5_stateChanged(int arg1)
{
    inputData();
}

void MainWindow::on_comboBox_2_activated(const QString &arg1)
{

}
