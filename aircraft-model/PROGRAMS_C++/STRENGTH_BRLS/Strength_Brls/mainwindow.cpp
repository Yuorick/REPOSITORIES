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

extern const bool BEZ_SHUMOV = false;
//шаг фильрации
extern const double constDbl_mh = 1./20000.;

const double VAlIntegrStep =1. / 20000.;

extern const double CONST_J_ROTOR;

extern const double COeff_L ;




/// тело БРЛС
extern const int NumPartFurkCaske_2540 ;
extern const int NumPointsFurk_Caske  ;

extern const double arrPointsCaske_Furk_2540[];


extern const int PartsCaske_Furk_2540[] ;
///

// толстый текстолит
extern const int NumPartTextFurk_2540;
extern const int NumPointsTextFurk_2540 ;

extern const double arrPoints_TextFurk_2540[] ;
extern const int Parts_TextFurk_2540[] ;
///

// алюминиевый тавр на антенне
extern const int NumPartAlFurk_2540 ;
extern const int NumPointsAlFurk_2540  ;

extern const double arrPoints_AlFurk_2540[];
extern const int Parts_AlFurk_2540[] ;
//---------------------------------------------------------------------------
// алюминиевый Пластина на приводе
extern const int NumPartDrvFurk_2540 ;
extern const int NumPointsDrvFurk_2540;

extern const double arrPoints_DrvFurk_2540[];
extern const int Parts_DrvFurk_2540[] ;

extern const double ARrBeginPlot_Furk2540[];
extern const double QUantPlots_Furk2540 ;
///
extern const int QUantColsReport_Wind  ;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

     ui->setupUi(this);

   /*  mbtableWidget_2Init = false;
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
       QString qstr1 = QString::number(-42.);
       ptwi1 = new QTableWidgetItem(qstr1);
       ui->tableWidget_2->setItem(0,0,ptwi1);

       QTableWidgetItem* ptwi2 = nullptr;
       QString qstr2 = QString::number(-40.);
       ptwi2 = new QTableWidgetItem(qstr2);
       ui->tableWidget_2->setItem(0,1,ptwi2);

       QTableWidgetItem* ptwi3 = nullptr;
       QString qstr3 = QString::number(-41.);
       ptwi3 = new QTableWidgetItem(qstr3);
       ui->tableWidget_2->setItem(0,2,ptwi3);

       QTableWidgetItem* ptwi30 = nullptr;
       QString qstr30 = QString::number(-43.);
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
  QString qstr21 = QString::number(5.);
  ptwi21 = new QTableWidgetItem(qstr21);
  ui->tableWidget_6->setItem(1,0,ptwi21);

  ui->tableWidget_6->item(0,0)->setText(QString::number(0.));
  ui->tableWidget_6->item(1,0)->setText(QString::number(5.));

  mMeasurmentImitator = QMeasurmentImitator ();

  mFiltr = QFiltr();
*/

}

MainWindow::~MainWindow()
{

    delete ui;
}


//-----------------------------------------------------------------------------------


void MainWindow:: inputData()
{
    /*TURPolygon *pplgarr = new TURPolygon[100];
     pplgarr[0] = TURPolygon( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                    ,arrPointsCaske_Furk_2540);


     QBeamPlot BeamPlot(0);

     QBeamPlot marrPlots[4];
    // marrPlots[0] = BeamPlot;
     const double arrRo [1] = {0.};
     const double arrSigma[1]= {0.};
     const double arrE [1]= {0.};
     marrPlots[0].marrPlgProfile[0] = pplgarr[0];

     marrPlots[0].marrPlgProfile[0] = TURPolygon ( NumPartFurkCaske_2540,NumPointsFurk_Caske,PartsCaske_Furk_2540
                                       ,arrPointsCaske_Furk_2540);


     // толстый текстолит
     marrPlots[0].marrPlgProfile[1] = TURPolygon ( NumPartTextFurk_2540,NumPointsTextFurk_2540,Parts_TextFurk_2540
                                      ,arrPoints_TextFurk_2540);

     // алюминиевый тавр на антенне
     marrPlots[0].marrPlgProfile[2] = TURPolygon ( NumPartAlFurk_2540,NumPointsAlFurk_2540,Parts_AlFurk_2540
                                      ,arrPoints_AlFurk_2540);

     // алюминиевая пластина на приводе
     marrPlots[0].marrPlgProfile[3] = TURPolygon ( NumPartDrvFurk_2540,NumPointsDrvFurk_2540,Parts_DrvFurk_2540
                                      ,arrPoints_DrvFurk_2540);



   // marrPlots[0] = QBeamPlot(  1, pplgarr, arrRo
                     //     ,arrSigma, arrE);
    QBeamPlot BeamPlot0;
    QBeamPlot::createBeamPlotFurke2540(0, BeamPlot0);
    /* marrPlots[1] = QBeamPlot(1);
     marrPlots[2] = QBeamPlot(2);
     marrPlots[3] = QBeamPlot(3);*/
   // delete []pplgarr;
  //  return;*/

    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

    // угловая скорость
    mOmega = (ui->doubleSpinBox_2->value()) / 180. * M_PI;

    // скорость ветра
     mWindV = ui->doubleSpinBox->value();

     mAccellerOmega = (ui->doubleSpinBox_3->value()) / 180. * M_PI;

     mPIce = ui->doubleSpinBox_11->value()*9.81;

     mHailDiam = ui->doubleSpinBox_20->value()* 0.001;

     mHailZ = ui->doubleSpinBox_21->value()* 0.01;

     mcs = 0.0005;


    int iii=0;

}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\АМЕТИСТ_2019\\БРЛС_ПРОЧНОСТЬ");
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
  inputData();
    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);

    /*
     * QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;
    */

    QSegment SgmBarrier;
    TPlane Plane;

    QBrlsBody BrlsBody( mCv, mCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                            ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane, mcs,0);


    BrlsBody.createBodyGraph(wchOutPutFold0);
    BrlsBody.createSideView(1.,1.,wchOutPutFold0);

/*
    ///
    switch(ui->comboBox->currentIndex())
    {
    case 0: // управление по скорости вращения без быстродействия
        mModelCtrlVelo = QCtrlVelo (mDriverModel.mElectMotor,mDriverModel.mpLoad, marrSpreadParams, mMomOut
                                    ,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpModelCtrl = &mModelCtrlVelo;

        mRealCtrlVelo = QCtrlVelo (mDriverReal.mElectMotor,mDriverReal.mpLoad
                                   ,marrSpreadParams, mRealMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpRealCtrl = &mRealCtrlVelo;
        break;
    case 1:// управление по скорости вращения + быстродействие
        mModelCtrlVeloQuiq = QCtrlVeloQuiq(mDriverModel.mElectMotor,mDriverModel.mpLoad
                              ,marrSpreadParams, mMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective
                                           ,mQuiqMaxU0,constDbl_mh);
        mpModelCtrl = &mModelCtrlVeloQuiq;

        mRealCtrlVeloQuiq = QCtrlVeloQuiq (mDriverReal.mElectMotor,mDriverReal.mpLoad
                      ,marrSpreadParams, mRealMomOut, mTettaBegin, mOmegaBegin, 0.
                                           ,marrObjective,mQuiqMaxU0,constDbl_mh);
        mpRealCtrl = &mRealCtrlVeloQuiq;
        break;

    case 2:
        mModelCtrlPos= QCtrlPos(mDriverModel.mElectMotor,mDriverModel.mpLoad
                     ,marrSpreadParams, mMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpModelCtrl = &mModelCtrlPos;

        mRealCtrlPos= QCtrlPos(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpRealCtrl = &mRealCtrlPos;
        break;
    case 3:
        mModelCtrlFollow = QCtrlFollow(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpModelCtrl = &mModelCtrlFollow;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpRealCtrl = &mRealCtrlFollow;
        break;
    case 4:
        mModelCtrlDnmkCorrect_M_R = QCtrlDnmkCorrect_M_R(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh
                         ,mCmpArrLamb, mTPeriodParamsRenew, 0.);
        mpModelCtrl = &mModelCtrlDnmkCorrect_M_R;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpRealCtrl = &mRealCtrlFollow;


        break;
    case 5:
        mModelCtrlDnmkCorrect_W = QCtrlDnmkCorrect_W(mDriverModel.mElectMotor,mDriverModel.mpLoad
                        ,marrSpreadParams, mMomOut, mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh
                        , mTPeriodParamsRenew, 0.);
        mpModelCtrl = &mModelCtrlDnmkCorrect_W;

        mRealCtrlFollow = QCtrlFollow(mDriverReal.mElectMotor,mDriverReal.mpLoad
                     ,marrSpreadParams, mRealMomOut,mTettaBegin, mOmegaBegin, 0.,marrObjective,constDbl_mh);
        mpRealCtrl = &mRealCtrlFollow;
        break;
    default:
        break;
    };


    QStatSolutionParams StatSolutionParams;
    mpModelCtrl->calcStationaryParams(marrObjective,mCmpArrLamb, &StatSolutionParams);
    mpModelCtrl->mStatSolutionParams = StatSolutionParams;
    (*mpRealCtrl).mStatSolutionParams = StatSolutionParams;


    TComp cmparrRealEigenValues[4];
    mpRealCtrl->CalcEigenValues(marrObjective[0], marrObjective[1],  StatSolutionParams.marrGears,cmparrRealEigenValues );
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

   mDriveMoveImit = QDriveMoveImit (mpModelCtrl,mDriverReal,mMeasurmentImitator );

   //const double VAlIntegratorTime = 2;

   mDriveMoveImit.move(0., mMovingT, VAlIntegrStep
                    ,mTimeIntegrator ,arrBuff, &quantDoneSteps);

   ui->tableWidget->item(0,1)->setText(QString::number(mDriveMoveImit.mpCtrl->mElectMotor.mResist));
   ui->tableWidget->item(0,3)->setText(QString::number(mDriveMoveImit.mpCtrl->mMomOut));


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
   double arrSyst[QVARS] = {0.};
   // веутор диспесий M(x-xsyst)*(x-xsyst)
   double arrDisp[QVARS] = {0.};
   // ПЕРИОД гармонической составляющей остат момента
   const double VAl_TPeriod = 2. * M_PI /(mDriverReal.mElectMotor.mZp * marrObjective[1])*2. ;
           //(mDriverReal.mElectMotor.mZp * marrObjective[1] / 2./M_PI) * 3.;
   mDriveMoveImit.StatisticProcess(arrBuff, quantDoneSteps,VAl_TPeriod, arrSyst, arrDisp);


   ui->tableWidget_4->item(1,6)->setText(QString::number((arrSyst[3] * 1000.)));
   ui->tableWidget_4->item(1,7)->setText(QString::number((sqrt(arrDisp[3]) * 1000.)));

   ui->tableWidget_4->item(1,8)->setText(QString::number((arrSyst[0] * 1000.)));
   ui->tableWidget_4->item(1,9)->setText(QString::number((sqrt(arrDisp[0]) * 1000.)));

   ///
  // анализ разбросов траектории в стационарном режиме аналитическим методом
   double arrW[2] = {0.};
   double arr_dx_po_dAlf[QVARS * QPARAMS] = {0.}, arrAHarmSquare[QVARS ] = {0.}, arrSystErr[QVARS ] = {0.}
         , arrAHarm[QVARS ] = {0.}, arrSumSig[QVARS ] = {0.}, arrRandomSig[QVARS ] = {0.};

   mDriveMoveImit.calcAnaliticSumTrajScatters(arrW, arr_dx_po_dAlf, arrAHarm, arrSystErr
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


 // построение графика ошибки по Omega
 double *parrx = new double[quantDoneSteps];
 double *parry = new double[quantDoneSteps];


 for (int i = 0; i < quantDoneSteps; ++i)
 {
     double *p = &(arrBuff[i * QUantColsCSVReport ]);
     parrx[i]= p[0] * 100.;
     parry[i] = (p[1]  - marrObjective[1]) * 1000000.;

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
     parry[i] = (p[4]  - marrObjective[0] - marrObjective[1] *(p[0] - valTStart)) * 1000000.;

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

*/
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





void MainWindow::on_pushButton_4_clicked()
{
      inputData();
      wchar_t wchOutPutFold0[400] = {0};
      wcscpy(wchOutPutFold0, mwchOutPutFold);

      QSegment SgmBarrier;
      TPlane Plane;

      QBrlsBody BrlsBody( mCv, mCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                              ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane,mcs, 0);

      double valMass = BrlsBody.calcMass();
      double valJY0 = BrlsBody.calcJY0();
      ui->doubleSpinBox_5->setValue(valMass);
      ui->doubleSpinBox_4->setValue(valJY0);


      const double VAlStep = 0.001;
      const int QUantRows = int(BrlsBody.mLength/ 2./VAlStep ) + 1;

      double *arrBuff = new double[QUantColsReport_Wind  * QUantRows];
      memset(arrBuff, 0, sizeof(double) * QUantColsReport_Wind *(QUantRows ));
      int quantDoneSteps = -1;

      double valFSum = 0., valMSum = 0.;
      BrlsBody.analysisWind(mOmega, mAccellerOmega, mWindV,VAlStep, &valFSum, &valMSum,arrBuff, &quantDoneSteps);

      ui->doubleSpinBox_7->setValue(valFSum);
      ui->doubleSpinBox_6->setValue(valMSum);
      // вывод графиков
      int iLenName = 30;
      wchar_t *pwcharrColNames = new wchar_t[QUantColsReport_Wind * iLenName];
      memset(pwcharrColNames, 0, QUantColsReport_Wind * iLenName*sizeof(wchar_t));

      wcscpy(pwcharrColNames, L"z");
      wcscpy(&pwcharrColNames[iLenName   ],  L"q");
      wcscpy(&pwcharrColNames[iLenName *2],  L"Qq");
      wcscpy(&pwcharrColNames[iLenName *3],  L"Mq");
      wcscpy(&pwcharrColNames[iLenName *4],  L"MaxSig_0_DIV_[Sig]");
      wcscpy(&pwcharrColNames[iLenName *5],  L"MaxSig_1_DIV_[Sig]");
      wcscpy(&pwcharrColNames[iLenName *6],  L"MaxSig_2_DIV_[Sig]");
      wcscpy(&pwcharrColNames[iLenName *7],  L"MaxSig_3_DIV_[Sig]");
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

      arrScale[4] =arrScale[5] =arrScale[6] =arrScale[7] = 10000.;

    for (int i =1; i < QUantColsReport_Wind; i++)
    {
      TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                      ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                      ,QUantColsReport_Wind  // - к-во переменных о корорых накоплена информация в буфере
                                      ,quantDoneSteps //  - к-во точек
                                      ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                      ,iLenName // максимальная длина имени переменной
                                      ,0  //  номер переменной по оси X
                                      ,i  //  номер переменной по оси Y
                                      ,1000.  //  масштаб по оси X
                                      ,arrScale[i]  // масштаб по оси Y
                                       ) ;
    }
    BrlsBody.createSideView(1000. ,1000. ,wchOutPutFold0);

    double  valZDangerous = -1.,  valCoeff = -1.;
    int iDangerousSection = -1;
    for(int i =0; i < quantDoneSteps; ++i)
    {
        for (int j =0; j < 4; ++j)
        {
            if (arrBuff[i * QUantColsReport_Wind + 4 + j] > valCoeff)
            {
              valCoeff =  arrBuff[i * QUantColsReport_Wind + 4 +j];
              valZDangerous = arrBuff[i * QUantColsReport_Wind];
              iDangerousSection = j;
            }
        }

    }

    ui->doubleSpinBox_9->setValue(valZDangerous);
    ui->doubleSpinBox_8->setValue(valCoeff);
    ui->spinBox->setValue(iDangerousSection);
    delete []arrBuff;
    delete []pwcharrColNames;

}

void MainWindow::on_pushButton_5_clicked()
{
    inputData();
    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);

    QSegment SgmBarrier;
    TPlane Plane;

    QBrlsBody BrlsBody( mCv, mCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                            ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane, mcs,0);



    const double VAlStep = 0.001;
    const int QUantRows = int(BrlsBody.mLength/ 2./VAlStep ) + 1;

    double *arrBuff = new double[QUantColsReport_Wind  * QUantRows];
    memset(arrBuff, 0, sizeof(double) * QUantColsReport_Wind *(QUantRows ));
    int quantDoneSteps = -1;

    double valFSum = 0., valMSum = 0.;
    BrlsBody.analysisIce(mPIce,VAlStep,  &valFSum, &valMSum,arrBuff, &quantDoneSteps);

    ui->doubleSpinBox_13->setValue(valFSum);
    ui->doubleSpinBox_12->setValue(valMSum);
    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[QUantColsReport_Wind * iLenName];
    memset(pwcharrColNames, 0, QUantColsReport_Wind * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"z");
    wcscpy(&pwcharrColNames[iLenName   ],  L"q");
    wcscpy(&pwcharrColNames[iLenName *2],  L"Qq");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Mq");
    wcscpy(&pwcharrColNames[iLenName *4],  L"MaxSig_0_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *5],  L"MaxSig_1_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *6],  L"MaxSig_2_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *7],  L"MaxSig_3_DIV_[Sig]");
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
    arrScale[1] = -10.;
    arrScale[2] = arrScale[3] = 10.;
    arrScale[4] =arrScale[5] =arrScale[6] =arrScale[7] = 100000.;

  for (int i =1; i < QUantColsReport_Wind; i++)
  {
    TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,QUantColsReport_Wind  // - к-во переменных о корорых накоплена информация в буфере
                                    ,quantDoneSteps //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1000.  //  масштаб по оси X
                                    ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }
  BrlsBody.createSideView(1000. ,1000. ,wchOutPutFold0);

  double  valZDangerous = -1.,  valCoeff = -1.;
  int iDangerousSection = -1;
  for(int i =0; i < quantDoneSteps; ++i)
  {
      for (int j =0; j < 4; ++j)
      {
          if (arrBuff[i * QUantColsReport_Wind + 4 + j] > valCoeff)
          {
            valCoeff =  arrBuff[i * QUantColsReport_Wind + 4 +j];
            valZDangerous = arrBuff[i * QUantColsReport_Wind];
            iDangerousSection = j;
          }
      }

  }

  ui->doubleSpinBox_15->setValue(valZDangerous);
  ui->doubleSpinBox_14->setValue(valCoeff);
  ui->spinBox_2->setValue(iDangerousSection);
  delete []arrBuff;
  delete []pwcharrColNames;
}

void MainWindow::on_pushButton_6_clicked()
{
    inputData();
    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);

    QSegment SgmBarrier;
    TPlane Plane;

    QBrlsBody BrlsBody( mCv, mCx, VAL_BRLS_MASS, VAL_BRLS_LENGTH, VAL_BRLS_HEIGHT
                            ,VAL_BRLS_dOm_po_dt, SgmBarrier,Plane,mcs,0);

     double valHail_Mass =  QBrlsBody::calc_Hail_Mass(mHailDiam);
     double valHail_V =  QBrlsBody::calc_Hail_V(mHailDiam);
     double valHail_H = QBrlsBody::calc_Hail_EffectH(mHailDiam);
     ui->doubleSpinBox_17->setValue(valHail_V);
     ui->doubleSpinBox_16->setValue(valHail_H );
     ui->doubleSpinBox_22->setValue(valHail_Mass * 1000. );



    const double VAlStep = 0.01;
    const int QUantRows = int(BrlsBody.mLength/ 2./VAlStep ) + 1;

    double *arrBuff = new double[QUantColsReport_Wind  * QUantRows];
    memset(arrBuff, 0, sizeof(double) * QUantColsReport_Wind *(QUantRows ));
    int quantDoneSteps = -1;

    double valFSum = 0., valMSum = 0.;
    BrlsBody.analysisHail_(mHailDiam,mHailZ, VAlStep,  &valFSum, &valMSum,arrBuff, &quantDoneSteps);

    ui->doubleSpinBox_23->setValue(valFSum);
    ui->doubleSpinBox_24->setValue(valMSum);
    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[QUantColsReport_Wind * iLenName];
    memset(pwcharrColNames, 0, QUantColsReport_Wind * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"z");
    wcscpy(&pwcharrColNames[iLenName   ],  L"Mq");
    wcscpy(&pwcharrColNames[iLenName *2],  L"Omega2");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Progib");
    wcscpy(&pwcharrColNames[iLenName *4],  L"MaxSig_0_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *5],  L"MaxSig_1_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *6],  L"MaxSig_2_DIV_[Sig]");
    wcscpy(&pwcharrColNames[iLenName *7],  L"MaxSig_3_DIV_[Sig]");
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
    arrScale[3] =arrScale[2] = 100000.;
    arrScale[4] =arrScale[5] =arrScale[6] =arrScale[7] = 1000.;

  for (int i =1; i < QUantColsReport_Wind; i++)
  {
    TYrWriteShapeFile::WriteOneReport(wchOutPutFold0// путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,QUantColsReport_Wind  // - к-во переменных о корорых накоплена информация в буфере
                                    ,quantDoneSteps //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1000.  //  масштаб по оси X
                                    ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }
  BrlsBody.createSideView(1000. ,1000. ,wchOutPutFold0);

  double  valZDangerous = -1.,  valCoeff = -1.;
  int iDangerousSection = -1;
  for(int i =0; i < quantDoneSteps; ++i)
  {
      for (int j =0; j < 4; ++j)
      {
          if (arrBuff[i * QUantColsReport_Wind + 4 + j] > valCoeff)
          {
            valCoeff =  arrBuff[i * QUantColsReport_Wind + 4 +j];
            valZDangerous = arrBuff[i * QUantColsReport_Wind];
            iDangerousSection = j;
          }
      }

  }

  ui->doubleSpinBox_19->setValue(valZDangerous);
  ui->doubleSpinBox_18->setValue(valCoeff);
  ui->spinBox_3->setValue(iDangerousSection);
  delete []arrBuff;
  delete []pwcharrColNames;
}

void MainWindow::on_pushButton_7_clicked()
{
    inputData();
    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mwchOutPutFold);


    const TURPointXY pointCentre(0.,0.);
    const double valR = 0.2;
    const double valFi0 = 0.;//M_PI /2.;
    const double valFi1 = -3./4.*M_PI - M_PI /2.;
    const int NPoints = 1000;
    TURPolyLine plnSEctor = TURPolyLine::fncCreateSector( pointCentre, valR,
                         valFi0, valFi1,NPoints) ;



    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  wchOutPutFold0);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.) ;

     wcscpy(  wchAxesFileName0,  wchOutPutFold0);
     wcscat(wchAxesFileName0, L"\\Sector.shp");
     plnSEctor.WriteSetSHPFiles(wchAxesFileName0, &plnSEctor, 1);

     const TURPointXY pointBegin(valR * cos(valFi1), valR * sin(valFi1));
     double valFi2 = valFi1 -M_PI/2.;
     double vald = 0.05;
     const TURPointXY pointEnd(pointBegin.X + vald *cos(valFi2),pointBegin.Y + vald *sin(valFi2));
             const double valLength = 0.05;
     const double valAng = M_PI/10.;

     TURPolyLine plnArrow = TURPolyLine::fncCreateArrow(pointBegin,  pointEnd
                             , valLength,valAng);
     wcscpy(  wchAxesFileName0,  wchOutPutFold0);
     wcscat(wchAxesFileName0, L"\\Arrow.shp");
     plnArrow.WriteSetSHPFiles(wchAxesFileName0, &plnArrow, 1);
}
