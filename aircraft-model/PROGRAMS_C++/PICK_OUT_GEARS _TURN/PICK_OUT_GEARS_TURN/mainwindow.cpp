#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <float.h>
#include <wchar.h>
#include <math.h>
#include <string.h>
#include <QFileDialog>
#include <QMessageBox>
//#include <QDir>
#include <stdlib.h>

#include "YrRead.h"
//#include "Plane.h"
//#include "Operation.h"
//#include "OperationNew.h"

//#include "Blade.h"
#include "Environment.h"
#include "YrWriteShapeFile.h"
#include "URPolygon.h"
#include "UrPointXY.h"
#include "dir.h"
#include "Comp.h"
#include "Equations.h"
#include "Helic.h"
#include "MatrixProccess.h"
#include "YrWrite.h"
#include "HelicConstants.h"
#include "BallanceCalc.h"
#include "PartHelicTraj.h"



#include <engine.h>
#include "liprog.h"

#define NOT_POSSIBLE_VALUE -1000000000.

extern long double constHelicMass;

// НОМЕР 6000

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    /*mQuantActualTaskPoints = 2;
    memset(marrCoord, 0, QUANT_FLY_TASK_POINTS * 3 * sizeof(double));
    marrCoord[3]= 3000.;
    marrCoord[5]= 1000.;
    for (int i = 6; i < QUANT_FLY_TASK_POINTS * 3; i++)
    {
     marrCoord[i] = NOT_POSSIBLE_VALUE;
    }

*/

 ui->progressBar ->setValue ( 1 );
 ui->progressBar->setOrientation(Qt::Horizontal);
     ui->progressBar->setRange(0, 100); // Let's say it goes from 0 to 100
     ui->progressBar->setValue(0); // With a current value of 10




    // начальная широта и долгота
  //  mLong0 = 30.; // долгота
  //  mLat0 = 30.;

    // полетное время
    mTFly = 1000.;


    ///



  // атмосфера
   mWindHor = 0.; // горизонтьальная скорость ветра
   mWindCourse = 0.;// курсовой угол вектора горизонтьальной скорости ветра
   mWindVert = 0.;// вертикальная скорость ветра полож направление вверх
   mTemperature0 = 18.;

   mparrCandidates = NULL;
  ui->comboBox->setCurrentIndex(0);
  ui->comboBox->setCurrentIndex(3);

}

MainWindow::~MainWindow()
{
    if(mparrCandidates)
    {
        free (mparrCandidates);
    }
    delete ui;
}






//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
void MainWindow::on_pushButton_clicked()
{
    inputData();


    int quantVars = -1;
    // mpPartHelicTraj->get_QuantOfControlledVarsTang(&quantVars);
   int maxQuant = 1000000;
   double *arrDataBuf = new double [maxQuant * (QUantCurNZSKVarsVS - 1)*4];

 memset(arrDataBuf, 0, maxQuant * (QUantCurNZSKVarsVS - 1) * 4 *sizeof( double));


    int quantRows= 0;




  mpPartHelicTraj->selectGearsTang( mz1 ,  mz2, mval_A,arrDataBuf ,maxQuant, &quantRows);

 // запись массива информации в  CSV файл
 // INPUT:
 // FileName
 // parrBuff[ iNumRows * iNumCols] - массив с информацией
 // iNumRows- к-во строк массива
 // iNumCols - к-во столбцов массива
 // pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
 // pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
 TYrWrite::WriteMassiveInFIleSCV(mpwchPeredChislaSCVFIle,arrDataBuf, quantRows, (QUantCurNZSKVarsVS - 1)*4
                              ,NULL,NULL, 30);

 delete []arrDataBuf ;
 ui->doubleSpinBox_4  ->setValue(((double)quantRows));
}

//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------

void MainWindow::on_doubleSpinBoxVHor_valueChanged(double arg1)
{
    mWindHor = arg1;
}


void MainWindow::on_doubleSpinBoxCourseAng_valueChanged(double arg1)
{
    mWindCourse = arg1* M_PI/ 180.;
}

void MainWindow::on_doubleSpinBoxVVert_valueChanged(double arg1)
{
    mWindVert = arg1;
}



void MainWindow::on_doubleSpinBoxTemperature0_valueChanged(double arg1)
{
    mTemperature0 = arg1;
}

void MainWindow::inputData()
{
    TEnvironment mEnvironment (0, 0., 0., 18. );
    const bool bFullType = false,  bTypeofFullRotor = false;
    THelic mHelic( bFullType, bTypeofFullRotor,constHelicMass);


    mvalYaw =  ui->doubleSpinBox_2->value();
    mvalRad =  ui->doubleSpinBox->value();
    mvalVx = ui->doubleSpinBox_19->value();
    mvalY = ui->doubleSpinBox_20->value();
    m_dPsi_po_dt = ui->doubleSpinBox_3->value();
    //mPsi0 = ui->doubleSpinBox_5->value();





   long double arrRotCentre_CurNZSK[3] ={0.};
   const long double  TBegin = 0.;



    switch (this->ui->comboBox->currentIndex())
    {
    case 0: // прямолин
         mLineMove = TLineMove( mHelic, mEnvironment , mvalY
                                        , mvalVx,  mvalYaw, 0.);
         mpPartHelicTraj = &mLineMove;
         memset(mpPartHelicTraj->marrPhaseVect, 0, sizeof(long double)* (QUantCurNZSKVarsVS-1));
         memcpy(mpPartHelicTraj->marrPhaseVect, mLineMove.marrSteadySolution0,  sizeof(long double)* (QUantCurNZSKVarsVS-1));

        break;

    case 1: // вираж
        arrRotCentre_CurNZSK[0] = 0.;
        arrRotCentre_CurNZSK[1] = mvalY;
        arrRotCentre_CurNZSK[2] = mvalRad;
        mTurnMove = TTurnMove(mHelic, mEnvironment,mvalY
                   ,mvalVx,  mvalYaw, mvalRad,marrRotCentre_CurNZSK
                   ,TBegin, M_PI/2.);

        mpPartHelicTraj = &mTurnMove;
        memset(mpPartHelicTraj->marrPhaseVect, 0, sizeof(long double)* (QUantCurNZSKVarsVS-1));
        memcpy(mpPartHelicTraj->marrPhaseVect, mTurnMove.marrSteadySolution0,  sizeof(long double)* (QUantCurNZSKVarsVS-1));


        break;

    case 2: // висение

        mHover = THover( mHelic, mEnvironment , mvalY
                         ,   10.,TBegin);
        mpPartHelicTraj = &mHover;
        memset(mpPartHelicTraj->marrPhaseVect, 0, sizeof(long double)* (QUantCurNZSKVarsVS-1));
        memcpy(mpPartHelicTraj->marrPhaseVect, mHover.marrSteadySolution0,  sizeof(long double)* (QUantCurNZSKVarsVS-1));

        break;

    case 3: // вращение

        mRotating = TRotating(mHelic, mEnvironment,mvalY
                   ,m_dPsi_po_dt,TBegin, 100., 0.);

        mpPartHelicTraj = &mRotating;
        memset(mpPartHelicTraj->marrPhaseVect, 0, sizeof(long double)* (QUantCurNZSKVarsVS-1));
        memcpy(mpPartHelicTraj->marrPhaseVect, mRotating.marrSteadySolution0,  sizeof(long double)* (QUantCurNZSKVarsVS-1));

        break;

    default:
        break;
    }

    mpPartHelicTraj->doAirDensity();







    // параметры лемнитскаты (овала Кассини):

    mLemn_a = ui->doubleSpinBox_a->value();
    mLemn_c = ui->doubleSpinBox_c->value();
    mLemn_l = ui->doubleSpin_L->value();
    ///
    /// \brief wchOutPutFold0




     mz1 = -mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l;
     mz2 = mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l;
     mval_A = mLemn_a * mLemn_a;


}





void MainWindow::on_doubleSpinBoxVHor_valueChanged(const QString &arg1)
{

}

void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{

}

void MainWindow::on_pushButton_4_clicked()
{   


    int quanTitleRows = 0;
    int quantCols = TYrRead::calcColCountFromCSV(mpwchPickedOutCandidates );
    mnumCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;




    mparrCandidates = (double*)malloc(mnumCandidates * quantCols * sizeof(double));
    int nrows = 0;
      TYrRead::YrReadTabCSV_1(mpwchPickedOutCandidates// файл с таблицей
          ,quanTitleRows // к-во заголов строк
          ,0 // к-во заголов cтолбцов
          ,mnumCandidates // к-во строк содержательной части таблицы
          ,quantCols// к-во столбцов содержательной части таблицы
          ,&nrows  // к-во прочитанных строк
          ,mparrCandidates  // массив в который принимается информация
          )  ;
      ///
 quantCols = (QUantCurNZSKVarsVS -1) * 4;

inputData();
   int iNumVars = -1;
   int iNumRuls = -1;

   int iarrNumsControlledVars[QUantCurNZSKVarsVS -1] = {0};

   int iarrNumsControls[4] = {0};


   if(mbTang)
   {    // канал тангажа

       mpPartHelicTraj->get_arrayOfControlledVarsTang(&iNumVars, iarrNumsControlledVars);
       iarrNumsControls[0] = 0;
       iarrNumsControls[1] = 1;
       iNumRuls = 2;
   }
   else
   {
       mpPartHelicTraj->get_arrayOfControlledVars(&iNumVars, iarrNumsControlledVars);

       iarrNumsControls[0] = 0;
       iarrNumsControls[1] = 1;
       iarrNumsControls[2] = 2;
       iarrNumsControls[3] = 3;
       iNumRuls = 4;

   }


  // ФОРМИРОВАНИЕ МАТРИЦ ДЛЯ ПОЛЛНОЙ ЗАДАЧИ

 long double *arr_dF_po_dW = (long double *)malloc(iNumVars * iNumRuls *sizeof(long double ));
 long double *arr_dF_po_dx = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));
 long double *arr_C = (long double *)malloc(iNumRuls * iNumVars *sizeof(long double ));

 long double *arrT0 = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));
 long double *arrT1 = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));



 if( mbTang)
 {
   mpPartHelicTraj->fill_df_po_px_and_df_po_dW_TangCanale( arr_dF_po_dx, arr_dF_po_dW);


 }
 else
 {
   mpPartHelicTraj->fill_df_po_px_and_df_po_dW( arr_dF_po_dx, arr_dF_po_dW);

 }
 int maxQuant = 1000000;
 int quantRows = 0;

 double *arrDataBuf = (double *)malloc(maxQuant * quantCols * 4 *sizeof(double));
 for (int i =0; i < mnumCandidates; ++i)
 {

     TPartHelicTraj::extractLocalGearMtrx_From_TotalGearMtrx(&(mparrCandidates [i * quantCols]),iNumVars, iarrNumsControlledVars
                                                                 ,iNumRuls, iarrNumsControls
                                                                 ,arr_C );


     MtrxMultMatrx(arr_dF_po_dW, iNumVars, iNumRuls, arr_C, iNumVars, arrT0) ;

     MtrxSumMatrx(arr_dF_po_dx,arrT0 ,iNumVars, iNumVars, arrT1) ;

     if(TPartHelicTraj::IsStability_( mz1 ,  mz2, mval_A,arrT1, iNumVars))
     {
         memcpy(&(arrDataBuf [ quantCols * quantRows]),&(mparrCandidates [i * quantCols]), sizeof(double) *quantCols);

         quantRows++;
         if (quantRows== maxQuant)
         {
             break;
         }
     }

     int iprogress =(100. * ((double)i) /((double)mnumCandidates)) ;
     ui->progressBar->setValue(iprogress );
 }

 if (0 == quantRows)
 {
     free(arrDataBuf) ;
     free(arr_dF_po_dW);
     free(arr_dF_po_dx);
     free(arr_C);

     free(arrT0);
     free(arrT1);

     free(mparrCandidates);
     mparrCandidates = NULL;
     QString  qstr= QString::number(quantRows);
     ui->lineEdit_7->setText(qstr);
     return;
 }


 // запись массива информации в  CSV файл
 // INPUT:
 // FileName
 // parrBuff[ iNumRows * iNumCols] - массив с информацией
 // iNumRows- к-во строк массива
 // iNumCols - к-во столбцов массива
 // pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
 // pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
 TYrWrite::WriteMassiveInFIleSCV(mpwchPeredChislaSCVFIle,arrDataBuf, quantRows, quantCols
                              ,NULL,NULL, 30);

 free(arrDataBuf) ;
 free(arr_dF_po_dW);
 free(arr_dF_po_dx);
 free(arr_C);

 free(arrT0);
 free(arrT1);




 free(mparrCandidates);
 mparrCandidates = NULL;

 QString  qstr= QString::number(quantRows,10);
 ui->lineEdit_7->setText(qstr);

}

void MainWindow::on_doubleSpinBoxMass_valueChanged(const QString &arg1)
{

}

void MainWindow::on_lineEdit_cursorPositionChanged(int arg1, int arg2)
{

}

void MainWindow::on_pushButton_5_clicked()
{
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с ПЧ малой1 задачи","D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY" ,"*.csv");
    this->ui->lineEdit_3->setText(strFold);

    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
    int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;
    ui->doubleSpinBox_21  ->setValue((double)numCandidates );

    /*
     QString strFold = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с кандидатами  в ПЧ", "D:\\REPOSITORIES\\aircraft-model\\СИСТЕМА_ДИФ_УРАВНЕНИЙ_ВЕРТОЛЕТА\\InputFile.csv","*.csv");
    this->ui->lineEdit_3->setText(strFold);

    // wchar_t array [400] = {0};
    // strFold.toWCharArray(array);
    // array[strFold.length()] = 0;
  //  wcscpy(mpwchPickedOutCandidates, array);

     mnumCandidates  = TYrRead::YrCalcRows(strFold ) -1;
    if(mparrCandidates)
    {
        free (mparrCandidates);
    }
  mparrCandidates = (double*)malloc(mnumCandidates * 5 * sizeof(double));
  // int numRows = TYrRead::YrCalcRows(wchar_t *FileName)
  int nrows = 0;
    TYrRead::YrReadTabCSV_1(strFold// файл с таблицей
        ,1 // к-во заголов строк
        ,0 // к-во заголов cтолбцов
        ,mnumCandidates // к-во строк содержательной части таблицы
        ,5// к-во столбцов содержательной части таблицы
        ,&nrows  // к-во прочитанных строк
        ,mparrCandidates  // массив в который принимается информация
        )  ;

  */
}

void MainWindow::on_pushButton_6_clicked()
{
   inputData();

    ///
    wchar_t wchOutPutFold0[400];
    wcscpy(wchOutPutFold0, mpwchPeredChislaSCVFIle);











    //-------------------------------------------------------------

    mnumCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;


   int quantCols = (QUantCurNZSKVarsVS-1) *4;
 mparrCandidates = (double*)malloc(mnumCandidates * quantCols * sizeof(double));
 // int numRows = TYrRead::YrCalcRows(wchar_t *FileName)
 int nrows = 0;
   TYrRead::YrReadTabCSV_1(mpwchPickedOutCandidates// файл с таблицей
       ,0 // к-во заголов строк
       ,0 // к-во заголов cтолбцов
       ,mnumCandidates // к-во строк содержательной части таблицы
       ,quantCols// к-во столбцов содержательной части таблицы
       ,&nrows  // к-во прочитанных строк
       ,mparrCandidates  // массив в который принимается информация
       )  ;
   // заготовка матрицы передаточных чисел
   int numSet = int(ui->doubleSpinBox_22->value() + 0.1) -1;
   double arrC[(QUantCurNZSKVarsVS-1) *4] = {0.};
   for (int i =0; i < quantCols;++i)
   {
       arrC[i] = mparrCandidates[ i +numSet *quantCols];//?????????????????????????
   }

   free(mparrCandidates);
   mparrCandidates = NULL;
   //-----------------------------------------------------------------------------------------
   int quantVars = -1;
   // mpPartHelicTraj->get_QuantOfControlledVarsTang(&quantVars);
  int maxQuant = 1000000;
  double *arrDataBuf = new double [maxQuant * (QUantCurNZSKVarsVS - 1)*4];

memset(arrDataBuf, 0, maxQuant * (QUantCurNZSKVarsVS - 1) * 4 *sizeof( double));


   int quantRows= 0;




 mpPartHelicTraj->selectGears( mz1 ,  mz2, mval_A,arrC, arrDataBuf ,maxQuant, &quantRows);

// запись массива информации в  CSV файл
// INPUT:
// FileName
// parrBuff[ iNumRows * iNumCols] - массив с информацией
// iNumRows- к-во строк массива
// iNumCols - к-во столбцов массива
// pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
// pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
TYrWrite::WriteMassiveInFIleSCV(mpwchPeredChislaSCVFIle,arrDataBuf, quantRows, (QUantCurNZSKVarsVS - 1)*4
                             ,NULL,NULL, 30);

delete []arrDataBuf ;
QString  qstr= QString::number(quantRows,10);
ui->lineEdit_5->setText(qstr);

/*
    // ФОРМИРОВАНИЕ МАТРИЦ ДЛЯ ПОЛНОЙ ЗАДАЧИ

   int iarrNumsControlledVars[QUantCurNZSKVarsVS -1] = {0};
   int iarrNumsControls[4] = {0};
   int iNumVars = -1;
   int iNumRuls = 4;
   mpPartHelicTraj->get_arrayOfControlledVars(&iNumVars, iarrNumsControlledVars);

   iarrNumsControls[0] = 0;
   iarrNumsControls[1] = 1;
   iarrNumsControls[2] = 2;
   iarrNumsControls[3] = 3;


   long double *arr_dF_po_dW = (long double *)malloc(iNumVars * iNumRuls *sizeof(long double ));
   long double *arr_dF_po_dx = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));
   long double *arr_C = (long double *)malloc(iNumRuls * iNumVars *sizeof(long double ));

   long double *arrT0 = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));
   long double *arrT1 = (long double *)malloc(iNumVars * iNumVars *sizeof(long double ));





     long double arr_BigA0[11 * 11] = {0.}, arr_BigB0[11*4] = {0.};;
   mpPartHelicTraj->fill_df_po_px_and_df_po_dW(arr_BigA0, arr_BigB0);


    // заготовка матрицы передаточных чисел
    //long double arrBigС0[4 * 11] = {0.};
    long double arrC00[44] = {0.};


    arrC00[2] = arrC[0];
    arrC00[7] = arrC[4];
    arrC00[9] = arrC[3];
    arrC00[11]= arrC[6];
    arrC00[14] =arrC[7] ;

    arrC00[3] = arrC[2];
    arrC00[0] = arrC[1];
    arrC00[13] = arrC[5];
    arrC00[20] = arrC[8];
    arrC00[18] = arrC[9];



    long double arrT00[121] = {0.}, arrT10[44] = {0.};



    int quantRows= 0;
    int maxQuant = 1000000;
    double *arrDataBuf = new double [maxQuant * 44];
    bool bend = false;

    //,0.000000		,-0.008833	,0.000000	,0.000000	,0.000000	,0.000000	,0.238000	,0.000000	,0.406500	,0.000000
   // ,-0.003680	,0.000000	,-0.016900	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000	,0.000000
  //  ,0.000000		,0.000000	,0.000000	,-0.101531	,-0.150000	,0.000000	,0.000000	,-0.500000	,0.000000	,0.000000
   // ,0.000000		,0.000000	,0.000000	,0.000000	,0.000000	,-0.050000	,0.000000	,0.000000	,0.000000	,-0.051000
    for (int i7 = 0; i7 < 10; ++i7)// 4
    {
          arrC00[39] = -0.0201 - 0.01 * ((long double) i7 ); // это OmY - DelFi
        // arrC00[39] = -0.01 - 0.01 * ((long double) i7 );
      for (int i6 =0; i6 < 20; ++i6)
      {
          arrC00[43] = -0.001- 0.005 * ((long double) i6 );// это Psi - DelFi
      for (int i0 = 0; i0 < 20; ++i0)// i0=8, 170,250, 0.1
      {
        arrC00[27] = -0.05 - 0.025 * ((long double) i0 ); // это OmX- Etta
         // arrC00[27] = -0.2 - 0.05 * ((long double) i0 );
         //  arrC00[27] = -0.02 - 0.01 * ((long double) i0 );
          for (int i1 = 0; i1< 40; ++i1)
          {
           // arrC00[30] = -0.05 - 0.05 * ((long double) i1 );
           arrC00[30] = -0.01 - 0.05 * ((long double) i1 ); // это Gamma -Etta

                     for (int i4 =0; i4 < 10; ++i4)
                     {


                       arrC00[26] = - 0.001 -0.06* ( static_cast<long double>( i4 )) ;   // это Vz - Etta
                       //  arrC00[26] = -0.001 - 0.002 * ((long double)i4 + 1.);
                       for (int i5 =0; i5 < 10; ++i5)
                       {


                          arrC00[23] =   - 0.001 -0.06* ( static_cast<long double>( i5)) ;   // это Z -Etta
                         // arrC00[23] =-0.0001 - 0.0002 *((long double) i5 );
                          MtrxMultMatrx(arr_BigB0, 11, 4, arrC00, 11, arrT00) ;

                          MtrxSumMatrx(arr_BigA0,arrT00 ,11, 11, arrT10) ;
                          if(TTurnMove::IsStability_( mz1 ,  mz2, mval_A,arrT10, 11))
                          {
                              for (int j=0; j < 44; ++j)
                              {
                                arrDataBuf [ 44 * quantRows +j] = arrC00[j] ;
                              }


                          quantRows++;
                          if (quantRows== maxQuant)
                          {

                          bend = true;
                          break;
                          }

                         }
                         else
                         {
                          int iii = 0;
                          }

                     }
                       if(bend)
                       {
                           break;
                       }
                  }



            if(bend)
            {
                break;
            }
          }
          if(bend)
          {
              break;
          }
      }
          if(bend)
          {
              break;
          }

      }
          if(bend)
          {
              break;
          }

   }



 // запись массива информации в  CSV файл
 // INPUT:
 // FileName
 // parrBuff[ iNumRows * iNumCols] - массив с информацией
 // iNumRows- к-во строк массива
 // iNumCols - к-во столбцов массива
 // pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
 // pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
 TYrWrite::WriteMassiveInFIleSCV(mpwchPeredChislaSCVFIle,arrDataBuf, quantRows, 44
                              ,NULL,NULL, 30);

 delete []arrDataBuf ;
 QString  qstr= QString::number(quantRows,10);
 ui->lineEdit_5->setText(qstr);
*/
//0.000000	0.000000	-0.011217	0.000000	0.000000	0.000000	0.000000	0.296557	0.000000	0.593113	0.000000	-0.001532	0.000000	0.000000	-0.011781	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	-0.009560	0.000000	0.000000	-0.050987	-0.050000	0.000000	0.000000	-0.300000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	0.000000	-0.050000	0.000000	0.000000	0.000000	-0.050000

}

void MainWindow::on_pushButton_7_clicked()
{
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с ПЧ малой1 задачи","D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY" ,"*.csv");
    this->ui->lineEdit_6->setText(strFold);

    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
    int quantCols = TYrRead::calcColCountFromCSV(mpwchPickedOutCandidates ) ;

    int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;

    ui->doubleSpinBox_23  ->setValue((double)numCandidates );

}

void MainWindow::on_pushButton_9_clicked()
{
 /*  int nvars = 5;
    double *f;
    int nrows = 1, nrows_eq;
    double *a, *b, *lb,*ub, *x,*a_eq, *b_eq, fval;
 int rez =    linprog11(  nvars,  f, nrows, a,  b
                 ,  nrows_eq,  a_eq,  b_eq
                  ,lb, ub
                , x, fval);*/


    // параметры лемнитскаты (овала Кассини):

    mLemn_a = ui->doubleSpinBox_a->value();
    mLemn_c = ui->doubleSpinBox_c->value();
    mLemn_l = ui->doubleSpin_L->value();
    ///
    /// \brief wchOutPutFold0

    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mpwchPeredChislaSCVFIle);

   // wchar_t *pch = wcsrchr(wchOutPutFold0, L'/');
   // wchar_t *pch = wcsrchr(wchOutPutFold0, L'\\');
     wchar_t *pch = wcsrchr(wchOutPutFold0, L'/');

    pch[0] = 0;
    wcscat(wchOutPutFold0, L"\\LEMNISCATA\\");

     _wmkdir(wchOutPutFold0);

    ui->doubleSpinBoxZ1  ->setValue(-mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l);
    ui->doubleSpinBoxZ2  ->setValue(mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l);
    ui->doubleSpinBoxYmax->setValue(sqrt(-mLemn_c * mLemn_c + mLemn_a * mLemn_a));
    TURPolygon PlgCassini0 = TURPolygon::createCassiniOval((double)mLemn_a, (double)mLemn_c ,500 );
    // сдвиг полигона на вектор  pntSdvig
    TURPointXY pntSdvig (-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l, 0.) ;
    TURPolygon  PlgCassini1 =  PlgCassini0.SdvigTransform( pntSdvig );


    // оси  координат
    wchar_t wchAxesFileName[300] ={0};
    wcscpy(  wchAxesFileName,  wchOutPutFold0);
    wcscat(wchAxesFileName, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-5000., 5000.
    ,-10000., 10000.,30.) ;
//
    wchar_t wchFileName[300] ={0};
    wcscpy(  wchFileName,  wchOutPutFold0);
    wcscat(wchFileName, L"\\Oval.shp");
    TURPolygon::WriteSetSHPFiles(wchFileName,&PlgCassini1, 1) ;
    ///

   long double arrA[]= {1.,2.,3.,4.};
    double z1 = -mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l;
    double z2 = mLemn_c-sqrt(mLemn_c * mLemn_c + mLemn_a * mLemn_a) - mLemn_l;
    double val_A = mLemn_a * mLemn_a;
}


void MainWindow::on_pushButton_10_clicked()
{
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с ПЧ малой1 задачи","D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY" ,"*.csv");
    this->ui->lineEdit_8->setText(strFold);

    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
    int quantCols = TYrRead::calcColCountFromCSV(mpwchPickedOutCandidates ) ;

    int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;

    ui->doubleSpinBox_24  ->setValue((double)numCandidates );

}

void MainWindow::on_pushButton_11_clicked()
{
   mqstrPeredChislaSCVFIle = QFileDialog::getSaveFileName(0,"Выбор .SCV файла с кандидатами  в ПЧ", "D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY","*.csv");
   this->ui->lineEdit_4->setText(mqstrPeredChislaSCVFIle);

     wchar_t array [400] = {0};
     mqstrPeredChislaSCVFIle.toWCharArray(array);
     array[mqstrPeredChislaSCVFIle.length()] = 0;
    wcscpy(mpwchPeredChislaSCVFIle, array);


}

void MainWindow::on_pushButton_8_clicked()
{


    inputData();
    mReMin = ui->doubleSpinBox_25->value();

    mReMax = ui->doubleSpinBox_27->value();

    mTgMax = ui->doubleSpinBox_26->value();

    // директия для графиков собственных чисел

    wchar_t wchOutPutFold0[400] = {0};
    wcscpy(wchOutPutFold0, mpwchPeredChislaSCVFIle);



   /*  int ii = wcslen(wchOutPutFold0);
    for(int j = (wcslen(wchOutPutFold0) -1); j > 0; --j)
    {
         wchar_t wch = wchOutPutFold0[j];
       if(wchOutPutFold0[j] == L'/')
       {
         wchOutPutFold0[j] = 0;

         break;
       }
    }
*/
    wchar_t *pwchr = wcsrchr(wchOutPutFold0, L'/');
    pwchr[0] = 0;
  // String wchFoldName = mwchOutFold;

   // wchar_t *pch = wcsrchr(wchOutPutFold0, L'\\');
   // pch[0] = 0;
   // String wchFileName = wchFoldName ;
          //  wchFileName += L"\\Image";
          //  wchFileName += i;
           // wchFileName += L".shp";
    wcscat(wchOutPutFold0, L"\\SobstvChisla\\");

     _wmkdir(wchOutPutFold0);

     ///
    const double VAlVx = ui->doubleSpinBox_19->value();
    const double VAlY = ui->doubleSpinBox_20->value();

    // массив с вариантами наборов передаточных чисел
   if(mparrCandidates)
   {
       free (mparrCandidates);
   }

    int quanTitleRows = 0;
    int quantCols = 10;


    mnumCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;




    // к-во переменных в задаче. В большой 11, в малой 5
    int iNumVars = 5;
    ///
    ///
    // к-во управлений в задаче. В большой 4, в малой 2
    int iNumRuls = 2;
    ///




    mparrCandidates = (double*)malloc(mnumCandidates * quantCols * sizeof(double));
    int nrows = 0;
      TYrRead::YrReadTabCSV_1(mpwchPickedOutCandidates// файл с таблицей
          ,0// к-во заголов строк
          ,0 // к-во заголов cтолбцов
          ,mnumCandidates // к-во строк содержательной части таблицы
          ,quantCols// к-во столбцов содержательной части таблицы
          ,&nrows  // к-во прочитанных строк
          ,mparrCandidates  // массив в который принимается информация
          )  ;
    /*  double *parrBuffer = new double [mnumCandidates * 10];
      memset(parrBuffer , 0., mnumCandidates * 10 * sizeof(double));
      for (int i =0; i < mnumCandidates; ++i)
      {
        parrBuffer[i * 10] = mparrCandidates[ 5 * i] ;
        parrBuffer[i * 10 +3] = mparrCandidates[ 5 * i+1] ;
        parrBuffer[i * 10 +4] = mparrCandidates[ 5 * i+2] ;
        parrBuffer[i * 10 +6] = mparrCandidates[ 5 * i+3] ;
        parrBuffer[i * 10 +7] = mparrCandidates[ 5 * i+4] ;
      }
      TYrWrite::WriteMassiveInFIleSCV(mpwchPickedOutCandidates,parrBuffer, mnumCandidates, 10
                                   ,NULL,NULL, 30);
      delete [] parrBuffer;
      return;*/
      ///







 long double arr_dF_po_dW [10] = {0.};
 long double arr_dF_po_dx [25] ={0};
 long double arr_C [10]  = {0.};

 long double arrT0 [25] ={0};
 long double arrT1 [25] ={0};



  mpPartHelicTraj->fill_df_po_px_and_df_po_dW_TangCanale(arr_dF_po_dx, arr_dF_po_dW);



 TURPointXY *pntarr = new TURPointXY[5];
 int maxQuant = 1000;
 int quantRows = 0;

 double *arrDataBuf = (double *)malloc(maxQuant * quantCols * sizeof(double));
 TComp  cmparrRoots[5];

 for (int i =0; i < mnumCandidates; ++i)

 {

     memset(arr_C, 0, iNumRuls * iNumVars *sizeof(long double ));
    for (int j = 0; j <quantCols; ++j)
    {
       arr_C[j] = mparrCandidates [i * quantCols + j];
    }
     //arr_C[0] = mparrCandidates [i * quantCols];
    // arr_C[3] = mparrCandidates [i * quantCols + 1];
    // arr_C[4] = mparrCandidates [i * quantCols + 2];
    // arr_C[6] = mparrCandidates [i * quantCols + 3];
    // arr_C[7] = mparrCandidates [i * quantCols + 4];

     MtrxMultMatrx(arr_dF_po_dW, iNumVars, iNumRuls, arr_C, iNumVars, arrT0) ;

     MtrxSumMatrx(arr_dF_po_dx,arrT0 ,iNumVars, iNumVars, arrT1) ;


     if((TPartHelicTraj::IsRootsSuit_dim5(arrT1,5, mReMin
                                        ,mReMax, mTgMax , cmparrRoots)))
     {

            memcpy(&(arrDataBuf [ quantCols * quantRows]),&(mparrCandidates [i * quantCols]), sizeof(double) *quantCols);

            quantRows++;
            if (quantRows== maxQuant)
            {
                break;
            }

     }
     else
     {
        // continue;
     }


  for (int j = 0; j < 5; ++j)
  {
       pntarr[j].X = cmparrRoots[j].m_Re;
       pntarr[j].Y = cmparrRoots[j].m_Im;
       if(pntarr[j].X >= 0.)
       {
           int qqq=1;
       }
   }
     wchar_t wchFileName[400] = {0};
     wcscpy(wchFileName, wchOutPutFold0);
     _wmkdir(wchFileName);
     wcscat(wchFileName, L"\\SobstvChisla.shp");

     TURPointXY::WriteSetSHPFiles(wchFileName,pntarr, 5) ;
     int iii0 = 0;
  }






 // запись массива информации в  CSV файл
 // INPUT:
 // FileName
 // parrBuff[ iNumRows * iNumCols] - массив с информацией
 // iNumRows- к-во строк массива
 // iNumCols - к-во столбцов массива
 // pwcharrRowNames[ iLenName* iNumRows] - имена строк массива
 // pwcharrColNames [ iLenName* iNumCols] - имена сьолбцов массива
 TYrWrite::WriteMassiveInFIleSCV(mpwchPeredChislaSCVFIle,arrDataBuf, quantRows, quantCols
                              ,NULL,NULL, 30);

 free(arrDataBuf) ;
// free(arr_dF_po_dW);
 //free(arr_dF_po_dx);
// free(arr_C);

 //free(arrT0);
 //free(arrT1);




 free(mparrCandidates);
 mparrCandidates = NULL;

 QString  qstr= QString::number(quantRows,10);
 ui->lineEdit_9->setText(qstr);




 // оси  координат
 wchar_t wchAxesFileName0[300] ={0};
 wcscpy(  wchAxesFileName0,  wchOutPutFold0);
 wcscat(wchAxesFileName0, L"\\AxesArr.shp");
 TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
 ,-10000., 10000.,30.) ;
 delete [] pntarr;

}

void MainWindow::on_comboBox_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:  // равном прямол
        ui->doubleSpinBox_19->setVisible(true);
        ui->doubleSpinBox_20->setVisible(true);
        ui->doubleSpinBox->setVisible(false);
        ui->doubleSpinBox_2->setVisible(false);
        ui->doubleSpinBox_3->setVisible(false);
        //ui->doubleSpinBox_5->setVisible(false);
        break;
    case 1:    // маневр
        ui->doubleSpinBox_19->setVisible(true);
        ui->doubleSpinBox_20->setVisible(true);
        ui->doubleSpinBox->setVisible(true);
        ui->doubleSpinBox_2->setVisible(true);
        ui->doubleSpinBox_3->setVisible(false);
       // ui->doubleSpinBox_5->setVisible(false);
        break;
    case 2:     // висение
        ui->doubleSpinBox_19->setVisible(false);
        ui->doubleSpinBox_20->setVisible(true);
        ui->doubleSpinBox->setVisible(false);
        ui->doubleSpinBox_2->setVisible(false);
        ui->doubleSpinBox_3->setVisible(false);
       // ui->doubleSpinBox_5->setVisible(true);
        break;
    case 3:      // вращение
        ui->doubleSpinBox_19->setVisible(false);
        ui->doubleSpinBox_20->setVisible(true);
        ui->doubleSpinBox->setVisible(false);
        ui->doubleSpinBox_2->setVisible(false);
        ui->doubleSpinBox_3->setVisible(true);
      //////////  ui->doubleSpinBox_5->setVisible(false);
    default: break;
    }
}



void MainWindow::on_progressBar_valueChanged(int value)
{

}

void MainWindow::on_comboBox_4_currentIndexChanged(int index)
{
    switch(index)
    {
    case 0:
      mbTang = true;
        break;
    case 1:
        mbTang = false;
        break;
    default: break;

    }

}

void MainWindow::on_doubleSpinBox_c_valueChanged(double arg1)
{
   mLemn_c = arg1;
   ui->doubleSpinBox_a->setValue(mLemn_c * sqrt(2.));
   mLemn_a = mLemn_c * sqrt(2.);
}
