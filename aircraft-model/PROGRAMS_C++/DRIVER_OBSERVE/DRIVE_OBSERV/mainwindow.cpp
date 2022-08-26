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






#define NOT_POSSIBLE_VALUE -1000000000.
extern const int QUantColsCSVReport;
extern const double CONST_ZP ;
extern const double GEARS[8];



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

   // ui->label_17->setVisible(true);
   // ui->doubleSpinBox_13->setVisible(true);
     //       ui->doubleSpinBox_12->setVisible(true);
          //  ui->doubleSpinBox_4->setVisible(true);
         //   ui->doubleSpinBox_6->setVisible(true);

// установка грида


    // полетное время
   // mTFly = 1000.;
  //  ui->doubleSpinBox_3->setValue(mTFly);


    ///



  // атмосфера



   mparrCandidates = NULL;



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







void MainWindow::inputData()
{
    int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) -2;
   if(mparrCandidates)
   {
       free (mparrCandidates);
   }
 mparrCandidates = (double*)malloc(numCandidates * 8 * sizeof(double));
 // int numRows = TYrRead::YrCalcRows(wchar_t *FileName)
 int nrows = 0;
   TYrRead::YrReadTabCSV_1(mpwchPickedOutCandidates// файл с таблицей
       ,2 // к-во заголов строк header
       ,0 // к-во заголов cтолбцов
       ,numCandidates // к-во строк содержательной части таблицы
       ,8// к-во столбцов содержательной части таблицы
       ,&nrows  // к-во прочитанных строк
       ,mparrCandidates  // массив в который принимается информация
       )  ;
   // заготовка матрицы передаточных чисел
   int numSet = int(ui->doubleSpinBox_17->value() + 0.1) -1;
   double arrC00[8] = {0.};
   for (int i =0; i < 8;++i)
   {

         arrC00[ i ] = mparrCandidates[numSet * 8 +  i ];


   }

   free(mparrCandidates);
   mparrCandidates = NULL;

    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;



    mInductL = (ui->doubleSpinBox_5->value())/1000.;
    mResist= ui->doubleSpinBox_11->value();
    mPsi_f= ui->doubleSpinBox_8->value();
    mJPayLoad  = ui->doubleSpinBox_10->value();
    mCx_om= ui->doubleSpinBox_7->value();
    mOmegaStat= ui->doubleSpinBox_25->value();
    mDriverModel = QElectDriver (mInductL , mResist ,
                                             mPsi_f , 170.
                                             , mJPayLoad, mCx_om,  CONST_ZP);

    mRealInductL = (ui->doubleSpinBox_15->value())/1000.;
    mRealResist= ui->doubleSpinBox_21->value();
    mRealPsi_f= ui->doubleSpinBox_19->value();
    mRealJPayLoad  = ui->doubleSpinBox_24->value();
    mRealCx_om= ui->doubleSpinBox_23->value();
    mRealOmegaStat= ui->doubleSpinBox_26->value();

    mDriverReal = QElectDriver (mRealInductL , mRealResist ,
                                             mRealPsi_f , 170.
                                             , mRealJPayLoad, mRealCx_om,  CONST_ZP);

    mDriveTrajReal = QDriveTraj(mDriverReal,arrC00);
    mTRotate = ui->doubleSpinBox_22->value();
}







void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\REPOSITORIES\\aircraft-model\\OUT_DRIVER");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}

//void MainWindow::on_pushButton_3_clicked()
//{
  //  mqstrInputFileName = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с ИД аэродинамических элементов", "D:\\REPOSITORIES\\aircraft-model\\СИСТЕМА_ДИФ_УРАВНЕНИЙ_ВЕРТОЛЕТА\\InputFile.csv","*.csv");
  // strFold =QString::fromWCharArray(L"D:\\REPOSITORIES\\aircraft-model\\СИСТЕМА_ДИФ_УРАВНЕНИЙ_ВЕРТОЛЕТА\\InputFile.csv",  -1);

   // this->ui->lineEdit_2->setText(mqstrInputFileName);





//}

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

    //wchar_t *pch = wcsrchr(wchOutPutFold0, L'/');
   // pch[0] = 0;
   //const double VAlT = 2.;
   const double VAlIntegrStep = 1. / 20000.;
   const int QUantRows = int(mTRotate/VAlIntegrStep) + 1;
   //const int iNumCols = 5;
   double *arrBuff = new double[QUantColsCSVReport * QUantRows];
   memset(arrBuff, 0, sizeof(double) * QUantColsCSVReport *(QUantRows ));
   int quantDoneSteps = -1;

   mDriveTrajReal.move( mTRotate,  VAlIntegrStep,mOmegaStat, arrBuff, &quantDoneSteps);

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

   arrScale[1] =arrScale[4] = 180./ M_PI;
   arrScale[2] = arrScale[3] = 1000.;

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
 double arrBuff1[13 *2] ={0.};
 arrBuff1[0] = 0.;
 arrBuff1[13] = mTRotate * 100.;
 double valCurrentIq = 0., valUd = 0., valUq = 0., valDelFi = 0.;
 mDriverModel.calcStationarySolution(mOmegaStat, &valCurrentIq, &valUd, &valUq, &valDelFi);
 arrBuff1[1] = arrBuff1[13 +1 ] =mOmegaStat;
 arrBuff1[3] = arrBuff1[13 +3 ] =valCurrentIq ;
 arrBuff1[5] = arrBuff1[13 +5 ] =valUd;
 arrBuff1[6] = arrBuff1[13 +6 ] =valUq;

  double valRealCurrentIq = 0., valRealUd = 0., valRealUq = 0., valRealDelFi = 0.;
 mDriverReal.calcStationarySolution(mOmegaStat, &valRealCurrentIq, &valRealUd, &valRealUq, &valRealDelFi);
 arrBuff1[6 + 1] = arrBuff1[13 + 6 +1 ] =mOmegaStat;
 arrBuff1[6 + 3] = arrBuff1[13 + 6 +3 ] =valRealCurrentIq ;
 arrBuff1[6 + 5] = arrBuff1[13 + 6 +5 ] =valRealUd;
 arrBuff1[6 + 6] = arrBuff1[13 + 6 +6 ] =valRealUq;

 int quantcols = 13;
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
 arrScale[1] =arrScale[4]=arrScale[7]=arrScale[10] = 180./ M_PI;
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

 //
 // вывод в окошки
 // расчетные с заложенной в мат модель
 ui->doubleSpinBox_6->setValue(0.);
 ui->doubleSpinBox_12->setValue(mOmegaStat / M_PI * 180.);
 ui->doubleSpinBox_9->setValue(0.);
 ui->doubleSpinBox_14->setValue(valCurrentIq * 1000000);
 ui->doubleSpinBox_13->setValue(valUd);
 ui->doubleSpinBox_28->setValue(valUq);

 // расчетные с реальной модели
 ui->doubleSpinBox_40->setValue(0.);
 ui->doubleSpinBox_35->setValue(mOmegaStat / M_PI * 180.);
 ui->doubleSpinBox_37->setValue(0.);
 ui->doubleSpinBox_38->setValue(valRealCurrentIq * 1000000);
 ui->doubleSpinBox_39->setValue(valRealUd);
 ui->doubleSpinBox_36->setValue(valRealUq);

 // с системы диф уравнений
 double *pp = &(arrBuff[(quantDoneSteps -1) *QUantColsCSVReport]);
 ui->doubleSpinBox_32->setValue(pp[4]);  // tetta
 ui->doubleSpinBox_29->setValue(pp[1] / M_PI * 180.);  // Om
 ui->doubleSpinBox_16->setValue(pp[2]);                  // Id
 ui->doubleSpinBox_30->setValue(pp[3] * 1000000.);// Iqu
 ui->doubleSpinBox_33->setValue(pp[5]);  //Ud
 ui->doubleSpinBox_31->setValue(pp[6]);  // Uqu


 delete[] pwcharrColNames;
 delete []pwcharrColNames1;
 delete [] arrBuff;


}

void MainWindow::on_pushButton_10_clicked()
{
    mqstrCandidatesSCVFIle = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с комплектами перед. чисел"
                         , "D:\\REPOSITORIES\\aircraft-model\\OUT_DRIVER");
   this->ui->lineEdit_2->setText(mqstrCandidatesSCVFIle);

     wchar_t array [400] = {0};
     mqstrCandidatesSCVFIle.toWCharArray(array);
     array[mqstrCandidatesSCVFIle.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
   int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) -2;

    ui->doubleSpinBox_18  ->setValue((double)numCandidates );




    TYrRead::ReadHdrFromGearFileForDriver(mpwchPickedOutCandidates
        , &mbVelo, &mInductL, &mResist,  &mPsi_f, &mJPayLoad
                              , &mCx_om , &mOmegaStat);
    if (mbVelo)
    {
    ui->lineEdit_3->setText("VELOCITY");
    }
    else
    {
        ui->lineEdit_3->setText("POSITION");
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
}




void MainWindow::on_pushButton_3_pressed()
{

}

void MainWindow::on_comboBox_currentIndexChanged(int index)
{


}

void MainWindow::on_comboBox_activated(const QString &arg1)
{

}

void MainWindow::on_comboBox_editTextChanged(const QString &arg1)
{

}
