//в настоящем пректе производится сравнение летных характеристик векртолета, приведенных в РЛЭ
// и модельного вертолета
// сравнение производится путем наложения графиков
//



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
#include "Plane.h"
#include "URPolyLine.h"
#include "UrPointXY.h"


#include "Blade.h"
#include "Environment.h"
#include "YrWriteShapeFile.h"
//#include "URPolygon.h"
//#include "UrPointXY.h"
#include "dir.h"
#include "Comp.h"
#include "Equations.h"
#include "Helic.h"
#include "MatrixProccess.h"
#include "YrWrite.h"
#include "HelicConstants.h"
#include "BallanceCalc.h"
#include "PartHelicTraj.h"
#include "LineMove.h"
#include "TurnMove.h"
#include "Hover.h"
#include "BallanceCalc.h"



#define NOT_POSSIBLE_VALUE -1000000000.



extern  long double constArrGearSets[];

extern TURPolyLine Pln_XB_XT125_p7_41;
extern TURPolyLine Pln_XB_XT40_p7_41;
extern TURPolyLine Pln_XB_XT80_p7_41;
extern TURPolyLine Pln_Nu_XT40_p7_41;
extern TURPolyLine Pln_Nu_XT80_p7_41;
extern TURPolyLine Pln_Nu_XT125_p7_41;
extern TURPolyLine Pln_Nu_XT80_p7_44;
extern TURPolyLine Pln_Xk_XT80_p7_43;
extern TURPolyLine Pln_XH_XT80_p7_43;
extern TURPolyLine Pln_Gamma_p7_45;

extern  long double constArrHelicPlanersData[8 * 13];
extern  long double ARrINertMtrxPrived[9];
extern const  long double constAlfaZaklNew;
extern   const  long double constXdistLow ;
extern   const  long double constDeltaXUp_Low;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->label_17->setVisible(true);
    ui->doubleSpinBox_13->setVisible(true);
            ui->doubleSpinBox_12->setVisible(true);
            ui->doubleSpinBox_4->setVisible(true);
            ui->doubleSpinBox_6->setVisible(true);

// установка грида
    ui->gridLayout->setSpacing(3);
    ui->gridLayout->addWidget(ui->label_3,0,0);
    ui->gridLayout->addWidget(ui->label_4,0,1);
    ui->gridLayout->addWidget(ui->label_5,0,2);
    ui->gridLayout->addWidget(ui->label_6,0,3);
    ui->gridLayout->addWidget(ui->doubleSpinBoxVHor,1,0);
    ui->gridLayout->addWidget(ui->doubleSpinBoxCourseAng,1,1);
    ui->gridLayout->addWidget(ui->doubleSpinBoxVVert,1,2);
    ui->gridLayout->addWidget(ui->doubleSpinBoxTemperature0,1,3);
///

    // полетное время
    mTFly = 1000.;
    ui->doubleSpinBox_3->setValue(mTFly);


    ///



  // атмосфера
   mWindHor = 0.; // горизонтьальная скорость ветра
   mWindCourse = 0.;// курсовой угол вектора горизонтьальной скорости ветра
   mWindVert = 0.;// вертикальная скорость ветра полож направление вверх
   mTemperature0 = 18.;
   ui->doubleSpinBoxTemperature0->setValue(mTemperature0);

   ui->doubleSpinBoxVVert->setValue(mWindVert);

   ui->doubleSpinBoxCourseAng->setValue(mWindCourse);

   ui->doubleSpinBoxVHor->setValue(mWindHor);

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

    mvalYaw =  (ui->doubleSpinBox_6->value()) * M_PI/ 180.;
    mvalRad =  ui->doubleSpinBox_4->value();
    mvalVx = ui->doubleSpinBox_12->value();
    mvalY = ui->doubleSpinBox_13->value();
    m_dPsi_po_dt = ui->doubleSpinBox->value();
    //mvalPsi0 = ui->doubleSpinBox_15->value();
    double valPsiBegin = ui->doubleSpinBox_20->value();
    mval_dPsi0 = ui->doubleSpinBox_2->value();

    QString strFold = this->ui->lineEdit->text();

    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

    mbFullGliders = false;
    mbFUllRotorModel = false;
    if (this->ui->comboBox_2->currentIndex() == 0)
    {
      mbFullGliders = false;
    }
    else
    {
       mbFullGliders = true;;
    }

    if (this->ui->comboBox_3->currentIndex() == 0)
    {
      mbFUllRotorModel = false;
    }
    else
    {
       mbFUllRotorModel = true;
    }
  // mHelic =  THelic(0);

 //   mHelic =  THelic(0);
    /// вертолет создан

    // создание окружающей среды (атмосферы)
    mEnvironment  = TEnvironment(mWindHor, mWindCourse, mWindVert, mTemperature0 );



}







void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY\\BigTask");
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

void MainWindow::on_tableWidget_2_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn)
{

}





void MainWindow::on_doubleSpinBoxVHor_valueChanged(const QString &arg1)
{

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






    //----------------------------------------------------------------------------------
// построение графиков балансировочных кривых со стр 7-41 РЛЭ
// строится зависимость от скорости.
// диапазон изменения скорости от 50 км/ч до 280 км/ч с шагом 10 км/ч
// всего точек на графике 24

const int QUantRows = 281 ;
int iNumCols = 9;
double *arrBuff = (double *)malloc(sizeof(double) * iNumCols *(QUantRows ));
memset(arrBuff, 0, sizeof(double) * iNumCols *(QUantRows ));


// графики строятся для 3 значений центровки 40мм, 80мм и 125 мм
long double arrCentrovka[3] = {-0.04,-0.08, -0.125};

// угол сколджения
const long double VAlPsi = 0.;
///

// скорость
long double valVx =0.;
///
long double VAlMass = 10800.;
long double arrCoordCentreMass[3] = {0.};
// решение ур-я балансировки
long double arrXRez[6] = {0.};
///

for(int i =0; i < 3; i++)
{


    THelic HelicCur( mbFUllRotorModel,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                                   ,  arrCoordCentreMass,constAlfaZaklNew, constXdistLow+ constDeltaXUp_Low
                                   , constXdistLow ,arrCentrovka[i]);

    for (int j = 0; j < QUantRows; ++j)
    {
       valVx =  ((long double)j)* 1000./ 3600.;
       long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
       bool brez = TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                                                    ,mvalY, valVx, VAlPsi
                                                                 ,arrX0, arrXRez );
       arrBuff[j * iNumCols] =  ((long double)j);
       arrBuff[j * iNumCols +1 + 2 * i    ] = arrXRez[2];
       arrBuff[j * iNumCols +1 + 2 * i + 1] = arrXRez[1];
       if (i==1)
       {
         arrBuff[j * iNumCols +7 ] = arrXRez[4] *100./HelicCur.mKappaTettaMax;
         arrBuff[j * iNumCols +8 ] = arrXRez[5]*50./HelicCur.mKappaTettaMax;
       }
    }




}



    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[iNumCols * iLenName];
    memset(pwcharrColNames, 0, iNumCols * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"V");
    wcscpy(&pwcharrColNames[iLenName   ],  L"Xв_XT=40");
    wcscpy(&pwcharrColNames[iLenName *2],  L"Nu_XT=40");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Xв_XT=80");
    wcscpy(&pwcharrColNames[iLenName *4],  L"Nu_XT=80");
    wcscpy(&pwcharrColNames[iLenName*5 ],  L"Xв_XT=125");
    wcscpy(&pwcharrColNames[iLenName *6],  L"Nu_XT=125");
    wcscpy(&pwcharrColNames[iLenName *7],  L"XK_XT=80");
    wcscpy(&pwcharrColNames[iLenName *8],  L"XH_XT=80");
    // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.) ;
    ///


    // график
    double arrScale[100];
    arrScale[0] = 1.;
    for (int i = 1; i < 7;++i)
    {
     arrScale[i] = 1000.;
    }
    arrScale[7] = 1.;
    arrScale[8] = 1.;


  for (int i =1; i < iNumCols; i++)
  {
    TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                    ,QUantRows //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }



    free(arrBuff);
    delete[] pwcharrColNames;

    //1
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_XB_XT80_p7_41.shp");
    Pln_XB_XT80_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_XB_XT80_p7_41, 1) ;

    //2
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_XB_XT40_p7_41.shp");
    Pln_XB_XT40_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_XB_XT40_p7_41, 1) ;

    //3
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_XB_XT125_p7_41.shp");
    Pln_XB_XT125_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_XB_XT125_p7_41, 1) ;

    // 4
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Nu_XT40_p7_41.shp");
    TURPolyLine Pln_Nu_XT40_p7_41_0 =  Pln_Nu_XT40_p7_41.MultScalar(1000. * M_PI/180. );
    Pln_Nu_XT40_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_Nu_XT40_p7_41_0, 1) ;

    //5
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Nu_XT80_p7_41.shp");
    TURPolyLine Pln_Nu_XT80_p7_41_0 =  Pln_Nu_XT80_p7_41.MultScalar(1000. * M_PI/180. );
    Pln_Nu_XT80_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_Nu_XT80_p7_41_0, 1) ;

    //6
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Nu_XT125_p7_41.shp");
    TURPolyLine Pln_Nu_XT125_p7_41_0 =  Pln_Nu_XT125_p7_41.MultScalar(1000. * M_PI/180. );
    Pln_Nu_XT125_p7_41.WriteSetSHPFiles(wchAxesFileName0,&Pln_Nu_XT125_p7_41_0, 1) ;

    //7
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Xk_XT80_p7_43.shp");
    TURPolyLine Pln_Xk_XT80_p7_43_0 =  Pln_Xk_XT80_p7_43;
    Pln_Xk_XT80_p7_43_0.WriteSetSHPFiles(wchAxesFileName0,&Pln_Xk_XT80_p7_43_0, 1) ;

    //8
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_XH_XT80_p7_43.shp");
    TURPolyLine Pln_XH_XT80_p7_43_0 =  Pln_XH_XT80_p7_43;
    Pln_XH_XT80_p7_43_0.WriteSetSHPFiles(wchAxesFileName0,&Pln_XH_XT80_p7_43_0, 1) ;

    //9
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Nu_XT80_p7_44.shp");
    TURPolyLine Pln_Nu_XT80_p7_44_0 =  Pln_Nu_XT80_p7_44.MultScalar(1000. * M_PI/180. );
    Pln_Nu_XT80_p7_44.WriteSetSHPFiles(wchAxesFileName0,&Pln_Nu_XT80_p7_44_0, 1) ;
  //  extern TURPolyLine Pln_Nu_XT80_p7_44







    //----------------------------------------------------------------------------------
// построение графиков балансировочных кривых со стр 7-45 РЛЭ
// строится зависимость гамма от скорости бокового ветра.
// диапазон изменения скорости от -30 м/с до 30 м/с с шагом 1 м/с
// всего точек на графике 24

const int QUantRows1 = 61 ;
int iNumCols1 = 4;
double *arrBuff1 = (double *)malloc(sizeof(double) * iNumCols1 *(QUantRows1 ));
memset(arrBuff1, 0, sizeof(double) * iNumCols1 *(QUantRows1 ));




// угол сколджения
const long double VAlPsi1 = 0.;
///

// скорость
long double valVx1 =0.;
///


///


    THelic HelicCur( mbFUllRotorModel,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                                   ,  arrCoordCentreMass,constAlfaZaklNew, constXdistLow+ constDeltaXUp_Low
                                   , constXdistLow ,0.);

    for (int j = 0; j < QUantRows1; ++j)
    {
       TEnvironment EnvironmentCur(-30. + ((double)j) , -M_PI/2. , 0.);
       long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
       bool brez = TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,EnvironmentCur
                                 ,mvalY, valVx, VAlPsi,arrX0, arrXRez );
       arrBuff1[j * iNumCols1] =  -30. + ((double)j);
       arrBuff1[j * iNumCols1 +1] = arrXRez[0] * 180. / M_PI;
       arrBuff1[j * iNumCols1 +2] = arrXRez[4];
       arrBuff1[j * iNumCols1 +3] = arrXRez[5];

    }








    // вывод графиков

    wchar_t *pwcharrColNames1 = new wchar_t[iNumCols1 * iLenName];
    memset(pwcharrColNames1, 0, iNumCols1 * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames1, L"WindV");
    wcscpy(&pwcharrColNames1[iLenName   ],  L"Gamma");
    wcscpy(&pwcharrColNames1[iLenName *2],  L"Xk");
    wcscpy(&pwcharrColNames1[iLenName *3],  L"Xв");
    // оси  координат



    // график
    double arrScale1[100];

    for (int i = 0; i < iNumCols1;++i)
    {
     arrScale[i] = 1.;
    }



  for (int i =1; i < iNumCols1; i++)
  {
    TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                    ,arrBuff1 // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,iNumCols1  // - к-во переменных о корорых накоплена информация в буфере
                                    ,QUantRows1 //  - к-во точек
                                    ,pwcharrColNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }



    free(arrBuff1);
    delete[] pwcharrColNames1;

    //1
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\Pln_Gamma_p7_45.shp");
    Pln_Gamma_p7_45.WriteSetSHPFiles(wchAxesFileName0,&Pln_Gamma_p7_45, 1) ;


}

void MainWindow::on_pushButton_10_clicked()
{
    mqstrCandidatesSCVFIle = QFileDialog::getOpenFileName(0,"Выбор .SCV файла с комплектами перед. чисел болшой задачи", "D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY\\PeredChislaBolshZadachi");
   this->ui->lineEdit_2->setText(mqstrCandidatesSCVFIle);

     wchar_t array [400] = {0};
     mqstrCandidatesSCVFIle.toWCharArray(array);
     array[mqstrCandidatesSCVFIle.length()] = 0;
    wcscpy(mpwchPickedOutCandidates, array);
   int numCandidates  = TYrRead::YrCalcRows(mpwchPickedOutCandidates ) ;

    ui->doubleSpinBox_18  ->setValue((double)numCandidates );



}




void MainWindow::on_pushButton_3_pressed()
{

}

void MainWindow::on_comboBox_currentIndexChanged(int index)
{
   switch(index)
   {
   case 0:  // равном прямол
       ui->doubleSpinBox_13->setVisible(true);
       ui->doubleSpinBox_12->setVisible(true);
       ui->doubleSpinBox_4->setVisible(false);
       ui->doubleSpinBox_6->setVisible(false);
       //ui->doubleSpinBox_15->setVisible(false);
       ui->doubleSpinBox_2->setVisible(false);

       break;
   case 1:    // маневр
       ui->doubleSpinBox_13->setVisible(true);
       ui->doubleSpinBox_12->setVisible(true);
       ui->doubleSpinBox_4->setVisible(true);
       ui->doubleSpinBox_6->setVisible(true);
       //ui->doubleSpinBox_15->setVisible(false);
       ui->doubleSpinBox_2->setVisible(false);


       break;
   case 2:     // висение
       ui->doubleSpinBox_13->setVisible(true);
       ui->doubleSpinBox_12->setVisible(false);
       ui->doubleSpinBox_4->setVisible(false);
       ui->doubleSpinBox_6->setVisible(false);
     //  ui->doubleSpinBox_15->setVisible(true);
       ui->doubleSpinBox_2->setVisible(false);

       break;
   case 3:      // вращение
       ui->doubleSpinBox_13->setVisible(true);
       ui->doubleSpinBox_12->setVisible(false);
       ui->doubleSpinBox_4->setVisible(false);
       ui->doubleSpinBox_6->setVisible(false);
      // ui->doubleSpinBox_15->setVisible(false);
       ui->doubleSpinBox_2->setVisible(true);


   default: break;
   }

}

void MainWindow::on_comboBox_activated(const QString &arg1)
{

}

void MainWindow::on_comboBox_editTextChanged(const QString &arg1)
{

}

void MainWindow::on_pushButton_clicked()
{
    // начало координат:
    TURPointXY pnt00(379.6, -189.3);

    TURPointXY pntSdvig(-pnt00.X, -pnt00.Y);

    TURPointXY pnt0_10(379.6, -63.5);

    TURPointXY pnt30_0(741.4, -189.3);
    ///

    //double kx = 30. / (pnt30_0.X - pnt00.X);
   // double ky = 10. / (pnt0_10.Y - pnt00.Y);

    /// \brief quantPoliline
    ///
    ///
    int quantPoliline = 1;
    TURPolyLine *pURPolyLine = (TURPolyLine *)malloc(quantPoliline * sizeof(TURPolyLine));
    TURPolyLine **ppurLine = &(pURPolyLine);
    TURPolyLine::ReadSHPFile(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\gamma_page_7-45.shp",ppurLine,  &quantPoliline) ;

    TURPolyLine pln0 = (*ppurLine)[0].SdvigTransform( pntSdvig );

    double  arrMtxPer [4] = {0.};
    arrMtxPer [0] = 30./(pnt30_0.X - pnt00.X);
    arrMtxPer [3] = 10./(pnt0_10.Y - pnt00.Y);

    TURPolyLine pln1 = pln0.fncLinTransform(arrMtxPer );

    pln1.WriteToASCII__(L"D:\\REPOSITORIES\\gamma_page_7-45.txt");
    pln1.WriteSetSHPFiles(L"D:\\REPOSITORIES\\RLE_GRAPHS\\gamma_page_7-45.shp",&pln1, 1);


   // Pln_XB_XT125_p7_41.WriteSetSHPFiles(L"D:\\REPOSITORIES\\RLE_GRAPHS\\XB_XT_125_p7_41__.shp",&Pln_XB_XT125_p7_41, 1);
   // Pln_XB_XT80_p7_41.WriteSetSHPFiles(L"D:\\REPOSITORIES\\RLE_GRAPHS\\XB_XT_80_p7_41__.shp",&Pln_XB_XT80_p7_41, 1);
   // Pln_XB_XT40_p7_41.WriteSetSHPFiles(L"D:\\REPOSITORIES\\RLE_GRAPHS\\XB_XT_40_p7_41__.shp",&Pln_XB_XT40_p7_41, 1);
}

void MainWindow::on_pushButton_4_clicked()
{
    // на каринке 7.28 изоюбражен грацфик измения угла тангажа в зависимости от вкотости с центровкой 80 мм
    // из этой картинки видно, что при нулевой скорости угол тангажа равен 5 град
    // зафиксируем координаты основания НВ в сиситеме координат, в которой были получены параметры элементов планера.
    // предполагалось что это СвСК. но это оказалось не так, т к обнаружились большие расхождения с картинкой 7.25
    // ЦМ находится не там, где предполагалось. попробуем его найти.
    // ценнтр масс лежит на прямойц пароходжящей через точку основания вала винта и лежит на прямой являющейся его продолжением
    // обозначим через тау расстояние от основания НВ до центра масс. Попробует найти тау.
    // Есть 3 неизвестных - T, Kappa и tay. Есть 3 уравнения - сумма сил=0, сумма моментов =0, угол тангажа Nu = 5 grad.
    // решать систему не будем - будем подбирать. Зафиксируем тау.
    //пернсчитаем координаты точек начал координат СвСК элементов планера.
    //учтем балансировку (80мм) решим задачу балансировки.
    //подберем тау так, чтобы угол тангажа был бы равен 5 град
    inputData();

    long double VAlMass = 10000.;
    long double arrCoordCentreMass[3] = {0.};





    const long double valDist = 2.18;
    const long double valXT = 0.08;
    long double alfaZakl = TBallanceCalc::calcAlfaZakl(6./ 180. * M_PI, valXT, valDist);

   // THelic HelicCur(true, true, alfaZakl);

    long double arrX0[6]={0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};

    arrX0[1] = 6.* M_PI/ 180.;
    THelic HelicCur(true, true, VAlMass);
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                                          ,1000.,0.,0.
                                                          ,arrX0, arrXRez );



    HelicCur.changeCentreMass(-0.08,0.,0.);
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                                          ,1000.,0.,0.
                                                          ,arrX0, arrXRez );
    // точка верхнего винта
    TURPointXY pntRotorUp(HelicCur.marrRotor[0].mBasePLane.marrS0[0]* 1000., HelicCur.marrRotor[0].mBasePLane.marrS0[1]* 1000.) ;
    pntRotorUp.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\pntRotorUp.shp", &pntRotorUp,1);
    ///


    // точка нижнего винта
    TURPointXY pntRotorLow(HelicCur.marrRotor[1].mBasePLane.marrS0[0]* 1000., HelicCur.marrRotor[1].mBasePLane.marrS0[1]* 1000.) ;
    pntRotorLow.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\pntRotorLow.shp", &pntRotorLow,1);
    ///

    // построить тоску ЦМ
    long double arrT0[3] = {0.},arrT1[3] = {0.};
    MtrxSumMatrx(HelicCur.marrRotor[0].mBasePLane.marrS0, HelicCur.marrRotor[1].mBasePLane.marrS0,1, 3, arrT0) ;

    ///


    MatrxMultScalar(arrT0, 1, 3, 0.5,arrT1);
    TURPointXY pntRotorMid(arrT1[0]* 1000., arrT1[1]* 1000.) ;
    pntRotorMid.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\pntRotorMid.shp", &pntRotorMid,1);


    arrT1[0] -= valDist* sinl(HelicCur.marrRotor[0].mZaklinAng);
    arrT1[1] -= valDist * cosl(HelicCur.marrRotor[0].mZaklinAng);

    TURPointXY pntCentreMass0(arrT1[0]* 1000., arrT1[1]* 1000.) ;
    pntCentreMass0.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\pntCentreMass0.shp", &pntCentreMass0,1);

   /* long double tau0 = (HelicCur.marrRotor[0].mBasePLane.marrS0[1]- 1.22)/cos(HelicCur.marrRotor[0].mZaklinAng);

    long double nu0 = 6./180. * M_PI;

    long double taurez = 0.;
    long double temp = 0.;
    TURPointXY pntCentreMass;
    long double xrez;
    for (int i =0; i < 4000; ++i)
    {
      long double tau = -tau0 - ((long double)i) * 0.001;
      long double arrCM_BSK[3] = {0.}, arrCM_SvSK[3] = {0.};
      arrCM_BSK[1] = tau;
      HelicCur.marrRotor[0].mBasePLane.transform_xyzSKP_to_xyzSSK(arrCM_BSK,arrCM_SvSK);
      pntCentreMass.X = arrCM_SvSK[0]* 1000.;
      pntCentreMass.Y = arrCM_SvSK[1]* 1000.;
      arrCM_SvSK[0] += 0.08;
      THelic HelicCur1 (mbFullGliders,mbFUllRotorModel);
      long double parrRez[3] = {0.};
      MtrxMinusMatrx(HelicCur1.mRuleGlEl.mPlaneSvSK.marrS0, arrCM_SvSK,1, 3, parrRez);
      memcpy(HelicCur1.mRuleGlEl.mPlaneSvSK.marrS0, parrRez, 3 * sizeof(long double));
      for (int j = 0; j < 7; ++j)
      {
          MtrxMinusMatrx(HelicCur1.marrComGlEl[j].mPlaneSvSK.marrS0, arrCM_SvSK,1, 3, parrRez);
          memcpy(HelicCur1.marrComGlEl[j].mPlaneSvSK.marrS0, parrRez, 3 * sizeof(long double));

      }

      MtrxMinusMatrx(HelicCur1.marrRotor[0].mBasePLane.marrS0, arrCM_SvSK,1, 3, parrRez);
      memcpy(HelicCur1.marrRotor[0].mBasePLane.marrS0, parrRez, 3 * sizeof(long double));


      MtrxMinusMatrx(HelicCur1.marrRotor[1].mBasePLane.marrS0, arrCM_SvSK,1, 3, parrRez);
      memcpy(HelicCur1.marrRotor[1].mBasePLane.marrS0, parrRez, 3 * sizeof(long double));

     // HelicCur1.marrRotor[0].mBasePLane.marrS0[0] += 0.1;
      //HelicCur1.marrRotor[1].mBasePLane.marrS0[0] += 0.1;

      long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.}, arrXRez[6] ={0.};
      long double valVx = 0., VAlPsi = 0.;
      bool brez = TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur1,mEnvironment
                                                                   ,mvalY, valVx, VAlPsi
                                                                ,arrX0, arrXRez );

      long double temp0 = arrXRez[1] - nu0;
      xrez = arrXRez[1];
      if ( i >0)
      {
          if (temp * temp0 <=0.)
          {
             taurez = tau;
             break;
          }
      }
      temp = temp0;


    }
pntCentreMass.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\RLE\\pntCentreMass0.shp", &pntCentreMass,1);
 int iii=0;
*/

}

void MainWindow::on_pushButton_5_clicked()
{
    inputData();
    //const long double valDist = 2.18;
    const long double valXT = 0.08;
    long double VAlAlfaZakl = 3.3 * M_PI/180.;//3.3 * M_PI/180.;
    const long double VAlMass = 10800.;
    long double arrCoordCentreMass [3] = {0.};
/*
    long double arrxdist0[2] = {1.5, 2.8}, arr_dx[2] = {0.001, 0.001};;

    long double arrRight[2] = {6.* M_PI/ 180., 5.99 * M_PI /180.};
    long double arrF[2] = {0.}, arrJacF[4] = {0.},arrJacFInv[4] = {0.}, arrDelX[2] ={0.}, arrRez[2] = {0.};

    for (int i =0; i < 40; ++i)
    {
        THelic HelicCur( true,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                         ,arrCoordCentreMass,VAlAlfaZakl, arrxdist0[0], arrxdist0[1], -valXT);

        calcVectF_and_JacF(VAlAlfaZakl, arrxdist0,  valXT,arrRight,arrF, arrJacF);
        InverseMtrx2(arrJacF, arrJacFInv);
        MtrxMultMatrx(arrJacFInv,2, 2, arrF,1, arrDelX) ;
        long double del = NormVect2(arrDelX);
        MtrxMinusMatrx(arrxdist0, arrDelX,1, 1, arrRez);
        memcpy(arrxdist0, arrRez, 2 * sizeof(long double));
        if (del < 0.001)
        {
            break;
        }

    }
*/
   long double deltadist = 1.3, xdist0 = 1.5;
    long double arrX0[6]={0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};



    for (int i =0; i < 40; ++i)
    {
        THelic HelicCur( true,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                                       ,  arrCoordCentreMass,VAlAlfaZakl, xdist0 + deltadist
                                       , xdist0,-valXT);
        TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                               ,1000.,0.,0.,arrX0, arrXRez );
        long double valf = arrXRez[1] - 6.* M_PI/ 180.;
        long double temp =arrXRez[1];
        THelic HelicCur1( true,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                                       ,  arrCoordCentreMass,VAlAlfaZakl, xdist0 + 0.01 + deltadist
                                       , xdist0+ 0.01,-valXT);
        TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur1,mEnvironment
                                               ,1000.,0.,0.,arrX0, arrXRez );
        long double val_df = (arrXRez[1] - temp)/0.01;
        long double del = valf/ val_df;
        xdist0 -= del;
        if (fabsl(del) <= 0.001)
        {
                break;
        }

    }

    THelic HelicCur( true,  VAlMass, ARrINertMtrxPrived, constArrHelicPlanersData
                                   ,  arrCoordCentreMass,VAlAlfaZakl, xdist0 + deltadist
                                   , xdist0,-valXT);
    long double arrXRez1[6] = {0.};
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                           ,1000.,0.,0.,arrX0, arrXRez );
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                           ,1000.,10. * 1000./ 3600.,0.,arrX0, arrXRez1 );
    double temp = arrXRez[1]/ M_PI * 180.;

    int yyy=0;
}

void MainWindow::calcVectF_and_JacF(long double VAlAlfaZakl,long double* arrxdist0, long double valXT
                                    , long double *arrRight, long double *arrF, long double *arrJacF)
{
    long double arrCoordCentreMass[3] = {0.};
    THelic HelicCur( true,  10800., ARrINertMtrxPrived, constArrHelicPlanersData
                     ,arrCoordCentreMass,VAlAlfaZakl, arrxdist0[1], arrxdist0[0], -valXT);

    // 1. вычисление arrF
    long double arrX0[6]={0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez_V0[6] = {0.}, arrXRez_V5[6] = {0.};
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                           ,1000.,0.,0.,arrX0, arrXRez_V0 );
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur,mEnvironment
                                           ,1000.,5. * 1000./ 3600.,0.,arrX0, arrXRez_V5 );
    arrF[0] = arrXRez_V0[1] - arrRight[0];
    arrF[1] = arrXRez_V5[1] - arrRight[1];
    long double arrtemp0 [2] = {0.};
    arrtemp0[0] = arrXRez_V0[1];
    arrtemp0[1] = arrXRez_V5[1];


    long double arrJacFT[4] = {0.}, arrF1[2] = {0.}, arrTemp0[2] = {0.};

    THelic HelicCur1( true,  10800., ARrINertMtrxPrived, constArrHelicPlanersData
                     ,arrCoordCentreMass,VAlAlfaZakl, arrxdist0[1]+ 0.005, arrxdist0[0], -valXT);
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur1,mEnvironment
                                           ,1000.,0.,0.,arrX0, arrXRez_V0 );
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur1,mEnvironment
                                           ,1000.,5. * 1000./ 3600.,0.,arrX0, arrXRez_V5 );
    arrF1[0] = arrXRez_V0[1] ;
    arrF1[1] = arrXRez_V5[1] ;
    MtrxMinusMatrx(arrF1, arrtemp0, 2, 1, arrTemp0);
    MatrxDivideScalar(arrTemp0, 1, 2, 0.005,arrJacFT);





    THelic HelicCur2( true,  10800., ARrINertMtrxPrived, constArrHelicPlanersData
                     ,arrCoordCentreMass,VAlAlfaZakl, arrxdist0[0], arrxdist0[1] + 0.005, -valXT);
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur2,mEnvironment
                                           ,1000.,0.,0.,arrX0, arrXRez_V0 );
    TBallanceCalc::calcBallParamsForSteadyLineMoving(HelicCur2,mEnvironment
                                           ,1000.,5. * 1000./ 3600.,0.,arrX0, arrXRez_V5 );
    arrF1[0] = arrXRez_V0[1] ;
    arrF1[1] = arrXRez_V5[1] ;
    MtrxMinusMatrx(arrF1, arrtemp0, 2, 1, arrTemp0);
    MatrxDivideScalar(arrTemp0, 1, 2, 0.005,&(arrJacFT[2]));

    MatrTransp(arrJacFT, 2, 2, arrJacF);


}
