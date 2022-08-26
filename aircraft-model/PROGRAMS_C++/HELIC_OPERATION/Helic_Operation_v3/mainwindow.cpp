#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <float.h>
#include <wchar.h>
#include <math.h>
#include <QFileDialog>
#include <QMessageBox>
#include <stdlib.h>

#include "YrRead.h"
#include "Plane.h"


#include "Blade.h"
#include "Environment.h"
#include "YrWriteShapeFile.h"


#include "HelicConstants.h"
#include "BallanceCalc.h"




#define NOT_POSSIBLE_VALUE -1000000000.


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
// установка грида
    ui->gridLayout->setSpacing(10);
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
    mTFly = 2000.;
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

}

MainWindow::~MainWindow()
{
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

    mTFly = ui->doubleSpinBox_3->value();


     ///
     ///  создание окружающей среды (атмосферы)
     mEnvironment  = TEnvironment(mWindHor, mWindCourse, mWindVert, mTemperature0 );
     ///    

}



void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\STABILITY\\OutPutHelicOperation");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}




void MainWindow::on_pushButton_3_clicked()
{
    const long double VAl_OutPut_dt = 0.01; // шаг по времени в соответствии с которым будет происходить выдачап информации
    inputData();
  //  const double VAlStIntegr = 0.01; // шаг интегрирования модели

    const int QUantRows = static_cast<int>(mTFly /VAl_OutPut_dt) +1 ;
    int iNumCols = 10;
    double *arrBuff = (double *)malloc(sizeof(double) * iNumCols *(QUantRows ));

    if(!arrBuff)
    {
        return;
    }
    memset(arrBuff, 0, sizeof(double) * iNumCols *(QUantRows ));

    // заготовка массива для построения графиков выходной информации для врешней модели
   // const int QUantColsOutArr = 26;
   // double *arrOutSideInfo = (double *)malloc(sizeof(double) * QUantColsOutArr *(QUantRows ));
    // memset(arrOutSideInfo, 0, sizeof(double) * QUantColsOutArr *(QUantRows ));
    double *arrOutSideInfo = (double *)malloc(sizeof(double) * (1 + QUANT_COLS_OUT_VECT) *(QUantRows ));

    if(!arrOutSideInfo)
    {
        return;
    }
    memset(arrOutSideInfo, 0, sizeof(double) * (1 + QUANT_COLS_OUT_VECT) *(QUantRows ));
    ///
    int quantDoneSteps = -1;
    int iTypeOfFlyTask = 1;
    THelicTraj HelicTraj   (mEnvironment, iTypeOfFlyTask, false, false);

    HelicTraj.imitateMoving_and_collectData(mTFly, VAl_OutPut_dt,  arrBuff,iNumCols, &quantDoneSteps, arrOutSideInfo);



    // вывод графиков моих
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[iNumCols * iLenName];
    memset(pwcharrColNames, 0, iNumCols * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"t");
    wcscpy(&pwcharrColNames[iLenName   ],  L"X");
    wcscpy(&pwcharrColNames[iLenName *2],  L"Y");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Z");
    wcscpy(&pwcharrColNames[iLenName *4],  L"VX");
    wcscpy(&pwcharrColNames[iLenName *5],  L"VY");
    wcscpy(&pwcharrColNames[iLenName *6],  L"VZ");
    wcscpy(&pwcharrColNames[iLenName *7], L"Gamma_Kren");
    wcscpy(&pwcharrColNames[iLenName *8], L"Nu_Tangaj");
    wcscpy(&pwcharrColNames[iLenName *9], L"Psi_Riskan");





    // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-50000., 50000.
    ,-100000., 100000.,30.) ;
    // график
    double arrScale[100];
    for (int i = 0; i < 100;++i)
    {
     arrScale[i] = 1.;
    }
    arrScale[7] = 100.;
    arrScale[8] = 100.;
    arrScale[9] = 100.;

  for (int i =1; i < iNumCols; i++)
  {
    TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                    ,quantDoneSteps //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantDoneSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,1  //  номер переменной по оси X
                                  ,2  //  номер переменной по оси Y
                                  ,arrScale[1] //  масштаб по оси X
                                  ,arrScale[2]  // масштаб по оси Y
                                   ) ;

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantDoneSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,1  //  номер переменной по оси X
                                  ,3  //  номер переменной по оси Y
                                  ,arrScale[1] //  масштаб по оси X
                                  ,arrScale[3]  // масштаб по оси Y
                                   ) ;

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantDoneSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,3  //  номер переменной по оси X
                                  ,2  //  номер переменной по оси Y
                                  ,arrScale[3] //  масштаб по оси X
                                  ,arrScale[2]  // масштаб по оси Y
                                   ) ;
  /// конец вывода моих графиков
/*

  // к-во столбцов вектора внешней модели, используемое для плстроения графиков
 // в первом столбце ( с нулевым номером) хранится время, далее хранится QUANT_COLS_OUT_VECT
 // элементов вектора
 const int  QUantColsOutArr = 1 + QUANT_COLS_OUT_VECT;

  //Вывод графиков с информацией для внешней модели
  wchar_t *pwcharrColNames1 = new wchar_t[ QUantColsOutArr * iLenName];
  memset(pwcharrColNames1, 0,  QUantColsOutArr * iLenName*sizeof(wchar_t));
  // составляющие воздушной скорости в связанной системе координат; 3
   // составляющие путевой скорости в ГСК (восточная, северная, модуль); 3
  // приборная скорость;                                          1

  // бароинерциальная вертикальная скорость;                      1
  // геометрическая высота;                                       1
  // абсолютная барометрическая высота;                           1
  // относительная барометрическая высота;                         1
  // геодезическая высота;                                        1
  // бароинерциальная высота;                                     1
  // угол сноса;                                                  1
  // угол атаки;                                                  1
  // угол скольжения;                                              1
  // скорость и направление ветра;                                2
  // крен;                                                        1
  // тангаж;                                                      1
  // курс истинный;                                                1
  //магнитное склонение;                                          1
  // температура наружного воздуха;                               1
  // текущие координаты местоположения объекта в географической системе координат. 2
  // модуль воздушной скорости                         1
  // вектор путевой скорости в СвСК                    3
  // боковая пергрузка
  wcscpy(pwcharrColNames1, L"t");
  wcscpy(&pwcharrColNames1[iLenName   ],  L"OUT_VXa");
  wcscpy(&pwcharrColNames1[iLenName *2],  L"OUT_VYa");
  wcscpy(&pwcharrColNames1[iLenName *3],  L"OUT_VZa");
  wcscpy(&pwcharrColNames1[iLenName *4],  L"OUT_V_GSK_East");
  wcscpy(&pwcharrColNames1[iLenName *5],  L"OUT_V_GSK_North");
  wcscpy(&pwcharrColNames1[iLenName *6],  L"OUT_Mod_V");
  wcscpy(&pwcharrColNames1[iLenName *7], L"OUT_V_pribor");// приборная скорость;
  wcscpy(&pwcharrColNames1[iLenName *8], L"OUT_V_baroVert");  // бароинерциальная вертикальная скорость;
  wcscpy(&pwcharrColNames1[iLenName *9], L"OUT_Hgeom");// геометрическая высота;

  wcscpy(&pwcharrColNames1[iLenName *10], L"OUT_H_absBarom");// абсолютная барометрическая высота;

  wcscpy(&pwcharrColNames1[iLenName *11], L"OUT_H_otnBarom"); // относительная барометрическая высота;

  wcscpy(&pwcharrColNames1[iLenName *12], L"OUT_H_geodes"); // геодезическая высота;
  wcscpy(&pwcharrColNames1[iLenName *13], L"OUT_H_baroinert");// бароинерциальная высота;
  wcscpy(&pwcharrColNames1[iLenName *14], L"OUT_Ugol_Snosa"); // угол сноса;
  wcscpy(&pwcharrColNames1[iLenName *15], L"OUT_Ugol_Ataki");  // угол атаки;
  wcscpy(&pwcharrColNames1[iLenName *16], L"OUT_Ugol_Skolgenia"); // угол скольжения;
  wcscpy(&pwcharrColNames1[iLenName *17], L"OUT_V_Wind"); // скорость  ветра;
  wcscpy(&pwcharrColNames1[iLenName *18], L"OUT_Ugol_Wind"); //  и направление ветра;
  wcscpy(&pwcharrColNames1[iLenName *19], L"OUT_Ugol_Krena"); // крен;
  wcscpy(&pwcharrColNames1[iLenName *20], L"OUT_Ugol_Tangaga"); // тангаж;
  wcscpy(&pwcharrColNames1[iLenName *21], L"OUT_Ugol_Course"); // курс истинный;
  wcscpy(&pwcharrColNames1[iLenName *22], L"OUT_MagDev"); //магнитное склонение;
  wcscpy(&pwcharrColNames1[iLenName *23], L"OUT_TAir_Out"); // температура наружного воздуха;
  wcscpy(&pwcharrColNames1[iLenName *24], L"OUT_Long"); // долгота
  wcscpy(&pwcharrColNames1[iLenName *25], L"OUT_Lat"); // широта

   wcscpy(&pwcharrColNames1[iLenName *27], L"OUT_ModVa_X"); // модуль воздушной скорости
   wcscpy(&pwcharrColNames1[iLenName *27], L"OUT_VSv_X"); // прод сост пут скор  в СвСК
   wcscpy(&pwcharrColNames1[iLenName *28], L"OUT_VSv_Y");
   wcscpy(&pwcharrColNames1[iLenName *29], L"OUT_VSv_Z");
   wcscpy(&pwcharrColNames1[iLenName *30], L"OUT_Peregruzka");







  for (int i = 0; i < 100;++i)
  {
   arrScale[i] = 1.;
  }
 // arrScale[7] = 100.;
  //arrScale[8] = 100.;
  //arrScale[9] = 100.;

for (int i =1; i < QUantColsOutArr; i++)
{
  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrOutSideInfo // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,QUantColsOutArr  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantDoneSteps //  - к-во точек
                                  ,pwcharrColNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,0  //  номер переменной по оси X
                                  ,i  //  номер переменной по оси Y
                                  ,1  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                   ) ;
}




  ///*/
    free(arrBuff);
    free(arrOutSideInfo);
    delete[] pwcharrColNames;
  //  delete []pwcharrColNames1;
}

void MainWindow::on_pushButton_clicked()
{
    /*TEnvironment Environment (10, M_PI/2., 0., 18. );
    const bool bFullType = false,  bTypeofRotor = false;
    THelic Helic( bFullType, bTypeofRotor);

    long double VAlY = 300., VAlVx = 40., VAlPsi = -0.24, VAlRad = 400.;
    long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};
  //  TBallanceCalc::calcBallParamsForSteadyLineMoving(Helic,Environment
     //               ,VAlY, VAlVx, VAlPsi,arrX0, arrXRez );
    for (int i =0; i < 200; ++i)
    {
     long double psi = -1. + 0.01 * static_cast<long double>(i);
     bool brez = TBallanceCalc::calcBallParamsForTurn(Helic,Environment
                     ,VAlY, VAlVx, psi,VAlRad, arrX0, arrXRez );
     if(brez)
     {
         int i0 = 0;
     }
    }


    TPartHelicTraj PartHelicTraj( Helic, Environment);
    PartHelicTraj.marrPhaseVect[1] = VAlY;
    PartHelicTraj.marrPhaseVect[3] = VAlVx;
    PartHelicTraj.marrPhaseVect[11] = VAlPsi;
    const long double VAlFlyTime = 2.;
    const long double VAlStIntegr = 0.001;
    int quantSteps = -1;

    const int QUantRows = static_cast<int>(VAlFlyTime /VAlStIntegr) +1 ;
    int iNumCols = 17;
    double *arrBuff = (double *)malloc(sizeof(double) * iNumCols *(QUantRows ));
    memset(arrBuff, 0, sizeof(double) * iNumCols *(QUantRows ));
    if(!PartHelicTraj.TryHorizontalTurn( VAlRad, VAlPsi, VAlFlyTime,  VAlStIntegr
                               , arrBuff, &quantSteps))
    {
        return;
    }















    // вывод графиков моих
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[iNumCols * iLenName];
    memset(pwcharrColNames, 0, iNumCols * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"t");
    wcscpy(&pwcharrColNames[iLenName   ],  L"X");
    wcscpy(&pwcharrColNames[iLenName *2],  L"Y");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Z");
    wcscpy(&pwcharrColNames[iLenName *4],  L"VX");
    wcscpy(&pwcharrColNames[iLenName *5],  L"VY");
    wcscpy(&pwcharrColNames[iLenName *6],  L"VZ");


    wcscpy(&pwcharrColNames[iLenName *7],  L"OmX");
    wcscpy(&pwcharrColNames[iLenName *8],  L"OmY");
    wcscpy(&pwcharrColNames[iLenName *9],  L"OmZ");

    wcscpy(&pwcharrColNames[iLenName *10], L"Gamma_Kren");
    wcscpy(&pwcharrColNames[iLenName *11], L"Nu_Tangaj");
    wcscpy(&pwcharrColNames[iLenName *12], L"Psi_Riskan");

    wcscpy(&pwcharrColNames[iLenName *13], L"Fi");
    wcscpy(&pwcharrColNames[iLenName *14], L"Kappa");
    wcscpy(&pwcharrColNames[iLenName *15], L"Etta");
    wcscpy(&pwcharrColNames[iLenName *16], L"DelFi");





    // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,  mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-50000., 50000.
    ,-100000., 100000.,30.) ;
    // график
    double arrScale[100];
    for (int i = 0; i < 100;++i)
    {
     arrScale[i] = 1.;
    }
    arrScale[7] = 100.;
    arrScale[8] = 100.;
    arrScale[9] = 100.;

  for (int i =1; i < iNumCols; i++)
  {
    TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                    ,quantSteps //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1  //  масштаб по оси X
                                  ,arrScale[i]  // масштаб по оси Y
                                     ) ;
  }

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,1  //  номер переменной по оси X
                                  ,2  //  номер переменной по оси Y
                                  ,arrScale[1] //  масштаб по оси X
                                  ,arrScale[2]  // масштаб по оси Y
                                   ) ;

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,1  //  номер переменной по оси X
                                  ,3  //  номер переменной по оси Y
                                  ,arrScale[1] //  масштаб по оси X
                                  ,arrScale[3]  // масштаб по оси Y
                                   ) ;

  TYrWriteShapeFile::WriteOneReport(mwchOutPutFold // путь к папке
                                  ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                  ,iNumCols  // - к-во переменных о корорых накоплена информация в буфере
                                  ,quantSteps //  - к-во точек
                                  ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                  ,iLenName // максимальная длина имени переменной
                                  ,3  //  номер переменной по оси X
                                  ,2  //  номер переменной по оси Y
                                  ,arrScale[3] //  масштаб по оси X
                                  ,arrScale[2]  // масштаб по оси Y
                                   ) ;




  ///
    free(arrBuff);

    delete[] pwcharrColNames;


  */
}
