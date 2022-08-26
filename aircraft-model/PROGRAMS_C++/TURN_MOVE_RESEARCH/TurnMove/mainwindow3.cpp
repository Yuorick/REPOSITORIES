#include "mainwindow3.h"
#include "ui_mainwindow.h"
#include <float.h>
#include <wchar.h>
#include <math.h>
#include <QFileDialog>
#include <QMessageBox>
#include <stdlib.h>

#include "YrRead.h"
#include "Plane.h"
//#include "Operation.h"
//#include "OperationNew.h"

#include "Blade.h"
#include "Environment.h"
#include "YrWriteShapeFile.h"


#include "HelicConstants.h"
#include "BallanceCalc.h"


extern const double constHelicMass;

#define NOT_POSSIBLE_VALUE -1000000000.


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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

    //mQuantActualTaskPoints = QUantActualTaskPoints0;
   // for (int i = 0 ; i < QUANT_FLY_TASK_POINTS * 3; i++)
   // {
   //  marrCoord[i] = NOT_POSSIBLE_VALUE;
   // }
   // for (int i = 0 ; i < QUANT_FLY_TASK_POINTS ; i++)
   // {
   //     marrV[i]      =   NOT_POSSIBLE_VALUE;
   //     marrCourse[i] =   NOT_POSSIBLE_VALUE;
   // }

   /* double arrCourseTemp[] = {0.,M_PI, -M_PI/2.,0., M_PI/2., -M_PI,-M_PI,-M_PI,-M_PI,-M_PI };

    double arrVTemp [] =     {0., 40.,   40.,   50., 50.,     50. , 40.,   40.,  10.,   0.};

    double arrCoordTemp[] = {    0.,   0.,    0.
                            ,-2500., 300.,    0.// 1.взлет с (0,0,0 ) на  скорость 40 на запад
                            ,-2900., 300., -3000.//2. поворот на 90 скорость 40 на север
                            , 4600., 300., -3400. //3. поворот на 90 скорость 50 на восток
                            , 5000., 300.,   -400.//4.поворот на 90 скорость 50 на юг
                            , 3000., 300.,     0. // 5.поворот на 90 скорость 50 на запад
                            , 1000., 300.,     0. //6. торможение до 40
                            , 500.,  100.,     0.  // 7.спуск на высоту 100 скорость 40 **
                            , 70.,  20.,     0. // 8.спуск на высоты 100м скорости 40  на высоту 20 и скорость 10
                            ,    0.,   0.,     0. // 9.посадка
                            };
    memcpy(marrCoord , arrCoordTemp,  mQuantActualTaskPoints * 3 * sizeof(double));
    memcpy(marrV     ,arrVTemp     ,  mQuantActualTaskPoints  * sizeof(double));
    memcpy(marrCourse,arrCourseTemp,  mQuantActualTaskPoints  * sizeof(double));
*/

   // memcpy(marrCoord , ARRCoord0,  mQuantActualTaskPoints * 3 * sizeof(double));
   // memcpy(marrV     ,ARRV0      ,  mQuantActualTaskPoints  * sizeof(double));
   // memcpy(marrCourse,ARRCourse0,  mQuantActualTaskPoints  * sizeof(double));



 // создание таблицы полетного задания
  //  ui->tableWidget->setColumnCount(5);
   // ui->tableWidget->setRowCount(QUANT_FLY_TASK_POINTS );


   // QHeaderView *pHeader = ui->tableWidget->horizontalHeader();
    //pHeader->setVisible(true);
   //pHeader = ui->tableWidget->verticalHeader();
   // pHeader->setVisible(true);



 /*   for (int i=0; i< ui->tableWidget->rowCount(); i++)
    for (int j = 0; j < ui->tableWidget->columnCount(); j++)
    {
    marrDblSpinBox[ i * 5 + j].setDecimals(4);
    double minimum = -DBL_MAX;
    double maximum = DBL_MAX;
    marrDblSpinBox[ i * 5 + j].setRange(minimum, maximum);
   // marrDblSpinBox[ i * 5 + j].setValue(marrCoord[i * 3 + j]);

    QAbstractSpinBox  *pabsSpin = &(marrDblSpinBox[ i * 5 + j]);
    (*pabsSpin).setButtonSymbols(QAbstractSpinBox::NoButtons) ;//= NoButtons;

  //  ui->tableWidget->setCellWidget(i,j,&(marrDblSpinBox[ i * 5 + j]));
    }

    for (int i=0; i< ui->tableWidget->rowCount(); i++)
    {
        marrDblSpinBox[ i * 5 ]   .setValue(marrCoord[i * 3 ]);
        marrDblSpinBox[ i * 5 + 1].setValue(marrCoord[i * 3 + 1]);
        marrDblSpinBox[ i * 5 + 2].setValue(marrCoord[i * 3 + 2]);
        marrDblSpinBox[ i * 5 + 3].setValue(marrV[i]);
        if ( i < mQuantActualTaskPoints)
        {
        marrDblSpinBox[ i * 5 + 4].setValue(marrCourse[i] / M_PI * 180.);
        }
        else
        {
         marrDblSpinBox[ i * 5 + 4].setValue(marrCourse[i]);
        }
    }
    for (int i=0; i< ui->tableWidget->rowCount(); i++)
    for (int j = 0; j < ui->tableWidget->columnCount(); j++)
    {
        ui->tableWidget->setCellWidget(i,j,&(marrDblSpinBox[ i * 5 + j]));
    }
    ///

    // начальная широта и долгота
    mLong0 = 30.; // долгота
    mLat0 = 30.;
    ui->doubleSpinBox->setValue(mLong0);
    ui->doubleSpinBox_2->setValue(mLat0);*/
    // полетное время
    mTFly = 1000.;
    ui->doubleSpinBox_3->setValue(mTFly);


    ///

    // создание  лопасти винта
       //  mBladeR = constBladeR; // радиус ометаемой винтом площади
        // mRadHorizHsarnir = constRadHorizHsarnir; //расстояние от центра гориз шарнира до оси вращения вала винта marrDblSpinBoxBlade[0]
        // mPofile_d0 = constPofile_d0; // высота профиля у оси вала винта marrDblSpinBoxBlade[1]
        // mPofile_d1 = constPofile_d1; // высота профиля на конце marrDblSpinBoxBlade[2]
       //  mBlade_b = constBlade_b; // хорда     marrDblSpinBoxBlade[3]
       // mBladeM = constBladeM;
      //   mBladeCX0 = constBladeCX0;


      // создание   винтов  !!!!!!!!!!!!!!!!!!!!!
    //  mQuantBlades = constQuantBlades;
      // максимально допустимый общий шаг НВ
    //  mFiMax = constFiMax;

    //  mKappaTettaMax = mFiMax ;
      // базовачя система координат привязанная к оси вращения вала  винта
      // ось Z направлена по оси вращения вала

      // координаты точки основания оси вала винта в СвСК
    //  mShaftAxeX = constShaftAxeX; // -0.37; ЭТО ДЛЯ ОТЛАДКИ!!!! // у Володко для К26 максимум 0,28 м !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //  mShaftAxeY = constShaftAxeY;
      ///
      // угол заклинения, рад., это угол между осью OZ СвСК вертолета и осью НВ
    //  mZaclinAng = constZaclinAng; //  0.05905;ЭТО ДЛЯ ОТЛАДКИ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111111111111111111

      // расстояние от основания вала винта до центра вращения верхнего винта
     // mShaftUpperL = constShaftUpperL;

      // расстояние от основания вала винта до центра вращения нижнего винта
    //  mShaftLowerL = constShaftLowerL;


      // масса ветролета
// mHelicMass = constHelicMass;



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


     ///  создание окружающей среды (атмосферы)
    // mEnvironment  = TEnvironment(mWindHor, mWindCourse, mWindVert, mTemperature0 );
     ///




}




void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\REPOSITORIES\\aircraft-model\\PROGRAMS_C++\\HELIC_OPERATION\\OutPut");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;

}




void MainWindow::on_pushButton_clicked()
{
    TEnvironment Environment (0, 0., 0., 18. );
    const bool bFullType = false,  bTypeofFullRotor = false;
    THelic Helic( bFullType, bTypeofFullRotor,constHelicMass);

    long double VAlY = 300., VAlVx = 40., VAlPsi = -0.24
            , VAlRad = 675.;


    long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};
  //  TBallanceCalc::calcBallParamsForSteadyLineMoving(Helic,Environment
     //               ,VAlY, VAlVx, VAlPsi,arrX0, arrXRez );

    int i = 0;
    long double psi =0.;
   // bool brez = TBallanceCalc::calcBallParamsForTurn(Helic,Environment
       //             ,VAlY, VAlVx, psi,VAlRad, arrX0, arrXRez );

    long double valRad = 0.;
    for ( i =0; i < 20; ++i)
    {
     //psi = -1.5 + 0.001 * static_cast<long double>(i);
    // bool brez = TBallanceCalc::calcBallParamsForTurn(Helic,Environment
        //             ,VAlY, VAlVx, psi,VAlRad, arrX0, arrXRez );

     valRad = 50. + 50. * static_cast<long double>(i);
     bool brez = TBallanceCalc::calcBallParamsForTurn(Helic,Environment
                     ,VAlY, VAlVx, psi,valRad , arrX0, arrXRez );
     if(brez)
     {
         int i0 = 0;
     }
    }


  /*  TPartHelicTraj PartHelicTraj( Helic, Environment);
    return;
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

    int iii = 0;
}

void MainWindow::on_pushButton_3_clicked()
{
    TEnvironment Environment (0, 0., 0., 18. );
    const bool bFullType = false,  bTypeofFullRotor = false;
    THelic Helic( bFullType, bTypeofFullRotor,constHelicMass);

    long double VAlY = 300.;


    long double arrX0[6] = {0.,0.,  0., 0.3, 0.,0.};
    long double arrXRez[6] = {0.};


    int i = 0;
    long double psi =0.;


    long double val_dPsi_po_dt = 0.2;
    bool brez0 = TBallanceCalc::calcBallParamsForRotation(Helic,Environment
                    ,VAlY, val_dPsi_po_dt, arrX0, arrXRez );
    /*for ( i =0; i < 20; ++i)
    {


     val_dPsi_po_dt = 2. * M_PI /10. + 2. * M_PI /10. * static_cast<long double>(i);
     bool brez = TBallanceCalc::calcBallParamsForRotation(Helic,Environment
                     ,VAlY, val_dPsi_po_dt, arrX0, arrXRez );
     if(brez)
     {
         int i0 = 0;
     }
    }*/

}

void MainWindow::on_lineEdit_cursorPositionChanged(int arg1, int arg2)
{

}
