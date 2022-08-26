#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "Receiver.h"
#include "YrWriteShapeFile.h"
#include <math.h>
#include "Comp.h"
//#include "RingedAnt.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include<QFileDialog>
#include "MatrixProccess.h"

extern bool BEZ_SHUMOV = false;;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->ui->lineEdit->setText(QString("D:\\AKIN\\ATENNA_STUDY"));
    ui->doubleSpinBox->setValue(4.8);
    ui->doubleSpinBox_2->setValue(1.4);
    ui->doubleSpinBox_4->setValue(18);
    ui->doubleSpinBox_3->setValue(7);
    ui->doubleSpinBox_5->setValue(0.);
    ui->doubleSpinBox_6->setValue(10.);
    ui->doubleSpinBox_7->setValue(11.);


}

MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow:: inputData()
{
    mQstrOutPutFold = this->ui->lineEdit->text();
    mQstrOutPutFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[mQstrOutPutFold.size()] = 0;

    mr = ui->doubleSpinBox_4->value() * 0.01/2.;
    mQuant = int(ui->doubleSpinBox_3->value() +0.1);
    ma= ui->doubleSpinBox_2->value() * 0.01;
    mlamb= ui->doubleSpinBox->value() * 0.01;
    mScan_q0= ui->doubleSpinBox_5->value() * M_PI/180.;
    mScan_e0= ui->doubleSpinBox_6->value() * M_PI/180.;

    mReceiver  = QReceiver (ma);

    mAnt = QRingedAnt (mQuant,  ma, mr) ;

    m_e_True = M_PI/2. - ui->doubleSpinBox_7->value() * M_PI/180.;
    m_e_Zv = 0.;
    m_qu_True = ui->doubleSpinBox_8->value() * M_PI/180.;
    m_qu_Zv = 0.;


}
void MainWindow::on_pushButton_clicked()
{



    wcscpy( mwchOutPutFold,  L"D:\\AKIN\\ATENNA_STUDY");

    double x = 0.;
    double step = 0.001;
    int QUantRows = 8./step;

    int QUantColsReport0 = 5;
    double *arrBuff = new double[QUantColsReport0 * QUantRows];
    memset(arrBuff, 0, sizeof(double) * QUantColsReport0 *(QUantRows ));

    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[QUantColsReport0 * iLenName];
    memset(pwcharrColNames, 0, QUantColsReport0 * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"x");
    wcscpy(&pwcharrColNames[iLenName   ],  L"Bessel");
    wcscpy(&pwcharrColNames[iLenName *2],  L"dBessel");
    wcscpy(&pwcharrColNames[iLenName *3],  L"Bx");

    wcscpy(&pwcharrColNames[iLenName *4],  L"dBx");
    for (int i =0; i < QUantRows; ++i)
    {
        double x = ((double)i) * step;
        arrBuff[i *QUantColsReport0 ] = x;;
        arrBuff[i *QUantColsReport0 +1] = fncBessel1(x);
        arrBuff[i *QUantColsReport0 +2] = fnc_dBessel1_po_dx(x);
        arrBuff[i *QUantColsReport0 +3] = fncBx(x);
        arrBuff[i *QUantColsReport0 +4] = fnc_dBx_po_dx(x);

    }




    // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,   mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.) ;
    // график
    double arrScale[100];
    for (int i = 0; i < 100;++i)
    {
     arrScale[i] = 1.;
    }




  for (int i =1; i < QUantColsReport0; i++)
  {
    TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                    ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                    ,QUantColsReport0  // - к-во переменных о корорых накоплена информация в буфере
                                    ,QUantRows //  - к-во точек
                                    ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                    ,iLenName // максимальная длина имени переменной
                                    ,0  //  номер переменной по оси X
                                    ,i  //  номер переменной по оси Y
                                    ,1. //  масштаб по оси X
                                    ,1.  // масштаб по оси Y
                                     ) ;
  }

delete []arrBuff;
delete []pwcharrColNames;

}

void MainWindow::on_pushButton_2_clicked()
{
    inputData();
    double step = 0.001;
    int QUantRows = M_PI/step;

    int QUantColsReport0 = 2;
    double *arrBuff = new double[QUantColsReport0 * QUantRows];
    memset(arrBuff, 0, sizeof(double) * QUantColsReport0 *(QUantRows ));

    double *arrBuffPolar = new double[QUantColsReport0 * QUantRows];
    memset(arrBuffPolar, 0, sizeof(double) * QUantColsReport0 *(QUantRows ));

    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[QUantColsReport0 * iLenName];
    memset(pwcharrColNames, 0, QUantColsReport0 * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"tet");
    wcscpy(&pwcharrColNames[iLenName   ],  L"Diagr");


    wchar_t *pwcharrColNames1 = new wchar_t[QUantColsReport0 * iLenName];
    memset(pwcharrColNames1, 0, QUantColsReport0 * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames1, L"tet");
    wcscpy(&pwcharrColNames1[iLenName   ],  L"PolarDiagr");


    mReceiver.collectData_forDiagrGraph(mlamb, QUantRows
                                        , step, arrBuff,arrBuffPolar);

    // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,   mwchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.) ;
    // график
    double arrScale[100];
    for (int i = 0; i < 100;++i)
    {
     arrScale[i] = 1.;
    }

  for (int i =1; i < QUantColsReport0; i++)
  {
      TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                      ,arrBuff // массив с информацией - матрица nBuffRows x nBuffCols
                                      ,QUantColsReport0  // - к-во переменных о корорых накоплена информация в буфере
                                      ,QUantRows //  - к-во точек
                                      ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                      ,iLenName // максимальная длина имени переменной
                                      ,0  //  номер переменной по оси X
                                      ,i  //  номер переменной по оси Y
                                      ,1. //  масштаб по оси X
                                      ,1.  // масштаб по оси Y
                                       ) ;

      TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                      ,arrBuffPolar // массив с информацией - матрица nBuffRows x nBuffCols
                                      ,QUantColsReport0  // - к-во переменных о корорых накоплена информация в буфере
                                      ,QUantRows //  - к-во точек
                                      ,pwcharrColNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
                                      ,iLenName // максимальная длина имени переменной
                                      ,0  //  номер переменной по оси X
                                      ,i  //  номер переменной по оси Y
                                      ,1. //  масштаб по оси X
                                      ,1.  // масштаб по оси Y
                                       ) ;
  }

delete []arrBuff;
  delete []arrBuffPolar;
delete []pwcharrColNames;
  delete []pwcharrColNames1;
}

void MainWindow::on_pushButton_3_clicked()
{
   inputData();
   /*   TComp arrQ[7];
      double q = 0., q0 = 0.;
      double e0 = 70./180. * M_PI;
      double e = e0 +0.01;
      mAnt.calcVectMeasures( mlamb,1., q, e , arrQ );
      double qv = 0., ev =0., qzv, ezv;
     mAnt.calcAmpRaznMeth( mlamb,arrQ, q0, e0, &qzv, &ezv);

     TURPolyLine pln_eTrue(  TURPointXY(e, 1000.), TURPointXY(e, 0.)) ;
     pln_eTrue.WriteSetSHPFiles(L"D:\\AKIN\\ATENNA_STUDY\\Centre\\pln_eTrue.shp",&pln_eTrue, 1);

     TURPolyLine pln_eZv(  TURPointXY(ezv, 1000.), TURPointXY(ezv, 0.)) ;
     pln_eZv.WriteSetSHPFiles(L"D:\\AKIN\\ATENNA_STUDY\\Centre\\pln_eZv.shp",&pln_eZv, 1);

return;*/

      double step = 0.001;
      int QUantRows = M_PI /step;

      int QUantCols= 3;
      double *arrBuff_q = new double[QUantCols * QUantRows];
       memset(arrBuff_q, 0, sizeof(double) * QUantCols *(QUantRows ));

       double *arrBuff_e = new double[QUantCols * QUantRows];
        memset(arrBuff_e, 0, sizeof(double) * QUantCols *(QUantRows ));
        double *arrBuffPolar_e = new double[2 * QUantRows];
         memset(arrBuffPolar_e, 0, sizeof(double) * 2 *(QUantRows ));

      // вывод графиков
      int iLenName = 30;
      wchar_t *pwcharrColNames = new wchar_t[QUantCols * iLenName];
      memset(pwcharrColNames, 0, QUantCols * iLenName*sizeof(wchar_t));

      wcscpy(pwcharrColNames, L"tet");
      wcscpy(&pwcharrColNames[iLenName   ],  L"Amp");
     wcscpy(&pwcharrColNames[iLenName  *2 ],  L"Phase");

     wchar_t *pwcharrColNames1 = new wchar_t[QUantCols * iLenName];
     memset(pwcharrColNames1, 0, QUantCols * iLenName*sizeof(wchar_t));

     wcscpy(pwcharrColNames1, L"tet");
     wcscpy(&pwcharrColNames1[iLenName   ],  L"Amp_e");
    wcscpy(&pwcharrColNames1[iLenName  *2 ],  L"Phase_e");

    wchar_t pwcharrColNames2 [100] ={0};
    wcscpy(pwcharrColNames2,  L"e");
   wcscpy(&pwcharrColNames2[iLenName   ],  L"Ant_R");






      mAnt.collectData_forDiagrGraph_e(mlamb, mScan_q0, mScan_e0 , QUantRows
                                       ,  step,arrBuff_e,arrBuffPolar_e);


     // оси  координат
      wchar_t wchAxesFileName0[300] ={0};
      wcscpy(  wchAxesFileName0,   mwchOutPutFold);
      wcscat(wchAxesFileName0, L"\\AxesArr.shp");
      TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
      ,-10000., 10000.,30.) ;
      // график
      double arrScale[100];
      for (int i = 0; i < 100;++i)
      {
       arrScale[i] = 1.;
      }


      TComp Comp = mAnt.calcSumNormalizedDiagr(mlamb,mScan_q0, M_PI/2. - mScan_e0
                                      ,mScan_q0,M_PI/2.);

    for (int i =1; i < QUantCols; i++)
    {
       /* TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                        ,arrBuff_q // массив с информацией - матрица nBuffRows x nBuffCols
                                        ,QUantCols  // - к-во переменных о корорых накоплена информация в буфере
                                        ,QUantRows //  - к-во точек
                                        ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                        ,iLenName // максимальная длина имени переменной
                                        ,0  //  номер переменной по оси X
                                        ,i  //  номер переменной по оси Y
                                        ,1. //  масштаб по оси X
                                        ,1.  // масштаб по оси Y
                                         ) ;*/
        TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                        ,arrBuff_e // массив с информацией - матрица nBuffRows x nBuffCols
                                        ,QUantCols  // - к-во переменных о корорых накоплена информация в буфере
                                        ,QUantRows //  - к-во точек
                                        ,pwcharrColNames1 //матрица с именаими переменных - матрица nBuffCols x lenName
                                        ,iLenName // максимальная длина имени переменной
                                        ,0  //  номер переменной по оси X
                                        ,i  //  номер переменной по оси Y
                                        ,1. //  масштаб по оси X
                                        ,1.  // масштаб по оси Y
                                         ) ;

    }
    for (int i =1; i < 2; i++)
    {


        TYrWriteShapeFile::WriteOneReport( mwchOutPutFold// путь к папке
                                        ,arrBuffPolar_e // массив с информацией - матрица nBuffRows x nBuffCols
                                        ,2  // - к-во переменных о корорых накоплена информация в буфере
                                        ,QUantRows //  - к-во точек
                                        ,pwcharrColNames2 //матрица с именаими переменных - матрица nBuffCols x lenName
                                        ,iLenName // максимальная длина имени переменной
                                        ,0  //  номер переменной по оси X
                                        ,i  //  номер переменной по оси Y
                                        ,1. //  масштаб по оси X
                                        ,1.  // масштаб по оси Y
                                         ) ;

    }

  delete []arrBuff_q;
    delete []arrBuff_e;
  delete []pwcharrColNames;
    delete []pwcharrColNames1;
    delete []arrBuffPolar_e;

    double valSignalAmp =1.;
    TComp *arrQ = new TComp[mQuant];
    mAnt.calcVectMeasures(mlamb,valSignalAmp, m_qu_True, m_e_True, arrQ );

    mAnt.calcAmpRaznMeth_e(mlamb,arrQ, mScan_q0, M_PI/2.- mScan_e0, &m_e_Zv);

    mAnt.calcAmpRaznMeth_q(mlamb,arrQ, mScan_q0, M_PI/2.- mScan_e0, &m_qu_Zv);

    ui->doubleSpinBox_9->setValue( (M_PI/2.- m_e_Zv)* 180./ M_PI);

    ui->doubleSpinBox_10->setValue( ( m_qu_Zv)* 180./ M_PI);
    delete []arrQ ;
}

void MainWindow::on_pushButton_4_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\AKIN\\ОТЧЕТ_ПОЗИЦИОНИРОВАНИЕ\\REZ_REPORT");
    this->ui->lineEdit->setText(strFold);

}

void MainWindow::on_pushButton_5_clicked()
{
    inputData();
    int NUM = 1000;
    float *Ztemp = new float [ NUM *7 *2] ;
    memset(Ztemp, 0,NUM *7 *2 * sizeof(float));
    if (Ztemp == NULL)
    {
    //ShowMessage (L"TYrRead::ReadFltFile\nThere are not memory for Ztemp" ) ;

    }

    FILE *fr ;
    if (	(fr=_wfopen(L"D:\\AKIN\\ОТ АБЮРАМОВИЧА 20-06-2022\\20220612_142423.proc.usbl14.x.dat",L"rb")) == NULL)
    {
    //String s_22 =  L"TYrRead::ReadFltFile\nNot possible to open file" ;
    //ShowMessage (s_22) ;

    }


    fread(Ztemp,sizeof(float),NUM *7 *2,fr) ;
    fclose(fr);
    int num = 0;
    bool bStop = false;
    for (int i =0; i < NUM; ++i)
    {
        float *p = &Ztemp[i * 14];
        for (int j =0; j < 14; ++j)
        {
           if( NormVect(p, 14) < 0.001)
               bStop= true;
               break;;
        }
        if (bStop)
        {
            break;
        }
        num++;
    }

    TComp *arrcmp = new TComp[num *7];

    for (int i =0; i < num *7; ++i)
    {
      float *p =  &Ztemp[2 * i];
      arrcmp[i] = TComp(((double)Ztemp[2 * i]),((double)Ztemp[2 * i +1]) ) ;

    }
    delete []Ztemp;

   //--------------------------------------------------



    int QUantCols= 15;
    int QuantRows =num;
    double *arrBuff = new double[QUantCols * QuantRows];
     memset(arrBuff, 0, sizeof(double) * QUantCols * num);



    // вывод графиков
    int iLenName = 30;
    wchar_t *pwcharrColNames = new wchar_t[QUantCols * iLenName];
    memset(pwcharrColNames, 0, QUantCols * iLenName*sizeof(wchar_t));

    wcscpy(pwcharrColNames, L"n");
    wcscpy(&pwcharrColNames[iLenName    ],   L"Amp0");
    wcscpy(&pwcharrColNames[iLenName  *2 ],  L"Amp1");
    wcscpy(&pwcharrColNames[iLenName  *3 ],  L"Amp2");
    wcscpy(&pwcharrColNames[iLenName  *4 ],  L"Amp3");
    wcscpy(&pwcharrColNames[iLenName  *5 ],  L"Amp4");
    wcscpy(&pwcharrColNames[iLenName  *6 ],  L"Amp5");
    wcscpy(&pwcharrColNames[iLenName  *7 ],  L"Amp6");
    wcscpy(&pwcharrColNames[iLenName  *8 ],  L"Arg0");
    wcscpy(&pwcharrColNames[iLenName  *9 ],  L"Arg1");
    wcscpy(&pwcharrColNames[iLenName  *10],  L"Arg2");
    wcscpy(&pwcharrColNames[iLenName  *11],  L"Arg3");
    wcscpy(&pwcharrColNames[iLenName  *12],  L"Arg4");
    wcscpy(&pwcharrColNames[iLenName  *13],  L"Arg5");
    wcscpy(&pwcharrColNames[iLenName  *14],  L"Arg6");

    for (int i =0; i < QuantRows; ++i)
    {
        double *p = &arrBuff[i *QUantCols ];
        TComp *pcmp = &arrcmp[7 * i];
        p[0] = ((double)i);

        for (int j = 0; j < 7; ++j)
        {
          p[1 + j] =  pcmp[j].modul();
          p[8 + j] =  pcmp[j].phase();
          if (p[8 + j] < 0.)
          {
             p[8 + j] += 2. * M_PI;
          }
        }
        double *pp =&p[8];
        for (int j = 1; j < 7;++j)
        {
           pp[j] -= pp[0];
           if (pp[j] < 0.)
           {
               pp[j] += M_PI * 2.;
           }
        }
        pp[0] = 0.;
    }





    wchar_t wchOutPutFold[400] = {0};
    wcscpy(  wchOutPutFold,   mwchOutPutFold);
    wcscat(wchOutPutFold, L"\\DIAGR");
    _wmkdir(wchOutPutFold);

   // оси  координат
    wchar_t wchAxesFileName0[300] ={0};
    wcscpy(  wchAxesFileName0,   wchOutPutFold);
    wcscat(wchAxesFileName0, L"\\AxesArr.shp");
    TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
    ,-10000., 10000.,30.) ;
    // график
    double arrScale[100];
    for (int i = 0; i < 100;++i)
    {
     arrScale[i] = 1.;
    }
    for (int i = 1; i < 8;++i)
    {
     arrScale[i] = 10.;
    }




  for (int i =1; i < QUantCols; i++)
  {
      TYrWriteShapeFile::WriteOneReport( wchOutPutFold// путь к папке
                                      ,arrBuff// массив с информацией - матрица nBuffRows x nBuffCols
                                      ,QUantCols  // - к-во переменных о корорых накоплена информация в буфере
                                      ,QuantRows //  - к-во точек
                                      ,pwcharrColNames //матрица с именаими переменных - матрица nBuffCols x lenName
                                      ,iLenName // максимальная длина имени переменной
                                      ,0  //  номер переменной по оси X
                                      ,i  //  номер переменной по оси Y
                                      ,1. //  масштаб по оси X
                                      ,arrScale[i] // масштаб по оси Y
                                       ) ;

  }


delete []arrBuff;

delete []pwcharrColNames;

  delete []arrcmp;
  double e = asin(2./2.25);
  TComp arrQ[7];
  mAnt.calcVectMeasures(mlamb,1., 0., e, arrQ );
  double arrArg[7] ={0.};
  for (int i = 0; i< 7; ++i)
  {
     arrArg[i] = arrQ[i].phase();
     if (arrArg[i] < 0.)
     {
         arrArg[i] += M_PI * 2;
     }
  }
  int i0 = 0;
  for (int i = 1; i< 7; ++i)
  {
      int icur = (i0 + i)%7;
     arrArg[icur] -= arrArg[i0];
     if (arrArg[icur]< 0.)
     {
         arrArg[icur] += M_PI * 2;
     }
  }
  arrArg[i0] = 0.;

  TURPolyLine parrPln[7];
  for (int i =0; i< 7; ++i)
  {
    parrPln[i] = TURPolyLine(TURPointXY(0.,arrArg[i]), TURPointXY(1000.,arrArg[i]));
  }

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args0.shp");
  parrPln[0].WriteSetSHPFiles(wchAxesFileName0, &parrPln[0], 1);



  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args1.shp");
  parrPln[1].WriteSetSHPFiles(wchAxesFileName0, &parrPln[1], 1);

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args2.shp");
  parrPln[2].WriteSetSHPFiles(wchAxesFileName0, &parrPln[2], 1);

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args3.shp");
  parrPln[3].WriteSetSHPFiles(wchAxesFileName0, &parrPln[3], 1);

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args4.shp");
  parrPln[4].WriteSetSHPFiles(wchAxesFileName0, &parrPln[4], 1);

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args5.shp");
  parrPln[5].WriteSetSHPFiles(wchAxesFileName0, &parrPln[5], 1);

  wcscpy(  wchAxesFileName0,   wchOutPutFold);
  wcscat(wchAxesFileName0, L"\\ideal_Args6.shp");
  parrPln[6].WriteSetSHPFiles(wchAxesFileName0, &parrPln[6], 1);

}
