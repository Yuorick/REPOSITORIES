#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

#include "MatrixProccess.h"
#include "Gauss.h"



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{    
    ui->setupUi(this);

    setupFirst();
    fillWindow();
    fillCalculatorArraysPrevious();
    fillCalculatorFrame();

}
void MainWindow:: inputData()
{

    // 1
    for (int i = 0; i < NUM_GADG_PARAMS; ++i)
    {
     marrXTrue_RLK_PSK[i] =  ui->tableWidget->item(0,i)->text().toDouble();
     marrXZv_RLK_PSK[i] =  ui->tableWidget->item(1,i)->text().toDouble();
    }

    // 2
    for (int i = 0; i < NUM_GADG_PARAMS; ++i)
    {
     marrXTrue_IU_PSK[i] =  ui->tableWidget_2->item(0,i)->text().toDouble();
     marrXZv_IU_PSK[i] =  ui->tableWidget_2->item(1,i)->text().toDouble();
    }





}



MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::setupFirst()
{

    // вектор истинных параметров позиционирования РЛК в ПСК
    marrXTrue_RLK_PSK[0] = -5.5; // X
    marrXTrue_RLK_PSK[1] = 2.5;  // Y
    marrXTrue_RLK_PSK[2] = 20.5;  // Z
    marrXTrue_RLK_PSK[3] = M_PI/4.+0.003;  // Bet
    marrXTrue_RLK_PSK[4] =10./ 180. * M_PI + 0.003;  // Eps
    marrXTrue_RLK_PSK[5] = 0.002;  // Alf


   // вектор первичных оценок параметров позиционирования РЛК в ПСК
    //marrXZv_RLK_PSK[0] = 1.;
    //marrXZv_RLK_PSK[1] = 10.;
    //marrXZv_RLK_PSK[2] = 25.;
    marrXZv_RLK_PSK[0] = -5.;
    marrXZv_RLK_PSK[1] = 2.;
    marrXZv_RLK_PSK[2] = 20.;
    marrXZv_RLK_PSK[3] = M_PI/4.;
    marrXZv_RLK_PSK[4] = 10./ 180. * M_PI;
    marrXZv_RLK_PSK[5] = 0.;

    // вектор истинных параметров позиционирования ИУ в ПСК
    marrXTrue_IU_PSK[0] = -0.5;
    marrXTrue_IU_PSK[1] = 14.5;
    marrXTrue_IU_PSK[2] = 19.5;
    marrXTrue_IU_PSK[3] = - M_PI/4.-0.003;
    marrXTrue_IU_PSK[4] = 10./ 180. * M_PI - 0.003;
    marrXTrue_IU_PSK[5] = -0.005;


   // вектор первичных оценок параметров позиционирования ИУ в ПСК
    //marrXZv_IU_PSK[0] = 5.;
    //marrXZv_IU_PSK[1] = 16.;
    //marrXZv_IU_PSK[2] = 4.;
    marrXZv_IU_PSK[0] = 0.;
    marrXZv_IU_PSK[1] = 15.;
    marrXZv_IU_PSK[2] = 20.;
    marrXZv_IU_PSK[3] = - M_PI/4.;
    marrXZv_IU_PSK[4] = 10./ 180. * M_PI;
    marrXZv_IU_PSK[5] = 0.;

  /*  // вектор истинных параметров позиционирования РЛК в ПСК
    marrXTrue_RLK_PSK[0] = 2.;
    marrXTrue_RLK_PSK[1] = 11;
    marrXTrue_RLK_PSK[2] = 26.;
    marrXTrue_RLK_PSK[3] = 0.001;
    marrXTrue_RLK_PSK[4] = 0.003;


   // вектор первичных оценок параметров позиционирования РЛК в ПСК
    marrXZv_RLK_PSK[0] = 2.;
    marrXZv_RLK_PSK[1] = 11.;
    marrXZv_RLK_PSK[2] = 26.;
    marrXZv_RLK_PSK[3] = 0.;
    marrXZv_RLK_PSK[4] = 0.;

    // вектор истинных параметров позиционирования ИУ в ПСК
    marrXTrue_IU_PSK[0] = 6.;
    marrXTrue_IU_PSK[1] = 15.;
    marrXTrue_IU_PSK[2] = 3.;
    marrXTrue_IU_PSK[3] = 0.;
    marrXTrue_IU_PSK[4] = 0.;


   // вектор первичных оценок параметров позиционирования ИУ в ПСК
    marrXZv_IU_PSK[0] = 6.;
    marrXZv_IU_PSK[1] = 15.;
    marrXZv_IU_PSK[2] = 3.;
    marrXZv_IU_PSK[3] = -0.003;
    marrXZv_IU_PSK[4] = -0.002;
*/





}
//-----------------------------------------------------

void MainWindow::fillWindow()
{

    //1
   int ia =0;
   double temp = 0.;

   for (int i=0; i<  ui->tableWidget->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                 {
                     QTableWidgetItem* item = new QTableWidgetItem;
                     item->setText(QString::number(0.));

                     ui->tableWidget->setItem(i,j,item);
                 }

    for (int j =0; j < ui->tableWidget->columnCount(); ++j)
    {
         ia = ((int)(marrXTrue_RLK_PSK[j] * 1000.));
         temp = ((double)ia)/1000.;
     ui->tableWidget->item(0,j)->setText(QString::number(temp));

     ia = ((int)(marrXZv_RLK_PSK[j] * 1000.));
     temp = ((double)ia)/1000.;
     ui->tableWidget->item(1,j)->setText(QString::number(temp));
    }

    ///

    // 3
    for (int i=0; i<  ui->tableWidget_2->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = nullptr;
                     ptwi0 = new QTableWidgetItem(QString::number(0.));
                     ui->tableWidget_2->setItem(i,j,ptwi0);
                 }

    for (int j =0; j < ui->tableWidget_2->columnCount(); ++j)
    {
        ia = ((int)(marrXTrue_IU_PSK[j] * 1000.));
        temp = ((double)ia)/1000.;
     ui->tableWidget_2->item(0,j)->setText(QString::number(temp));

     ia = ((int)(marrXZv_IU_PSK[j] * 1000.));
     temp = ((double)ia)/1000.;
     ui->tableWidget_2->item(1,j)->setText(QString::number(temp));
    }
    ///

}



// замполнение массивов для фрема калькулятора превоначальное
void MainWindow::fillCalculatorArraysPrevious()
{
    // палубные углы - угол курса, килевой качки, бортовой качки
    memset(marrDeckAngles, 0,3 * sizeof(double));

    // массив данных пеленга цели - угол пеленга, дальность, высота
    marrBaring[0] = marrDeckAngles[0] + marrXTrue_RLK_PSK[3];
    marrBaring[1] = 5000.;
    marrBaring[2] = 0.;

    // массив углов наведения истинный - Betta, Eps
    memset(marrCntrlAngTrue_V1, 0,2 * sizeof(double));

    // массив углов наведения оценка первоначальная - Betta, Eps
    memset(marrCntrlAngZv_V1, 0,2 * sizeof(double));

    // матрица частных производных углов наведения по параметрам позиционирования
    memset(marr_dVIU_po_dX_V1, 0,2 * LENX_V1 * sizeof(double));
    ///
    // массив углов наведения истинный - Betta, Eps
    memset(marrCntrlAngTrue_V2, 0,2 * sizeof(double));

    // массив углов наведения оценка первоначальная - Betta, Eps
    memset(marrCntrlAngZv_V2, 0,2 * sizeof(double));

    // матрица частных производных углов наведения по параметрам позиционирования
    memset(marr_dVIU_po_dX_V2, 0,2 * NUM_GADG_PARAMS * sizeof(double));
}
//---------------------------------------------
//заполнение фрейма калькулятора
void MainWindow::fillCalculatorFrame()
{
    int ia = 0;
    double temp = 0.;
    //1 палубные углы

      for (int j = 0; j < ui->tableWidget_5->columnCount(); j++)
        {
        QTableWidgetItem* ptwi0 = new QTableWidgetItem;
        ptwi0 ->setText(QString::number(0.));
        ui->tableWidget_5->setItem(0,j,ptwi0);
        }

        for (int j =0; j < ui->tableWidget_5->columnCount(); ++j)
        {
           double coeff  = 180. / M_PI;
           ui->tableWidget_5->item(0,j)->setText(QString::number(marrDeckAngles[j] * coeff));
        }

        ///

        //2 таблица параметров пеленга цели


        for (int j = 0; j < 3; j++)
          {
          QTableWidgetItem* ptwi0 = new QTableWidgetItem;
          ptwi0 ->setText(QString::number(0.));
          ui->tableWidget_9->setItem(0,j,ptwi0);
          }

          for (int j =0; j < 3; ++j)
          {
             double coeff = 1.;//;(j != 0)?1:1./ M_PI * 180.;
             ui->tableWidget_9->item(0,j)->setText(QString::number(marrBaring[j] *coeff));
          }

          // ПО ВАРИАНТУ 1
            //3  таблица углов наведения вычисленная
          for (int i =0; i < 3; ++i)
            for (int j = 0; j < ui->tableWidget_7->columnCount(); j++)
            {
            QTableWidgetItem* ptwi0 = new QTableWidgetItem;
            ptwi0 ->setText(QString::number(0.));
            ui->tableWidget_7->setItem(i,j,ptwi0);
            }

            for (int j =0; j < ui->tableWidget_7->columnCount(); ++j)
            {
                ia = (10000.*marrCntrlAngTrue_V1[j]);
                temp = (double)(ia/ 10000.);
            ui->tableWidget_7->item(0,j)->setText(QString::number(temp));


            ia = (10000.*marrCntrlAngZv_V1[j]);
            temp = (double)(ia/ 10000.);
            ui->tableWidget_7->item(1,j)->setText(QString::number(temp));

            ia = (10000.*(marrCntrlAngZv_V1[j] -marrCntrlAngTrue_V1[j]));
            temp = (double)(ia/ 10000.);
            ui->tableWidget_7->item(2,j)->setText(QString::number(temp ));
            }


            ///

            //3  таблица дифференциалов вычисленная
            for (int i =0; i < 2; ++i)
            for (int j = 0; j < ui->tableWidget_8->columnCount(); j++)
            {
            QTableWidgetItem* ptwi0 = new QTableWidgetItem;
            ptwi0 ->setText(QString::number(0.));
            ui->tableWidget_8->setItem(i,j,ptwi0);
            }

            for (int i =0; i < 2; ++i)
            for (int j = 0; j < ui->tableWidget_8->columnCount(); j++)
            {
                int ia = ((int)(marr_dVIU_po_dX_V1[i * LENX_V1 + j] * 1000.));
                double temp = ((double)ia)/1000.;
                //ui->tableWidget_8->item(i,j)->setText(QString::number(marr_dVIU_po_dX[i * 5 + j] ));
                ui->tableWidget_8->item(i,j)->setText(QString::number(temp ));
            }
            ///


            // ПО ВАРИАНТУ 2
              //3  таблица углов наведения вычисленная
            for (int i =0; i < 3; ++i)
              for (int j = 0; j < ui->tableWidget_10->columnCount(); j++)
              {
              QTableWidgetItem* ptwi0 = new QTableWidgetItem;
              ptwi0 ->setText(QString::number(0.));
              ui->tableWidget_10->setItem(i,j,ptwi0);
              }

              for (int j =0; j < ui->tableWidget_10->columnCount(); ++j)
              {
              ui->tableWidget_10->item(0,j)->setText(QString::number(marrCntrlAngTrue_V2[j] ));
              ui->tableWidget_10->item(1,j)->setText(QString::number(marrCntrlAngZv_V2[j] ));
              ui->tableWidget_10->item(2,j)->setText(QString::number(marrCntrlAngZv_V2[j] - marrCntrlAngTrue_V2[j]));
              }




              ///

              //3  таблица дифференциалов вычисленная
              for (int i =0; i < 2; ++i)
              for (int j = 0; j < ui->tableWidget_11->columnCount(); j++)
              {
              QTableWidgetItem* ptwi0 = new QTableWidgetItem;
              ptwi0 ->setText(QString::number(0.));
              ui->tableWidget_11->setItem(i,j,ptwi0);
              }

              for (int i =0; i < 2; ++i)
              for (int j = 0; j < ui->tableWidget_11->columnCount(); j++)
              {

                  int ia = ((int)(marr_dVIU_po_dX_V2[i * NUM_GADG_PARAMS + j] * 1000.));
                  double temp = ((double)ia)/1000.;
                  //ui->tableWidget_8->item(i,j)->setText(QString::number(marr_dVIU_po_dX[i * 5 + j] ));
                  ui->tableWidget_11->item(i,j)->setText(QString::number(temp ));
              }


}



void MainWindow::on_pushButton_2_clicked()
{
    inputData();
    calculatorInput();
    // 1.
    double arrTargS_KGSK[3] = {0.};
    arrTargS_KGSK[0] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * sin(marrBaring[0]);
    arrTargS_KGSK[1] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * cos(marrBaring[0]);
    arrTargS_KGSK[2] = marrBaring[2];
    ///

    // ОБСЧЕТ 2-ГО ВАРИАНТА
    double valR = 0., valr = 0.;

    QAdjustment::calc_VIU_and_dVIU_po_dX(arrTargS_KGSK ,marrDeckAngles
                ,marrXTrue_IU_PSK,valR,marrCntrlAngTrue_V2[0],marrCntrlAngTrue_V2[1],marr_dVIU_po_dX_V2);

    QAdjustment::calc_VIU(arrTargS_KGSK ,marrDeckAngles
                ,marrXZv_IU_PSK,valr,marrCntrlAngZv_V2[0],marrCntrlAngZv_V2[1]);

/////////////////////////////////////////////////////////////////////////////////////
   double  arrVZv_RLK[3] = {0.};
   QAdjustment::calc_VIU(arrTargS_KGSK ,marrDeckAngles
               ,marrXTrue_RLK_PSK,arrVZv_RLK[1],arrVZv_RLK[0],arrVZv_RLK[2]);

    // ОБСЧЕТ 1-ГО ВАРИАНТА
   QAdjustment::recalc_VIU_from_RLK_to_IU_and_dVIU_po_dX(arrVZv_RLK
        ,marrXTrue_RLK_PSK,marrXTrue_IU_PSK, marrCntrlAngTrue_V1[0], marrCntrlAngTrue_V1[1], marr_dVIU_po_dX_V1);



    QAdjustment::recalc_VIU_from_RLK_to_IU(arrVZv_RLK
      ,marrXZv_RLK_PSK,marrXZv_IU_PSK, marrCntrlAngZv_V1[0], marrCntrlAngZv_V1[1]);

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            marr_dVIU_po_dX_V1[ i * LENX_V1 + j] *= 1000.;
            marr_dVIU_po_dX_V2[ i * NUM_GADG_PARAMS + j] *= 1000.;
        }
    }

  fillCalculatorFrame();

//--------------------
  double valQc = 0., valEc = 0.;
  QAdjustment::calc_VIU(arrTargS_KGSK ,marrDeckAngles
              ,marrXTrue_IU_PSK,valr,valQc,valEc);
  double temmp0 = -cos(marrXTrue_IU_PSK[5] * cos(marrXTrue_IU_PSK[4])) +tan(valEc)
       *(cos(valQc)*sin(marrXTrue_IU_PSK[4]) - sin(valQc)* cos(marrXTrue_IU_PSK[4])* sin(marrXTrue_IU_PSK[5]));

  double temp1 = -cos(valQc)* cos(marrXTrue_IU_PSK[4])*sin(marrXTrue_IU_PSK[5])-
          sin(valQc)*sin(marrXTrue_IU_PSK[4]);
  //
  double temmp2 = -tan(valEc)*cos(valQc);
  double temp3 = sin(valQc);

  //
  double temmp4 =sin(marrXTrue_IU_PSK[5]) -sin(valQc)*tan(valEc)*cos(marrXTrue_IU_PSK[5]);
  double temp5 = -cos(valQc)*cos(marrXTrue_IU_PSK[5]);


  int yy=0;
}
//---------------------------------------------
void MainWindow::calculatorInput()
{
    for (int i = 0; i < 3; ++i)
    {
     double coeff  = 180. / M_PI;
     marrDeckAngles[i] =  (ui->tableWidget_5->item(0,i)->text().toDouble())/coeff;

    }
    for (int i = 0; i < 3; ++i)
    {
     double coeff = 1.;//(i != 0)?1:1./ M_PI * 180.;
     marrBaring[i] = ( ui->tableWidget_9->item(0,i)->text().toDouble())/coeff;

    }

}
