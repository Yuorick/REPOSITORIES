#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>

#include "MatrixProccess.h"
#include "Gauss.h"
#include <QErrorMessage>

extern  bool BEZ_SHUMOV = true;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    mbFirstTime = true;
    ui->setupUi(this);
    setupFirst();
    fillWindow();
    mbFirstTime = false;

  //  mbCalculationIsDone = false;



    mbStar = false;



}
void MainWindow:: inputData()
{

    //0 калькулятор установка
    //mbCalculationIsDone = false;




    // 1
    for (int i = 0; i < QUANT_GDG_PRMS; ++i)
    {
     marrXTrue_IU_PSK[i] =  ui->tableWidget->item(0,i)->text().toDouble();
     marrXZv_IU_PSK[i] =  ui->tableWidget->item(1,i)->text().toDouble();
    }

    // 2


      // Кол-во контрольных объектов (КО)
    mNumCtrlObj = (int)(ui->doubleSpinBox_2->value() + 0.1);

    // высота КО
    mCtrlObjHeight =ui->doubleSpinBox->value();

    // массив исходных данных КО из входной таблицы
    for (int i = 0; i <mNumCtrlObj; ++i)
        for (int j =0; j < 5; ++j)
        {
          double coeff = (j == 1)?1:1./ M_PI * 180.;
          mparrDataCtrlObj[i * 5 + j]= (ui->tableWidget_3->item(i,j)->text().toDouble())/  coeff;
        }

    // СКЗ ошибки по углу ОЭК
    mAngleSigIU = ui->doubleSpinBox_14->value() *0.001;

    // СКЗ ошибки дальномера
    mDistSigIU = 0.1;

    //признак случайных шумов измерений
    mbNoise= ui->checkBox_2->isChecked();

    // Коррел матрица ошибок первичных измерений РЛК
    memset(marrCorMtrxRLK, 0, 9 * sizeof(double));
    marrCorMtrxRLK[0] = marrCorMtrxRLK[8] = (0.00025)*(0.00025);
    marrCorMtrxRLK[4] = 100.;

    double valTSinsDelay = ui->doubleSpinBox_12->value() * 0.001;
    double valSigDelay = valTSinsDelay / sqrt(12.);

    double valSigQ = ui->doubleSpinBox_7->value()* 0.001;
    double valVSigQ = ui->doubleSpinBox_8->value()* 0.001;

    double valSigPsi = ui->doubleSpinBox_9->value()* 0.001;
    double valVSigPsi = ui->doubleSpinBox_10->value()* 0.001;

    double valSigTet = ui->doubleSpinBox_11->value()* 0.001;
    double valVSigTet = ui->doubleSpinBox_13->value()* 0.001;


    marrDispEilerCntrKP[0] = valSigQ * valSigQ + (valVSigQ * valSigDelay) *(valVSigQ * valSigDelay);
    marrDispEilerCntrKP[1] = valSigPsi * valSigPsi + (valVSigPsi * valSigDelay) *(valVSigPsi * valSigDelay);
    marrDispEilerCntrKP[2] = valSigTet * valSigTet + (valVSigTet * valSigDelay) *(valVSigTet * valSigDelay);

    mVessSKZ = ui->doubleSpinBox_5->value();

    mCntrlObjSKZ = ui->doubleSpinBox_6->value();



     marrOrtStar[0] = 0.;
     marrOrtStar[1] = cos(66./ 180. *M_PI );
     marrOrtStar[2] = sin(66./ 180. *M_PI );

     if (ui->checkBox_5->checkState() == 0)
     {
       mbStar = false;
     }
     else
     {
         mbStar = true;
     }


}

void MainWindow::on_pushButton_clicked()
{

    if(ui->checkBox_4->isChecked())
    {
        on_checkBox_4_stateChanged(2);
    }
    inputData();
    if (!checkObj())
    {
     (new QErrorMessage(this))->showMessage("ЗАПОЛНИТЕ ТАБЛИЦУ КОНТРОЛЬНЫХ ОБЪЕКТОВ");
        return;
    }
    double arrXTrue_RLK_PSK[QUANT_GDG_PRMS] = {0.},arrXZv_RLK_PSK[QUANT_GDG_PRMS] = {0.};
    // подготовка массива КО
    QContObj6 *pqtrlObjectArr = new QContObj6[mNumCtrlObj];

    if (mbStar)
    {
        mVessSKZ = 0.;
        mCntrlObjSKZ =0.;
        memset(marrXTrue_IU_PSK, 0, 3 * sizeof(double));
        memset(marrXZv_IU_PSK, 0, 3 * sizeof(double));
        for (int i = 0; i < mNumCtrlObj; ++i)
         {
             double arrS_KGSK[3] = {0.};
             double *p = &(mparrDataCtrlObj[5 * i]);
             arrS_KGSK[0] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * sin(p[0]);
             arrS_KGSK[1] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * cos(p[0]);
             arrS_KGSK[2] = mCtrlObjHeight;


             pqtrlObjectArr[i] =  QContObj6(&(p[2]), marrDispEilerCntrKP
                ,arrXTrue_RLK_PSK,marrXTrue_IU_PSK,arrXZv_RLK_PSK,marrXZv_IU_PSK
                ,marrOrtStar, marrCorMtrxRLK, mAngleSigIU, mDistSigIU
                     , mVessSKZ,mCntrlObjSKZ);
         }
    }
    else
    {
       for (int i = 0; i < mNumCtrlObj; ++i)
        {
            double arrS_KGSK[3] = {0.};
            double *p = &(mparrDataCtrlObj[5 * i]);
            arrS_KGSK[0] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * sin(p[0]);
            arrS_KGSK[1] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * cos(p[0]);
            arrS_KGSK[2] = mCtrlObjHeight;


            pqtrlObjectArr[i] =  QContObj6(&(p[2]), marrDispEilerCntrKP
               ,arrXTrue_RLK_PSK,marrXTrue_IU_PSK,arrXZv_RLK_PSK,marrXZv_IU_PSK
               ,arrS_KGSK, marrCorMtrxRLK, mAngleSigIU, mDistSigIU
                    , mVessSKZ,mCntrlObjSKZ);
        }
    }


    double rez =-1.;
    if (!mbStar)
    {

        if (ui->checkBox_3->isChecked())
        {
        rez = QContObj6::imitate_and_estimateAngs_Var2(pqtrlObjectArr,mNumCtrlObj
               , mbNoise, marrXEst_IU_PSK, marrMtrxK);
        }
        else
        {
        rez = QContObj6::imitate_and_estimateParams_Var2(pqtrlObjectArr,mNumCtrlObj
              , mbNoise, marrXEst_IU_PSK, marrMtrxK);
        }
    }
    else
    {
        rez = QContObj6::imitate_and_estimateStarProcessing(pqtrlObjectArr,mNumCtrlObj
               , mbNoise, marrXEst_IU_PSK, marrMtrxK);
    }


    int ia = 0;
    double temp = 0.;

    for (int j = 0; j < QUANT_GDG_PRMS; ++j)
    {
       QTableWidgetItem* item0 = new QTableWidgetItem;
        item0 ->setText(QString::number(marrXEst_IU_PSK[j]));
        ui->tableWidget->setItem(2,j,item0);

        QTableWidgetItem* item2 = new QTableWidgetItem;
        double coeff = 1.;
        if (j >2)
        {
            coeff = 1000.;
        }
         ia = int(coeff *(marrXEst_IU_PSK[j] -marrXTrue_IU_PSK[j]) * 10000.);
         temp = ia/ 10000.;
        item2 ->setText(QString::number(temp));
        ui->tableWidget->setItem(3,j,item2);



        QTableWidgetItem* item3 = new QTableWidgetItem;
        if (j < 3)
        {
            ia = int(sqrt(marrMtrxK[j * QUANT_GDG_PRMS+ j])*100.);
            temp = ia/ 100.;
        item3 ->setText(QString::number(temp));
        }
        else
        {
            ia = int(sqrt(marrMtrxK[j * QUANT_GDG_PRMS+ j]) * 1000. * 100.);
            temp = ia/ 100.;
         item3 ->setText(QString::number(temp));
        }
        ui->tableWidget->setItem(4,j,item3);


    }

    delete []pqtrlObjectArr;

  //  mbCalculationIsDone = true;



}



MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::setupFirst()
{

    // вектор истинных параметров позиционирования ИУ в ПСК
    marrXTrue_IU_PSK[0] = 2.;
    marrXTrue_IU_PSK[1] = 11;
    marrXTrue_IU_PSK[2] = 26.;
    marrXTrue_IU_PSK[3] = -0.5*M_PI/ 180;
    marrXTrue_IU_PSK[4] = 0.5*M_PI/ 180;
    marrXTrue_IU_PSK[5] = 0.005;


   // вектор первичных оценок параметров позиционирования ИУ в ПСК
    //marrXZv_IU_PSK[0] = 1.;
    //marrXZv_IU_PSK[1] = 10.;
    //marrXZv_IU_PSK[2] = 25.;
    marrXZv_IU_PSK[0] = 2.;
    marrXZv_IU_PSK[1] = 11.;
    marrXZv_IU_PSK[2] = 26.;
    marrXZv_IU_PSK[3] = 0.;
    marrXZv_IU_PSK[4] = 0.;
    marrXZv_IU_PSK[5] = 0.;




    // корреляционная матрица ошибок оценивания
    memset(marrMtrxK, 0, QUANT_GDG_PRMS*QUANT_GDG_PRMS * sizeof(double));


    // Кол-во контрольных объектов (КО)
    mNumCtrlObj = 4;



    //признак случайных шумов измерений
    mbNoise = false;

    // массив исходных данных КО из входной таблицы
    memset(mparrDataCtrlObj, 0, MAX_QUANT_CNTRL_OBJ *5 * sizeof(double));
    mparrDataCtrlObj[0] = 0.;
    mparrDataCtrlObj[1] = 6000.;
    mparrDataCtrlObj[3] = -10./ 180. * M_PI;
    mparrDataCtrlObj[5] = M_PI/3;
    mparrDataCtrlObj[6] = 6000.;
    mparrDataCtrlObj[11] = 8000.;
    mparrDataCtrlObj[15] = M_PI/2.;
    mparrDataCtrlObj[16] = 10000.;

    // высота КО
    mCtrlObjHeight = 50.;

    mVessSKZ = 3.;


    mCntrlObjSKZ = 1.;



}
//-----------------------------------------------------

void MainWindow::fillWindow()
{

    //1
    ui->tableWidget_3->setRowCount(MAX_QUANT_CNTRL_OBJ );
    for (int i=0; i<  ui->tableWidget_3->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget_3->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = new QTableWidgetItem;
                     ptwi0 ->setText(QString::number(0.));
                     ui->tableWidget_3->setItem(i,j,ptwi0);
                 }
    for (int i = 0; i < mNumCtrlObj; ++i)
        for (int j =0; j < ui->tableWidget_3->columnCount(); ++j)
        {
         double coeff = (j == 1)?1:1./ M_PI * 180.;
         ui->tableWidget_3->item(i,j)->setText(QString::number(mparrDataCtrlObj[i * 5 +j] * coeff));
        }
    //2
   // QTableWidget table(2,5) ;
   // (ui->tableWidget) = &table ;
   for (int i=0; i<  ui->tableWidget->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                 {
                     QTableWidgetItem* item = new QTableWidgetItem;
                     item->setText(QString::number(0.));

                     ui->tableWidget->setItem(i,j,item);
                 }

    for (int j =0; j < ui->tableWidget->columnCount(); ++j)
    {
     ui->tableWidget->item(0,j)->setText(QString::number(marrXTrue_IU_PSK[j]));
     ui->tableWidget->item(1,j)->setText(QString::number(marrXZv_IU_PSK[j]));
    }

    ///

    // 3


        //3
        ui->doubleSpinBox_2->setValue(mNumCtrlObj);
        ///

        //4
        ui->doubleSpinBox->setValue(mCtrlObjHeight);

        ///
        //5

        //6
        ui->checkBox_2->setChecked(false);


        ui->doubleSpinBox_5->setValue(mVessSKZ);

        ui->doubleSpinBox_6->setValue(mCntrlObjSKZ);




        //ui->comboBox->setCurrentIndex(mTypeOfTask);

}



void MainWindow::on_doubleSpinBox_2_valueChanged(const QString &arg1)
{

}



void MainWindow::on_tableWidget_2_cellChanged(int row, int column)
{
   if(!mbFirstTime)
   {
  inputData();

   }
}

void MainWindow::on_tableWidget_cellChanged(int row, int column)
{
    if(!mbFirstTime)
    {
   inputData();

    }
}

void MainWindow::on_checkBox_4_stateChanged(int arg1)
{
    if(arg1 == 2)
    {
      memset(mparrDataCtrlObj, 0, sizeof(double) *MAX_QUANT_CNTRL_OBJ * 5);
      for(int i =0 ; i < mNumCtrlObj; ++i )
      {

        double *p = &(mparrDataCtrlObj [ i * 5]);
        // курс корабля
        p[2] =  (getRand01( )- 0.5)*2. *(10./ 180. * M_PI); ;
        // пеленг цели
        p[0] = p[2] + (getRand01( )- 0.5)*2. * 0.9 * M_PI /2.;
        // дальность цели
        double valdmin = ui->doubleSpinBox_3->value();
        double valdmax = ui->doubleSpinBox_4->value();
        double valdMiddle = (valdmin + valdmax)/2.;
        p[1] = valdMiddle + (getRand01( )- 0.5)*2. * (valdmax - valdMiddle);

        // угол кормовой качки
        p[3] =  (getRand01( )- 0.5)*2. *(10./ 180. * M_PI); ;
        // угол бортовой  качки
        p[4] =  (getRand01( )- 0.5)*2. *(18./ 180. * M_PI);

      }

    }
    for (int i=0; i<  ui->tableWidget_3->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget_3->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = new QTableWidgetItem;
                     ptwi0 ->setText(QString::number(0.));
                     ui->tableWidget_3->setItem(i,j,ptwi0);
                 }
    for (int i = 0; i < ui->tableWidget_3->rowCount(); ++i)
        for (int j =0; j < ui->tableWidget_3->columnCount(); ++j)
        {
         double coeff = (j == 1)?1:1./ M_PI * 180.;
         ui->tableWidget_3->item(i,j)->setText(QString::number(mparrDataCtrlObj[i * 5 +j] * coeff));
        }
}

void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
 if (arg1 > MAX_QUANT_CNTRL_OBJ)
 {
   ui->doubleSpinBox_2->setValue(MAX_QUANT_CNTRL_OBJ +0.001);
 }
}

void MainWindow::on_doubleSpinBox_2_editingFinished()
{
    if(ui->checkBox_4->isChecked())
    {
        on_checkBox_4_stateChanged(2);
    }


    int numCtrlOld = mNumCtrlObj;
    mNumCtrlObj = (int)(ui->doubleSpinBox_2->value() + 0.1);
    if(mNumCtrlObj< numCtrlOld)
    {
        numCtrlOld = mNumCtrlObj;
    }


    for (int i=0; i<  ui->tableWidget_3->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget_3->columnCount(); j++)
                 {
                     QTableWidgetItem* ptwi0 = new QTableWidgetItem;
                     ptwi0 ->setText(QString::number(0.));
                     ui->tableWidget_3->setItem(i,j,ptwi0);
                 }
    for (int i = 0; i < numCtrlOld; ++i)
        for (int j =0; j < ui->tableWidget_3->columnCount(); ++j)
        {
         double coeff = (j == 1)?1:1./ M_PI * 180.;
         ui->tableWidget_3->item(i,j)->setText(QString::number(mparrDataCtrlObj[i * 5 +j] * coeff));
        }
   inputData();
}

void MainWindow::on_doubleSpinBox_3_valueChanged(double arg1)
{
    ui->checkBox_4->setChecked(2);
    if (ui->checkBox_4->isChecked())
    {
      on_checkBox_4_stateChanged(2);

    }
}

void MainWindow::on_doubleSpinBox_4_valueChanged(const QString &arg1)
{
    ui->checkBox_4->setChecked(2);
    if (ui->checkBox_4->isChecked())
    {
      on_checkBox_4_stateChanged(2);
    }
}

//---------------------------------------------



void MainWindow::on_checkBox_stateChanged(int arg1)
{

}



void MainWindow::on_checkBox_5_stateChanged(int arg1)
{

}

void MainWindow::on_checkBox_2_stateChanged(int arg1)
{

}
/*
void MainWindow::on_pushButton_2_clicked()
{
    const int N=1000;
     double valAver = 0., valAverSquare = 0.;
    for (int i =0; i < N; ++i)
    {
        ui->pushButton->click();
        CalcStatParams(marrXEst_IU_PSK[4],i +1 , valAver, valAverSquare);


        int gg=0;
    }
    double disp0 = valAverSquare - valAver*valAver;
    double disp = marrMtrxK[QUANT_GDG_PRMS * 4  +4];
    int iui=0;

}
*/
//-----------------------------------
bool MainWindow::checkObj()
{
    for (int i =0; i < mNumCtrlObj; ++i)
    {
        if (mparrDataCtrlObj[i * 5 + 1]< 10.)
        {
            return false;
        }
    }
    return true;
}

