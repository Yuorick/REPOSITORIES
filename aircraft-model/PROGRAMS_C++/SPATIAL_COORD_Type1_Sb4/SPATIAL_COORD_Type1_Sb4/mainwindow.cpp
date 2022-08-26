#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include "ContObj6.h"
#include "MatrixProccess.h"
#include "Gauss.h"
#include <QErrorMessage>

extern const double COeff;
extern  bool BEZ_SHUMOV = true;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{    
    ui->setupUi(this);
    mbFirstTime = true;
    setupFirst();
    fillWindow();
    mbFirstTime = false;
    mbCalculationIsDone = false;
    ui->frame_7->setHidden(true);
   //

}
void MainWindow:: inputData()
{
    mVessSKZ = 0.;
    mCntrlObjSKZ = 0.;


    //0 калькулятор установка
    mbCalculationIsDone = false;
    ui->frame_7->setHidden(true);
    fillCalculatorArraysPrevious();
    fillCalculatorFrame();
    ///

    // 1
    for (int i = 0; i < QUANT_GDG_PRMS; ++i)
    {
     marrXTrue_RLK_PSK[i] =  ui->tableWidget->item(0,i)->text().toDouble();
     marrXZv_RLK_PSK[i] =  ui->tableWidget->item(1,i)->text().toDouble();
    }

    // 2
    for (int i = 0; i < QUANT_GDG_PRMS; ++i)
    {
     marrXTrue_IU_PSK[i] =  ui->tableWidget_2->item(0,i)->text().toDouble();
     marrXZv_IU_PSK[i] =  ui->tableWidget_2->item(1,i)->text().toDouble();
    }



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
    mAngleSigIU =  0.001;

    // СКЗ ошибки дальномера
    mDistSigIU = 0.1;

    //признак случайных шумов измерений
    mbNoise= ui->checkBox_2->isChecked();

    // Коррел матрица ошибок первичных измерений РЛК
    memset(marrCorMtrxRLK, 0, 9 * sizeof(double));
    marrCorMtrxRLK[0] = marrCorMtrxRLK[8] = (0.0005)*(0.0005);
    marrCorMtrxRLK[4] = 100.;

    marrDispEilerCntrKP[0] = 0.0006 * 0.0006;
    marrDispEilerCntrKP[1] = 0.0004 * 0.0004;
    marrDispEilerCntrKP[2] = 0.0004 * 0.0004;



    //признак случайных шумов измерений
    mbNoise= ui->checkBox_2->isChecked();

    // Коррел матрица ошибок первичных измерений РЛК
    memset(marrCorMtrxRLK, 0, 9 * sizeof(double));
    marrCorMtrxRLK[0] = marrCorMtrxRLK[8] = (0.0005)*(0.0005);
    marrCorMtrxRLK[4] = 100.;

    marrDispEilerCntrKP[0] = 0.0006 * 0.0006;
    marrDispEilerCntrKP[1] = 0.0004 * 0.0004;
    marrDispEilerCntrKP[2] = 0.0004 * 0.0004;



    //fill_Partly_OutTable();

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
    //fill_Partly_OutTable();
    // подготовка массива КО
    QContObj6 *pqtrlObjectArr = new QContObj6[mNumCtrlObj];
    for (int i = 0; i < mNumCtrlObj; ++i)
    {
        double arrS_KGSK[3] = {0.};
        double *p = &(mparrDataCtrlObj[5 * i]);
        arrS_KGSK[0] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * sin(p[0]);
        arrS_KGSK[1] = sqrt(p[1]* p[1] - mCtrlObjHeight*mCtrlObjHeight) * cos(p[0]);
        arrS_KGSK[2] = mCtrlObjHeight;

        pqtrlObjectArr[i] =  QContObj6(&(p[2]), marrDispEilerCntrKP
                                       ,marrXTrue_RLK_PSK,marrXTrue_IU_PSK
                                       ,marrXZv_RLK_PSK,marrXZv_IU_PSK
                                       ,arrS_KGSK, marrCorMtrxRLK
                                       , mAngleSigIU, mDistSigIU,mVessSKZ,mCntrlObjSKZ);
    }


    double rez =-1.;
    if (ui->checkBox_3->isChecked())
    {
         rez = QContObj6::imitate_and_estimateAngles_Var1_(pqtrlObjectArr,mNumCtrlObj
                                            , mbNoise, marrXEst, marrMtrxK);
    }
    else
    {
        rez = QContObj6::imitate_and_estimateX_Var1(pqtrlObjectArr,mNumCtrlObj
               , mbNoise,  marrXEst, marrMtrxK);
    }

   double arrEilersRlkIdeal[3] = {0.};
    QContObj6::calc_IdealRelativeEilers_Var1(&(marrXTrue_RLK_PSK[3]),&(marrXTrue_IU_PSK[3])
                   ,&(marrXZv_IU_PSK[3]) ,&(marrXRLK_RelativeTrue[3]));
    MtrxSumMatrx(marrXTrue_RLK_PSK, marrXZv_IU_PSK,1, 3, marrXRLK_RelativeTrue);


        memcpy(marrXEst_RLK,marrXEst, QUANT_GDG_PRMS * sizeof(double));
        MtrxSumMatrx(marrXEst, marrXZv_IU_PSK,1, 3,marrXEst_RLK) ;

    int ia = 0;
    double temp = 0.;
    for (int j = 0; j < QUANT_GDG_PRMS; ++j)
    {
        QTableWidgetItem* item0 = new QTableWidgetItem;
        if ( j <4)
        {
            ia = (10000. *(marrXEst_RLK[j]));
            temp = ia/10000.;
        item0 ->setText(QString::number(temp));
        ui->tableWidget_4->setItem(0,j,item0);
        }
        else
        {
            ia = (10000. *marrXEst_RLK[j]);
            temp = ia/10000.;
            item0 ->setText(QString::number(temp));
            ui->tableWidget_4->setItem(0,j,item0);
        }


        QTableWidgetItem* item3 = new QTableWidgetItem;
        if (j < 3)
        {
            ia = (10000. *sqrt(marrMtrxK[j * QUANT_GDG_PRMS + j]));
            temp = ia/10000.;
        item3 ->setText(QString::number(temp));
        }
        else
        {
            ia = (10000. *sqrt(marrMtrxK[j * QUANT_GDG_PRMS + j])* 1000.);
            temp = ia/10000.;
         item3 ->setText(QString::number(temp));
        }
        ui->tableWidget_4->setItem(1,j,item3);


        //
        QTableWidgetItem* item4 = new QTableWidgetItem;
        if ( j <4)
        {
            ia = (10000. *(marrXRLK_RelativeTrue[j]));
            temp = ia/10000.;
        item4 ->setText(QString::number(temp));
        ui->tableWidget_4->setItem(2,j,item4);
        }
        else
        {
            ia = (10000. *marrXRLK_RelativeTrue[j]);
            temp = ia/10000.;
            item4 ->setText(QString::number(temp));
            ui->tableWidget_4->setItem(2,j,item4);
        }
    }

    delete []pqtrlObjectArr;



    // ТЕСТ
    double marrDeckAngles[3] = {0.};

    double marrBaring[3] ={0.78, 6000., 50.};
    double arrTargS_KGSK[3] = {0.};
    arrTargS_KGSK[0] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * sin(marrBaring[0]);
    arrTargS_KGSK[1] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * cos(marrBaring[0]);
    arrTargS_KGSK[2] = marrBaring[2];
    ///

    // 2
  QContObj6  qtrlObjectZv(marrDeckAngles, marrDispEilerCntrKP
                                   ,marrXTrue_RLK_PSK,marrXTrue_IU_PSK
                                   ,marrXZv_RLK_PSK,marrXZv_IU_PSK
                                   ,arrTargS_KGSK, marrCorMtrxRLK
                                   , mAngleSigIU, mDistSigIU, mVessSKZ, mCntrlObjSKZ);

    double arrVTrue_RLK[3] = {0.};
    QContObj6::recalcPositionFromKGSK_to_SphericalSK(marrDeckAngles,qtrlObjectZv.marrS_KGSK,marrXTrue_RLK_PSK,arrVTrue_RLK);

    double arrVTrue_IU[3] = {0.};
    QContObj6::recalcPositionFromKGSK_to_SphericalSK(marrDeckAngles,qtrlObjectZv.marrS_KGSK,marrXTrue_IU_PSK,arrVTrue_IU);

    ///
    double arrTrueViu[3] = {0.};
    const bool BDalnomer = false;
    double marrCntrlAngTrue[2] = {0.};
    double arrXTrue_RLK_PSK[QUANT_GDG_PRMS] = {0.};
    memcpy(arrXTrue_RLK_PSK, marrXTrue_RLK_PSK, sizeof(double )*QUANT_GDG_PRMS);
    MtrxMinusMatrx(arrXTrue_RLK_PSK, marrXTrue_IU_PSK,1, 3, arrXTrue_RLK_PSK);
    QContObj6::calc_Viu( arrXTrue_RLK_PSK,marrXTrue_IU_PSK, arrVTrue_RLK, marrCntrlAngTrue);
    ///

   // 3
   // double arrZvViu[3] = {0.};
   // QContObj6::calc_Viu( BDalnomer,arrXZv,marrXZv_IU_PSK, arrVTrue_RLK, marrCntrlAngZv);
    ///

    // 4

     QContObj6::calc_Viu( marrXEst,marrXZv_IU_PSK, arrVTrue_RLK, marrCntrlAngEst);
     ///
     /// \brief uu
     ///
     // 5

     QContObj6::calc_Viu( marrXZv_RLK_PSK,marrXZv_IU_PSK, arrVTrue_RLK, marrCntrlAngPrevEst);
     int uu =0;

   /*  double arrDifferZv[3] = {0.};
     MtrxMinusMatrx(arrTrueViu, arrZvViu,1, 3, arrDifferZv);

     double arrDifferEst[3] = {0.};
     MtrxMinusMatrx(arrTrueViu, arrEstViu,1, 3, arrDifferEst);


     // 5
     double arrC[2 * NUmParams] = {0.};
     qtrlObjectZv.calc_dViu_po_dX_Var1_(false,arrXTrue,marrXTrue_IU_PSK
                                , arrVTrue_RLK, marr_dVIU_po_dX);

     double valScal = 1000.;

     MatrxMultScalar(marrCntrlAngTrue, 2, 1, valScal, marrCntrlAngTrue);
     MatrxMultScalar(marrCntrlAngZv, 2, 1, valScal, marrCntrlAngZv);
     MatrxMultScalar(marrCntrlAngEst, 2, 1, valScal, marrCntrlAngEst);


     MatrxMultScalar(marr_dVIU_po_dX, 3, 1, valScal, marr_dVIU_po_dX);
     MatrxMultScalar(&(marr_dVIU_po_dX[6]), 3, 1, valScal, &(marr_dVIU_po_dX[6]));
    */

     mbCalculationIsDone = true;
     ui->frame_7->setHidden(false);
     //fillCalculatorFrame();
     fillCalculatorArraysPrevious();

     fillCalculatorFrame();
}

/*
void MainWindow::fill_Partly_OutTable()
{


    MtrxMinusMatrx(marrXTrue_RLK_PSK, marrXTrue_IU_PSK,1, QUANT_GDG_PRMS, marrXRLK_RelativeTrue);
    MtrxMinusMatrx(marrXZv_RLK_PSK, marrXZv_IU_PSK,1, QUANT_GDG_PRMS, marrXRLK_RelativeZv);
    double arr_delta[QUANT_GDG_PRMS ] = {0.};
    MtrxMinusMatrx(marrXRLK_RelativeZv, marrXRLK_RelativeTrue,1, QUANT_GDG_PRMS, arr_delta);

   int ia = 0;
   double temp = 0.;
    for (int j = 0; j < QUANT_GDG_PRMS; ++j)
    {
        QTableWidgetItem* item0 = new QTableWidgetItem;
        ia = 10000. *marrXRLK_RelativeTrue[j];
        temp = ia/10000.;
        item0 ->setText(QString::number(temp));
        ui->tableWidget_4->setItem(0,j,item0);

        QTableWidgetItem* item1 = new QTableWidgetItem;
        ia = 10000. *marrXRLK_RelativeZv[j];
        temp = ia/10000.;
        item1 ->setText(QString::number(temp));
        ui->tableWidget_4->setItem(1,j,item1);

        if(j <3)
        {
        QTableWidgetItem* item2 = new QTableWidgetItem;
        ia = 10000. *arr_delta[j];
        temp = ia/10000.;
        item2 ->setText(QString::number(temp));
        ui->tableWidget_4->setItem(3,j,item2);
        }
        else
        {
            QTableWidgetItem* item2 = new QTableWidgetItem;
            ia = 10000. *1000. *(arr_delta[j]);
            temp = ia/10000.;
            item2 ->setText(QString::number(temp));
            ui->tableWidget_4->setItem(3,j,item2);
        }

    }

}

*/
MainWindow::~MainWindow()
{
    delete ui;
}
void MainWindow::setupFirst()
{

    // вектор истинных параметров позиционирования РЛК в ПСК
    marrXTrue_RLK_PSK[0] = 0.;
    marrXTrue_RLK_PSK[1] = 0;
    marrXTrue_RLK_PSK[2] = 0.;
    marrXTrue_RLK_PSK[3] = M_PI/4.+0.003;
    marrXTrue_RLK_PSK[4] =10./ 180. * M_PI + 0.003;
    marrXTrue_RLK_PSK[5] = 0.003;


   // вектор первичных оценок параметров позиционирования РЛК в ПСК
    //marrXZv_RLK_PSK[0] = 1.;
    //marrXZv_RLK_PSK[1] = 10.;
    //marrXZv_RLK_PSK[2] = 25.;
    marrXZv_RLK_PSK[0] = 0.;
    marrXZv_RLK_PSK[1] = 0.;
    marrXZv_RLK_PSK[2] = 0.;
    marrXZv_RLK_PSK[3] = M_PI/4.;
    marrXZv_RLK_PSK[4] = 10./ 180. * M_PI;
    marrXZv_RLK_PSK[5] = 0.;

    // вектор истинных параметров позиционирования ИУ в ПСК
    marrXTrue_IU_PSK[0] = 0.;
    marrXTrue_IU_PSK[1] = 0.;
    marrXTrue_IU_PSK[2] = 0.;
    marrXTrue_IU_PSK[3] = - M_PI/4.-0.003;
    marrXTrue_IU_PSK[4] = 10./ 180. * M_PI - 0.003;
    marrXTrue_IU_PSK[5] = -0.003;


   // вектор первичных оценок параметров позиционирования ИУ в ПСК
    //marrXZv_IU_PSK[0] = 5.;
    //marrXZv_IU_PSK[1] = 16.;
    //marrXZv_IU_PSK[2] = 4.;
    marrXZv_IU_PSK[0] = 0.;
    marrXZv_IU_PSK[1] = 0.;
    marrXZv_IU_PSK[2] = 0.;
    marrXZv_IU_PSK[3] = - M_PI/4.;
    marrXZv_IU_PSK[4] = 10./ 180. * M_PI;
    marrXZv_IU_PSK[0]= 0.;


    // корреляционная матрица ошибок оценивания
    memset(marrMtrxK, 0, QUANT_GDG_PRMS * QUANT_GDG_PRMS * sizeof(double));

    // Кол-во контрольных объектов (КО)
    mNumCtrlObj = 4;





    //признак случайных шумов измерений
    mbNoise = false;

    // массив исходных данных КО из входной таблицы
    memset(mparrDataCtrlObj, 0, MAX_QUANT_CNTRL_OBJ *5 * sizeof(double));
    mparrDataCtrlObj[0] = 0.;
    mparrDataCtrlObj[1] = 6000.;
    mparrDataCtrlObj[5] = M_PI/3;
    mparrDataCtrlObj[6] = 6000.;
    mparrDataCtrlObj[11] = 8000.;
    mparrDataCtrlObj[15] = M_PI/2.;
    mparrDataCtrlObj[16] = 10000.;

    // высота КО
    mCtrlObjHeight = 50.;

}
//-----------------------------------------------------

void MainWindow::fillWindow()
{
  int ia = 0;
  double temp = 0.;
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
         ia = 10000.*mparrDataCtrlObj[i * 5 +j] * coeff;
         temp = ia/10000.;
         ui->tableWidget_3->item(i,j)->setText(QString::number(temp));
        }
    //2

   for (int i=0; i<  ui->tableWidget->rowCount() ; i++)
                 for (int j = 0; j < ui->tableWidget->columnCount(); j++)
                 {
                     QTableWidgetItem* item = new QTableWidgetItem;
                     item->setText(QString::number(0.));

                     ui->tableWidget->setItem(i,j,item);
                 }

    for (int j =0; j < ui->tableWidget->columnCount(); ++j)
    {
        double coeff = (j == 1)?1:1./ M_PI * 180.;
        ia = 10000.*marrXTrue_RLK_PSK[j];
        temp = ia/10000.;
     ui->tableWidget->item(0,j)->setText(QString::number(temp));

     ia = 10000.*marrXZv_RLK_PSK[j];
     temp = ia/10000.;
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
        ia = 10000.*marrXTrue_IU_PSK[j];
        temp = ia/10000.;
     ui->tableWidget_2->item(0,j)->setText(QString::number(temp));

     ia = 10000.*marrXZv_IU_PSK[j];
     temp = ia/10000.;
     ui->tableWidget_2->item(1,j)->setText(QString::number(temp));
    }
    ///

        //3
        ui->doubleSpinBox_2->setValue(mNumCtrlObj);
        ///

        //4
        ui->doubleSpinBox->setValue(mCtrlObjHeight);


        //6
        ui->checkBox_2->setChecked(false);

}




void MainWindow::on_tableWidget_2_cellChanged(int row, int column)
{
    //if(!mbFirstTime)
   // {
    //inputData();
    //fill_Partly_OutTable();
   // }
}

void MainWindow::on_tableWidget_cellChanged(int row, int column)
{
  //  if(!mbFirstTime)
   // {
   // inputData();
   // }
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

        // угол килевой качки
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



void MainWindow::on_doubleSpinBox_2_editingFinished()
{
   // if(ui->checkBox_4->isChecked())
    //{
    //    on_checkBox_4_stateChanged(2);
   // }


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
//-----------------------------------------
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
    memset(marrCntrlAngTrue, 0,2 * sizeof(double));

    // массив углов наведения оценка первоначальная - Betta, Eps
  //  memset(marrCntrlAngZv, 0,2 * sizeof(double));

    // массив углов наведения оценка уточненная - Betta, Eps
    memset(marrCntrlAngEst, 0,2 * sizeof(double));

    // матрица частных производных углов наведения по параметрам позиционирования
    memset(marr_dVIU_po_dX, 0,12 * sizeof(double));
}
//---------------------------------------------
//заполнение фрейма калькулятора
void MainWindow::fillCalculatorFrame()
{
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
            //3  таблица углов наведения вычисленная
          int ia = 0.;
          double temp = 0.;
          for (int i =0; i < 4; ++i)
            for (int j = 0; j < ui->tableWidget_7->columnCount(); j++)
            {
            QTableWidgetItem* ptwi0 = new QTableWidgetItem;
            ptwi0 ->setText(QString::number(0.));
            ui->tableWidget_7->setItem(i,j,ptwi0);
            }

            for (int j =0; j < ui->tableWidget_7->columnCount(); ++j)
            {
                 int ia = marrCntrlAngTrue[j] * 100.;
                 double temp = ((double)ia)/100.;
            ui->tableWidget_7->item(0,j)->setText(QString::number(temp ));
            }

            for (int j =0; j < ui->tableWidget_7->columnCount(); ++j)
            {
                int ia = marrCntrlAngEst[j] * 100.;
                double temp = ((double)ia)/100.;
            ui->tableWidget_7->item(1,j)->setText(QString::number(temp ));
            }

            for (int j =0; j < ui->tableWidget_7->columnCount(); ++j)
            {
                ia = int(100. *(marrCntrlAngEst[j]-marrCntrlAngTrue[j] ));
                temp = ia / 100.;
            ui->tableWidget_7->item(2,j)->setText(QString::number(temp ));
           }

            for (int j =0; j < ui->tableWidget_7->columnCount(); ++j)
            {
                int ia = marrCntrlAngPrevEst[j] * 100.;
                double temp = ((double)ia)/100.;
            ui->tableWidget_7->item(3,j)->setText(QString::number(temp ));
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
                int ia = ((int)(marr_dVIU_po_dX[i * QUANT_GDG_PRMS + j] * 1000.));
                double temp = ((double)ia)/1000.;
                //ui->tableWidget_8->item(i,j)->setText(QString::number(marr_dVIU_po_dX[i * 5 + j] ));
                ui->tableWidget_8->item(i,j)->setText(QString::number(temp ));
            }


}


void MainWindow::on_pushButton_3_clicked()
{
    calculatorInput();
    // 1.
    double arrTargS_KGSK[3] = {0.};
    arrTargS_KGSK[0] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * sin(marrBaring[0]);
    arrTargS_KGSK[1] = sqrt(marrBaring[1]* marrBaring[1] - marrBaring[2] * marrBaring[2]) * cos(marrBaring[0]);
    arrTargS_KGSK[2] = marrBaring[2];
    ///



    // 2
  QContObj6  qtrlObjectZv(marrDeckAngles, marrDispEilerCntrKP
                                   ,marrXTrue_RLK_PSK,marrXTrue_IU_PSK
                                   ,marrXZv_RLK_PSK,marrXZv_IU_PSK
                                   ,arrTargS_KGSK, marrCorMtrxRLK
                                   , mAngleSigIU, mDistSigIU, mVessSKZ, mCntrlObjSKZ);

    double arrVTrue_RLK[3] = {0.};
    QContObj6::recalcPositionFromKGSK_to_SphericalSK(marrDeckAngles,qtrlObjectZv.marrS_KGSK,marrXTrue_RLK_PSK,arrVTrue_RLK);

    double arrVTrue_IU[3] = {0.};
    QContObj6::recalcPositionFromKGSK_to_SphericalSK(marrDeckAngles,qtrlObjectZv.marrS_KGSK,marrXTrue_IU_PSK,arrVTrue_IU);

    ///


    double arrXTrue_RLK_PSK[QUANT_GDG_PRMS] = {0.};
    memcpy(arrXTrue_RLK_PSK, marrXTrue_RLK_PSK, sizeof(double )*QUANT_GDG_PRMS);
    MtrxMinusMatrx(arrXTrue_RLK_PSK, marrXTrue_IU_PSK,1, 3, arrXTrue_RLK_PSK);
    QContObj6::calc_Viu( arrXTrue_RLK_PSK,marrXTrue_IU_PSK, arrVTrue_RLK, marrCntrlAngTrue);
    ///



    // 4

     QContObj6::calc_Viu( marrXEst,marrXZv_IU_PSK, arrVTrue_RLK, marrCntrlAngEst);
     ///
  // 5
      QContObj6::calc_Viu( marrXZv_RLK_PSK,marrXZv_IU_PSK, arrVTrue_RLK, marrCntrlAngPrevEst);

     qtrlObjectZv.calc_dViu_po_dX_Var1_(marrXEst,marrXTrue_IU_PSK
                                , arrVTrue_RLK, marr_dVIU_po_dX);

     double valScal = 1000.;

     MatrxMultScalar(marrCntrlAngTrue, 2, 1, valScal, marrCntrlAngTrue);
     MatrxMultScalar(marrCntrlAngPrevEst, 2, 1, valScal,marrCntrlAngPrevEst);
     MatrxMultScalar(marrCntrlAngEst, 2, 1, valScal, marrCntrlAngEst);


     MatrxMultScalar(marr_dVIU_po_dX, 3, 1, valScal/COeff, marr_dVIU_po_dX);
     MatrxMultScalar(&(marr_dVIU_po_dX[QUANT_GDG_PRMS]), 3, 1, valScal/COeff, &(marr_dVIU_po_dX[QUANT_GDG_PRMS]));
     int uu =0;

     fillCalculatorFrame();

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



void MainWindow::on_doubleSpinBox_2_valueChanged(double arg1)
{
    if (arg1 > MAX_QUANT_CNTRL_OBJ)
    {
      ui->doubleSpinBox_2->setValue(MAX_QUANT_CNTRL_OBJ +0.001);
    }
}

/*
void MainWindow::on_pushButton_2_clicked()
{

    const int N=1000;
     double valAver = 0., valAverSquare = 0.;
    for (int i =0; i < N; ++i)
    {
        ui->pushButton->click();
        CalcStatParams(marrXEst[2],i +1 , valAver, valAverSquare);
        int gg=0;
    }
    double disp0 = valAverSquare - valAver*valAver;
    double disp = marrMtrxK[QUANT_GDG_PRMS *2 + 2];
    int iui=0;

}
*/
//-----------------------------------
bool MainWindow::checkObj()
{
    for (int i =0; i < mNumCtrlObj; ++i)
    {
        if (mparrDataCtrlObj[i * 5 +1]< 10.)
        {
            return false;
        }
    }
    return true;
}
