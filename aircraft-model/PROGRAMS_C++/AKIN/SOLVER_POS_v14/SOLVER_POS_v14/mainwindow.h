#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>

#include "Comp.h"
#include <QTableWidget>

#include "Environment.h"
#include "Table_1D.h"
#include "PeaceVess.h"
#include "PeaceSins.h"
#include "Platform.h"
#include "HidroRLS.h"

#include "Rls_Usbl2D.h"
#include "Rls_Usbl3D.h"
#include "Gps.h"
#include "PosSolver.h"
#include "LblSolver.h"

#include "SubWaterBeam.h"

#include "mythread.h"
#include "PntCloud.h"

class QDoubleSpinBox;




class QMeasurmentImitator;


class QWarShip;
class QWarSins;



class TEnvironment;
class TTable_1D;
class QPeaceVess;
class QPosSolver;
class QLblSolver;
class QUsbl_2D_Solver;
class QUsbl_3D_Solver;


#define nullptr 0

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow();


    // папка с графиками результатов
  //  wchar_t mwchOutPutFold[400];

    //
    QString mqstrOutPutFold;

    // файл спрофилем скорости
   // wchar_t mpwchPrflFIle[400];
    QString mqstrPrflFIle;

    // файл с DATA
   // wchar_t mpwchDataFIle[400];
    QString mqstrDataFIle;

    TTable_1D mtblEstPrfl;



    // радиус зоны приема
    double mZonaR;
// радиус рабочей зоны
    double mWorkR;

   // вектор параметров позиционироваия ГРЛС априорный
    double marrAntPosParams[6];



    // координаты маяка  априорные
  double marrBeaconPos[3];


  // вектор параметров позиционироваия GPS априорный
   double marrGPSPosParams[6];
  //
   double marrXRez [8];



   TYPE_OF_SOLVER_TASK mTypeOfSolverTask;

    // к-во измерений выкачанное изфайла
   int mQuantMeas0;

    QVector<QBigMeasure> mvctMeasures_work;

    bool mbtableWidget_4Init;
    bool mbtableWidget_6Init;
    bool mbtableWidget_8Init;

    bool mbtableWidget_5Init;
    bool mbtableWidget_7Init;

    bool mbOpenDialog_000;
    bool mbOpenDialog_acd;

    // время обработки сигнала на маяке
    double mTObrabotki;

    QPosSolver *mpPosSolver;

    QLblSolver mLblSolver;

    enumTypeOf_000 mType_of_000;

    //порядковый номер,решаемой задачи,нужен для ползунка
    int mNumOfTask;

    int mMaxQuantIter_LBL;

    int mMaxQuantIter_USBL;

    int mNum_LBL;
    double merrror_LBL;
    double mNeviaz_LBL;

    int mNum_USBL;
    double merrror_USBL;
    double mNeviaz_USBL;

    // минимальный угол места замера
    double mTetta_min;

    // максимальныйугол места замера
    double mTetta_max;

    //  ворота для перебора по углу
    double mGape;

    //  шаг перебора по углу
    double mAngStep;

    ALGOR_TYPE mALGOR_TYPE;

    double mLBLcoeff;

    double mUSBLcoeff;

    TYPE_OF_OUTPUT_FILE mTYPE_OF_OUTPUT_FILE;

    // номер кликания на кнопку Старт
    int mNumCkickStart;


    void inputData();

    //******** ЭТО ДЛЯ ПОЛЗУНКА *******************
private slots:
    void on_pushButtonStartTask_clicked();

    void on_pushButtonStopTask_clicked();

    void onTaskProgress(int numIter, QVector<double>vctNevSquare, QVector<double> vctMean
                        , QVector<double> vctDisp, int numpart ,bool bEndPart, QVector<double>vctX0);

    void onKleaningEnded(int numpart, QVector<int>vctGoodZamersNum, QVector<double>vctMean, QVector<double>Disp);

    void onTaskFinished();

    void onCalcFinished(int result);
//********  !ЭТО ДЛЯ ПОЛЗУНКА **************** !

private slots:

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();



    void on_tableWidget_6_cellChanged(int row, int column);

    void on_tableWidget_4_cellChanged(int row, int column);



    void on_pushButton_clicked();

    void on_pushButton_4_clicked();

    void on_comboBox_5_currentIndexChanged(int index);

    void on_comboBox_5_activated(const QString &arg1);

    void on_tableWidget_6_cellClicked(int row, int column);

    void read_tableWidget_6();

    void on_comboBox_9_currentIndexChanged(int index);

    void on_pushButton_5_clicked();





    void on_comboBox_6_activated(const QString &arg1);

    void on_comboBox_6_currentIndexChanged(int index);

    void on_comboBox_6_currentIndexChanged(const QString &arg1);

    void on_comboBox_currentIndexChanged(const QString &arg1);

    void on_comboBox_9_currentTextChanged(const QString &arg1);

    void on_comboBox_activated(const QString &arg1);

    void on_pushButton_6_clicked();

    void on_pushButton_NEW_clicked();



    void on_labelTaskInfo1_objectNameChanged(const QString &objectName);

private:
   void read_tableWidget_7(double *arr);

   void write_tableWidget_7(double *arr);

   void read_tableWidget_5(double *arr);

   void write_tableWidget_5(double *arr);

   void read_tableWidget_6(double *arrAprioriBeaconPos);


   //int createInputDataReport(wchar_t*FileName, const bool bHeader, QTargTrackingOperation &TargTrackingOperation);
//
 //  int createOutputDataReport(wchar_t*FileName, const bool bHeader
        //                      ,const double* arrHorSyst,const double* arrHorDisp
        //                      ,const double*  arrVertSyst,const double* arrVertDisp);



   void read_tableWidget_8(double *arr2);  



   void readHeaderDataFile(QString TargFileName,int &iN, double *arrP_Ant
                                       , double *arrP_Gps, double *arrP_Beacon, double &valTBeacon);

   static void createPanoramePictXY(wchar_t *wchOutPutFold, double *arrBeaconPos, double *arrAntPosParams, const QBigMeasure &Meas);

 //  void readDataFile_(QString TargFileName,const int QuantMeas
  //                                               ,QBigMeasure *parrMeas);

   //void readDataFile(QString TargFileName,const int QuantMeas
      //                                           ,QBigMeasure *parrMeas);

   void fill_tableWidget_4(double *arr);

   void fill_tableWidget_6(double *arr);

   void fill_tableWidget_8(double *arr);

   void fillFields();

   void on_comboBox_7_currentIndexChanged(int index);

   void write_tableWidget_6(double *arr);

   //bool estimateParams__( QPosSolver *pQPosSolver, double *arrX0,double *arrXRez, double &valNeviaz0);

   void setTable(QTableWidget& tbl, double *arr);

   void setTableZero(QTableWidget& tbl);

   void tuneGroupBoxUSBL(int i);

   void cleanInteface();

   QString createInputParamsConstructorString();

   QString createInputParamsString(int num);
private:
    Ui::MainWindow *ui;
    MyThread *myThread;
};

#endif // MAINWINDOW_H
