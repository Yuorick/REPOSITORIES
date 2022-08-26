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
#include "ImitMod.h"
#include "Rls_Usbl2D.h"
#include "Rls_Usbl3D.h"
#include "Gps.h"
#include "PosSolver.h"
#include "LblSolver.h"
#include "Usbl_2D_Solver.h"
#include "Usbl_3D_Solver.h"


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
    wchar_t mwchOutPutFold[400];

    // файл спрофилем скорости
    wchar_t mpwchPrflFIle[400];

    TTable_1D mtblEstPrfl;

    TTable_1D mtblRealPrfl;

    // радиус зоны приема
    double mZonaR;

    // макс глубина
    double mDeepthMax;

    TEnvironment mEnvironment;

    // вектор параметров позиционироваия ГРЛС истинный
   double marrTruePosParams[6];
   // вектор параметров позиционироваия ГРЛС априорный
    double marrAprioriPosParams[6];
    // вектор параметров позиционироваия ГРЛС оценка
    double marrEstPosParams[6];

     // координаты маяка  истинные
    double marrTrueHeadLightPos[3];
    // координаты маяка  априорные
  double marrAprioriHeadLightPos[3];
  // координаты маяка  оцененные
   double marrEstHeadLightPos[3];

   // вектор параметров позиционироваия GPS истинный
  double marrGPSTruePosParams[6];
  // вектор параметров позиционироваия GPS априорный
   double marrGPSAprioriPosParams[6];
  QGps mGPS;

    //корабль

    double mvalMaxQ;
    double mvalMaxTet;
    double mvalMaxPsi;


    // корабль
    QPeaceVess mVess;

    // СИНС
   QPeaceSins mSins;

   // платформа
   QPlatform mPlatform;

   // РЛС
   QHidroRLS mHidroRLS;

   QRls_Usbl2D mRls_Usbl2D;

   QRls_Usbl3D mRls_Usbl3D;

   QHidroRLS* mpHidroRLS;

   enumTypeOfRLS mTypeOfRLS;

    // темп выдачи информации СИНС
    double mTimeTempSins;

    // время начала моделирования
    double mT0;

    //тип траектории корабля
    enumTypeOfVessTraj mTypeOfVessTraj;

    // к-во измерений
    int mQuantMeas;

    QVector<QBigMeasure> mvctMeasures;

    // имит модель
    QImitMod mImitMod;

    bool mbtableWidget_4Init;
    bool mbtableWidget_6Init;
    bool mbtableWidget_8Init;

    //массив с текущей информацией выходной таблицы по маяку
    double marrMayjak_TblData[9];
    //массив с текущей информацией выходной таблицы по антенне
    double marrGRLS_TblData[18];

    // время обработки сигнала на маяке
    double mTObrabotki;

    QPosSolver *mpPosSolver;

    QLblSolver mLblSolver;
    QUsbl_2D_Solver mUsbl_2D_Solver;
    QUsbl_3D_Solver mUsbl_3D_Solver;

    void inputData();

private slots:

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

   void on_comboBox_4_currentIndexChanged(int index);

   void on_tableWidget_4_cellChanged(int row, int column);

   void read_tableWidget_4(double *arr);

   void write_tableWidget_4(double *arr);

   void on_tableWidget_6_cellChanged(int row, int column);

   void read_tableWidget_6(double *arr);

   void write_tableWidget_6(double *arr);

   void read_tableWidget_7(double *arr);

   void write_tableWidget_7(double *arr);

   void read_tableWidget_5(double *arr);

   void write_tableWidget_5(double *arr);

   void read_tableWidget_6();


   //int createInputDataReport(wchar_t*FileName, const bool bHeader, QTargTrackingOperation &TargTrackingOperation);
//
 //  int createOutputDataReport(wchar_t*FileName, const bool bHeader
        //                      ,const double* arrHorSyst,const double* arrHorDisp
        //                      ,const double*  arrVertSyst,const double* arrVertDisp);

   void on_pushButton_5_clicked();

   void on_tableWidget_7_itemActivated(QTableWidgetItem *item);

   void on_tableWidget_7_itemSelectionChanged();


   void on_tableWidget_5_itemActivated(QTableWidgetItem *item);

   void on_tableWidget_5_itemSelectionChanged();

   void read_tableWidget_8(double *arr1,double *arr2);

  // void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
