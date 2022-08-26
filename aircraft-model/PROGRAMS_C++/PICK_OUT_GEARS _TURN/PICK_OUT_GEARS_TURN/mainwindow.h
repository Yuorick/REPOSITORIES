#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QApplication>
#include <QMainWindow>
#include <QDoubleSpinBox>
#include "Helic.h"
//#include "CommonGliderElement.h"
//#include "RuleGliderElement.h"
//#include "Rotor.h"
//#include "Blade.h"
#include "Environment.h"
//#include "OperationNew.h"
//#include "FlyTask.h"
#include "TurnMove.h"
#include "LineMove.h"
#include "Rotating.h"
#include "Hover.h"


#define QUANT_FLY_TASK_POINTS 50 // макс число точек полетного задания



class QDoubleSpinBox;
class THelic;
//class TCommonGliderElement;
//class TRuleGliderElement;
//class TRotor;
//class TBlade;
class TEnvironment;

//class TFlyTask;
class TTurnMove;
class TLineMove;
class TRotating;
class THover;


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow();

    // полетное время
    long double mTFly ;


    // параметры лемнискаты:
        long double mLemn_a; // корень из произведения расстояний

        long double mLemn_c; // половина расстояния между фокусами
        long double mLemn_l; // запас устойчивости

   THelic mHelic;

   // Re MIN
     double mReMin;

    // Re Max
     double mReMax;
     // TgMax
     double mTgMax;

    long double mvalYaw ;
    long double  mvalRad ;
    long double mvalVx ;
    long double mvalY ;
    long double mPsi0 ;
    TTurnMove mTurnMove;
    TLineMove mLineMove;
    TRotating mRotating;
    THover mHover;
    long double marrRotCentre_CurNZSK[3];
    TPartHelicTraj *mpPartHelicTraj;
    double mz1 ;
    double mz2 ;
    double mval_A;
    bool mbTang;

    long double m_dPsi_po_dt;


   void inputData();





    // исходные данные по атмосфере
    double mWindHor; // горизонтьальная скорость ветра
    double mWindCourse;// курсовой угол вектора горизонтьальной скорости ветра
    double mWindVert; // вертикальная скорость ветра полож направление вверх
    double mTemperature0; // температура воздуха у поверхности земли
    TEnvironment mEnvironment;
    ///

    ///
    // папка с графиками результатов
  wchar_t mwchOutPutFold[400];

  // путь к файлу с данными аэродинамических элементов планера
  wchar_t mwchInputFileName[400];
  QString mqstrInputFileName;

  // мвссив с исходными данными аэродинамических характеристик 8 элементов планера
  double marrInpDataGliders[8*13];

 // полетная операция вертолета
  //TOperationNew mOperation;

  // файл с массивом кадидитов в передат числа
  double *mparrCandidates;
  int mnumCandidates;
  QString mqstrCandidatesSCVFIle;
  // путь к файлу с результатми отбора кагендидатов в перед числа
   wchar_t mpwchPickedOutCandidates[400];

  //
  QString mqstrPeredChislaSCVFIle;
  wchar_t mpwchPeredChislaSCVFIle[400];
private slots:  




    void on_pushButton_clicked();

    void on_doubleSpinBoxVHor_valueChanged(double arg1);

    void on_doubleSpinBoxCourseAng_valueChanged(double arg1);

    void on_doubleSpinBoxVVert_valueChanged(double arg1);

    void on_doubleSpinBoxTemperature0_valueChanged(double arg1);





    void on_doubleSpinBoxVHor_valueChanged(const QString &arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_pushButton_4_clicked();

    void on_doubleSpinBoxMass_valueChanged(const QString &arg1);

    void on_lineEdit_cursorPositionChanged(int arg1, int arg2);

    void on_pushButton_5_clicked();

    void on_pushButton_6_clicked();

    void on_pushButton_7_clicked();

    void on_pushButton_9_clicked();


    void on_pushButton_10_clicked();

    void on_pushButton_11_clicked();

    void on_pushButton_8_clicked();    

    void on_comboBox_currentIndexChanged(int index);

    void on_progressBar_valueChanged(int value);

    void on_comboBox_4_currentIndexChanged(int index);

    void on_doubleSpinBox_c_valueChanged(double arg1);

public:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
