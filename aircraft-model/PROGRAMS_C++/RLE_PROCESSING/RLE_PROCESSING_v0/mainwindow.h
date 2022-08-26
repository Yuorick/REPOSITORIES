#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "Helic.h"

#include "Environment.h"


#include "TurnMove.h"
#include "PartHelicTraj.h"
#include "LineMove.h"
#include "Hover.h"
#include "Rotating.h"


class QDoubleSpinBox;
class THelic;
class TEnvironment;

class TTurnMove;
class TPartHelicTraj;
class TLineMove;
class THover;
class TRotating;

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


       // long double mvalPsi0 ;
        long double mval_dPsi0 ;
        long double mvalYaw;
        long double  mvalRad ;
        long double mvalVx ;
        long double mvalY ;
        long double m_dPsi_po_dt;
      //  TTurnMove mTurnMove;
      //  TLineMove mLineMove;
      //  TRotating mRotating;
      //  THover mHover;
        long double marrRotCentre_CurNZSK[3];
      //  TPartHelicTraj *mpPartHelicTraj;

   //THelic mHelic;
        bool mbFullGliders;
        bool mbFUllRotorModel;

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



  // файл с массивом кадидитов в передат числа
  double *mparrCandidates;
 // int mnumCandidates;
  QString mqstrCandidatesSCVFIle;
  // путь к файлу с результатми отбора кагендидатов в перед числа
   wchar_t mpwchPickedOutCandidates[400];

  //
  QString mqstrPeredChislaSCVFIle;
  wchar_t mpwchPeredChislaSCVFIle[400];
private slots:  




    void on_doubleSpinBoxVHor_valueChanged(double arg1);

    void on_doubleSpinBoxCourseAng_valueChanged(double arg1);

    void on_doubleSpinBoxVVert_valueChanged(double arg1);

    void on_doubleSpinBoxTemperature0_valueChanged(double arg1);

    void on_pushButton_2_clicked(); 

    void on_tableWidget_2_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

    void on_doubleSpinBoxVHor_valueChanged(const QString &arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);



    void on_doubleSpinBoxMass_valueChanged(const QString &arg1);

    void on_lineEdit_cursorPositionChanged(int arg1, int arg2);






    void on_pushButton_3_clicked();

    void on_pushButton_10_clicked();





    void on_pushButton_3_pressed();

    void on_comboBox_currentIndexChanged(int index);

    void on_comboBox_activated(const QString &arg1);

    void on_comboBox_editTextChanged(const QString &arg1);

    void on_pushButton_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void calcVectF_and_JacF(long double VAlAlfaZakl,long double* arrxdist0, long double valXT, long double *arrRight, long double *arrF, long double *arrJacF);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
