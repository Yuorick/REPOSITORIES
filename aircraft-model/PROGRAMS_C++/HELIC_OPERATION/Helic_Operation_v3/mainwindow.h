#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "Helic.h"
//#include "CommonGliderElement.h"
//#include "RuleGliderElement.h"
//#include "Rotor.h"
//#include "Blade.h"
#include "Environment.h"

#include "FlyTask.h"
#include "HelicTraj.h"






class QDoubleSpinBox;
class THelic;
//class TCommonGliderElement;
//class TRuleGliderElement;
//class TRotor;
//class TBlade;
class TEnvironment;

class TFlyTask;
class THelicTraj;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow();
  // виджеты полетного задания
  // QDoubleSpinBox  marrDblSpinBox[QUANT_FLY_TASK_POINTS * 5];


   // массив координат полетного задания(зарезервирован на максимальное к-во точек полетного задания)
  // double marrCoord[QUANT_FLY_TASK_POINTS * 3];
  // double marrV[QUANT_FLY_TASK_POINTS];
  // double marrCourse[QUANT_FLY_TASK_POINTS];
   TFlyTask mFlyTask;

   // широта и долгота начальной точки полетного задания
  //  long double mLong0; // долгота
  //  long double mLat0; // широта
    // полетное время
    long double mTFly ;



   THelic mHelic;

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



  // траектория
  //THelicTraj mHelicTraj;
private slots:  

   //void on_tableWidget_cellChanged(int row, int column);

   // void on_tableWidget_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

  //  void on_pushButton_clicked();

    void on_doubleSpinBoxVHor_valueChanged(double arg1);



    void on_doubleSpinBoxCourseAng_valueChanged(double arg1);

    void on_doubleSpinBoxVVert_valueChanged(double arg1);


    void on_doubleSpinBoxTemperature0_valueChanged(double arg1);

    void on_pushButton_2_clicked();


    void on_pushButton_3_clicked();

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
