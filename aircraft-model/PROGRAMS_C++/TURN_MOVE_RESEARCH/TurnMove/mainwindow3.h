#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "Helic.h"

#include "Environment.h"

#include "HelicTraj.h"






class QDoubleSpinBox;
class THelic;
class TEnvironment;
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




    // полетное время
    long double mTFly ;


        QDoubleSpinBox  marrDblSpinBoxBlade[6];
        ///

        // виджеты винтов
        // для ввода
        QDoubleSpinBox  marrDblSpinBoxShafts[6];




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
  THelicTraj mHelicTraj;
private slots:  


    void on_doubleSpinBoxVHor_valueChanged(double arg1);



    void on_doubleSpinBoxCourseAng_valueChanged(double arg1);

    void on_doubleSpinBoxVVert_valueChanged(double arg1);


    void on_doubleSpinBoxTemperature0_valueChanged(double arg1);

    void on_pushButton_2_clicked();








    void on_pushButton_clicked();

    void on_pushButton_3_clicked();

    void on_lineEdit_cursorPositionChanged(int arg1, int arg2);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
