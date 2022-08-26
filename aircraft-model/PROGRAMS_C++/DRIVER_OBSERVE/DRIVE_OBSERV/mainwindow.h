#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "ElectDriver.h"
#include "DriveTraj.h"




class QDoubleSpinBox;
class QlectDriver;
class QDriveTraj;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow();

    // параметры  мат модели привода
       // индуктивность катушки статора
        double mInductL;
        // сопротивление катушки статора
         double mResist;

         // индуктивный поток
         double mPsi_f;

         //желаемая угловая скорость вращения
         double mOmegaStat;
         //
         double mJPayLoad;
         //
         double mCx_om;
         //
         bool mbVelo;
         QElectDriver mDriverModel;

         // параметры  реального привода
            // индуктивность катушки статора
             double mRealInductL;
             // сопротивление катушки статора
              double mRealResist;

              // индуктивный поток
              double mRealPsi_f;

              //желаемая угловая скорость вращения
              double mRealOmegaStat;
              //
              double mRealJPayLoad;
              //
              double mRealCx_om;
              //

              QElectDriver mDriverReal;

              //
              QDriveTraj mDriveTrajReal;
    // полетное время
        double mTRotate ;




   void inputData();





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






    void on_pushButton_2_clicked(); 




    void on_doubleSpinBox_2_valueChanged(double arg1);



    void on_doubleSpinBoxMass_valueChanged(const QString &arg1);

    void on_lineEdit_cursorPositionChanged(int arg1, int arg2);






    void on_pushButton_3_clicked();

    void on_pushButton_10_clicked();





    void on_pushButton_3_pressed();

    void on_comboBox_currentIndexChanged(int index);

    void on_comboBox_activated(const QString &arg1);

    void on_comboBox_editTextChanged(const QString &arg1);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
