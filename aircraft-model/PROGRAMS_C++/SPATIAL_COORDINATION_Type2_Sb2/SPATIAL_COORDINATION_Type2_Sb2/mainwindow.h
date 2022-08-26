#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "ContObj6.h"

#define MAX_QUANT_CNTRL_OBJ 10000


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    //

    // вектор истинных параметров пространственной ориентации грани РЛК
    // в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXTrue_IU_PSK[QUANT_GDG_PRMS];


    // вектор первичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXZv_IU_PSK[QUANT_GDG_PRMS];






    // вектор вторичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК
    double marrXEst_IU_PSK[QUANT_GDG_PRMS];

    // корреляционная матрица ошибок оценивания
    double marrMtrxK[QUANT_GDG_PRMS * QUANT_GDG_PRMS];

    // Кол-во контрольных объектов (КО)
    int mNumCtrlObj;

    // высота КО
    double mCtrlObjHeight;

    // массив исходных данных КО из входной таблицы
    double mparrDataCtrlObj[MAX_QUANT_CNTRL_OBJ * 5];

    // Коррел матрица ошибок первичных измерений РЛК
    double marrCorMtrxRLK[9];


    // СКЗ ошибки по углу ОЭК
    double mAngleSigIU;

    // СКЗ ошибки дальномера
    double mDistSigIU;

    //признак случайных шумов измерений
    bool mbNoise;

    // вектор дисперсий ИС - Q, Psi, Tetta
    double marrDispEilerCntrKP[3];

    // дисп матрица ошибок КО в КГСК
  //  double  marrCorMtrx_CO_KGSK[9];
  //
    bool mbFirstTime;

// по калькулятору
    // признак того, что вычисления закончены и результаты сохранены в таблице
  //  bool mbCalculationIsDone;

    //  СКЗ ошибки позиционирования корабля-носителя

    double mVessSKZ;
    //  СКЗ ошибки позиционирования КО

    double mCntrlObjSKZ;

    bool mbStar;

    double marrOrtStar[3];




    void setupFirst();

    void fillWindow();

    void inputData();



private slots:
    void on_doubleSpinBox_2_valueChanged(const QString &arg1);

    void on_pushButton_clicked();

    void on_tableWidget_2_cellChanged(int row, int column);

    void on_tableWidget_cellChanged(int row, int column);

    void on_checkBox_4_stateChanged(int arg1);

    void on_doubleSpinBox_2_valueChanged(double arg1);

    void on_doubleSpinBox_2_editingFinished();

    void on_doubleSpinBox_3_valueChanged(double arg1);

    void on_doubleSpinBox_4_valueChanged(const QString &arg1);

    void on_checkBox_stateChanged(int arg1);

    void on_checkBox_5_stateChanged(int arg1);

    void on_checkBox_2_stateChanged(int arg1);

    //void on_pushButton_2_clicked();

    bool checkObj();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
