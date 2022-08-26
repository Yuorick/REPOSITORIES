#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#define MAX_QUANT_CNTRL_OBJ 1000

#include "Adjustment.h"


#define LENX_V1 9


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
    double marrXTrue_RLK_PSK[NUM_GADG_PARAMS];

    // вектор истинных параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXTrue_IU_PSK[NUM_GADG_PARAMS];

    // вектор первичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXZv_RLK_PSK[NUM_GADG_PARAMS];

    // вектор первичных оценок параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXZv_IU_PSK[NUM_GADG_PARAMS];



    // массив ланных пеленга цели - угол пеленга, дальность, высота
    double marrBaring[3];

    // палубные углы - угол курса, килевой качки, бортовой качки
    double marrDeckAngles[3];



    // массив углов наведения истинный - Betta, Eps
    double marrCntrlAngTrue_V1[2];

    // массив углов наведения оценка первоначальная - Betta, Eps
    double marrCntrlAngZv_V1[2];



    // матрица частных производных углов наведения по параметрам позиционирования V1
     double marr_dVIU_po_dX_V1[2 * LENX_V1];

/////////////////////////////////////////////////////////////////////////////////////
     // массив углов наведения истинный - Betta, Eps
     double marrCntrlAngTrue_V2[2];

     // массив углов наведения оценка первоначальная - Betta, Eps
     double marrCntrlAngZv_V2[2];


     // матрица частных производных углов наведения по параметрам позиционирования V2
      double marr_dVIU_po_dX_V2[2 * NUM_GADG_PARAMS];



    void setupFirst();

    void fillWindow();

    void inputData();

   // void fill_Partly_OutTable();

private slots:

    void fillCalculatorArraysPrevious();

    void fillCalculatorFrame();

    void on_pushButton_2_clicked();

    void calculatorInput();


private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
