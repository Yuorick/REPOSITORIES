#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#define MAX_QUANT_CNTRL_OBJ 10000
#include "ContObj6.h"

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
    double marrXTrue_RLK_PSK[QUANT_GDG_PRMS];

    // вектор истинных параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXTrue_IU_PSK[QUANT_GDG_PRMS];

    // вектор первичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXZv_RLK_PSK[QUANT_GDG_PRMS];

    // вектор первичных оценок параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXZv_IU_PSK[QUANT_GDG_PRMS];


    // вектор истинных параметров   пространственной ориентации
    // грани РЛК относительно основания  ИУ
    // в ПСК-ИУ -вектор параллакса РЛК, угол BettaРЛК-BettaИУ, угол EpsРЛК-EpsИУ

   double marrXRLK_RelativeTrue[QUANT_GDG_PRMS];


    // вектор первичных оценок параметров пространственной ориентации
    // грани РЛК относительно основания  ИУ
    // в ПСК-ИУ -вектор параллакса РЛК, угол BettaРЛК-BettaИУ, угол EpsРЛК-EpsИУ
    double marrXRLK_RelativeZv[QUANT_GDG_PRMS];

    // вектор вторичных оценок параметров пространственной ориентации грани РЛК отнсительно основания ИУ
    // в ПСК-ИУ
    double marrXRLK_RelativeEst[QUANT_GDG_PRMS];



    double marrXEst[QUANT_GDG_PRMS];

    // вектор вторичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК
    double marrXEst_RLK[QUANT_GDG_PRMS];



    // корреляционная матрица ошибок оценивания
    double marrMtrxK[QUANT_GDG_PRMS *QUANT_GDG_PRMS];

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


    bool mbFirstTime;


     double mVessSKZ;
     double mCntrlObjSKZ;

     // по калькулятору
         // признак того, что вычисления закончены и результаты сохранены в таблице
         bool mbCalculationIsDone;

         // массив ланных пеленга цели - угол пеленга, дальность, высота
         double marrBaring[3];

         // палубные углы - угол курса, килевой качки, бортовой качки
         double marrDeckAngles[3];

         // массив углов наведения истинный - Betta, Eps
         double marrCntrlAngTrue[2];

         // массив углов наведения оценка первоначальная - Betta, Eps
        // double marrCntrlAngZv[2];

         // массив углов наведения оценка уточненная - Betta, Eps
         double marrCntrlAngEst[2];

         //массив углов наведения оценка предварительная - Betta, Eps
         double marrCntrlAngPrevEst[2];

         // матрица частных производных углов наведения по параметрам позиционирования
          double marr_dVIU_po_dX[2 * QUANT_GDG_PRMS];



    void setupFirst();

    void fillWindow();

    void inputData();

  //  void fill_Partly_OutTable();

private slots:


    void on_pushButton_clicked();

    void on_tableWidget_2_cellChanged(int row, int column);

    void on_tableWidget_cellChanged(int row, int column);

    void on_checkBox_4_stateChanged(int arg1);



    void on_doubleSpinBox_2_editingFinished();

    void on_doubleSpinBox_3_valueChanged(double arg1);

    void on_doubleSpinBox_4_valueChanged(const QString &arg1);

    void fillCalculatorFrame();

    void fillCalculatorArraysPrevious();

    void calculatorInput();

    void on_pushButton_3_clicked();

    void on_doubleSpinBox_2_valueChanged(double arg1);

    //void on_pushButton_2_clicked();

    bool checkObj();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
