#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "ElectDriver.h"
#include "Ctrl.h"
#include "Comp.h"
#include <QTableWidget>
#include "CtrlVeloQuiq.h"
#include "CtrlVelo.h"
#include "Ctrl.h"
#include "CtrlPos.h"
#include "CtrlFollow.h"
#include "MeasurmentImitator.h"
#include "Filtr.h"
#include "DriveMoveImit.h"
#include "CtrlDnmkCorrect_M_R.h"
#include "CtrlDnmkCorrect_W.h"
#include "CtrlTrad.h"
#include "Environment.h"
#include "AbstractTraj.h"
#include "PrimitiveCirleUniformTraj.h"


class QDoubleSpinBox;
class QlectDriver;
class QCtrl;
class QBrls;
class QCtrlVeloQuiq;
class QCtrlVelo;
class QCtrl;
class QCtrlPos;
class TComp;
class QCtrlFollow;
class QMeasurmentImitator;
class QFiltr;
class QDriveMoveImit;
class QCtrlDnmkCorrect_M_R;
class QCtrlDnmkCorrect_W;
class QCtrlTrad;
class TEnvironment;
class QAbstractTraj;
class QPrimitiveCirleUniformTraj;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);

    ~MainWindow();



    // индуктивность катушки статора
    double mInductL;
    // сопротивление катушки статора
    double mResist;

    // индуктивный поток
    double mPsi_f;
    // момент сопротивления эл двигателя при обесточенных обмотках, не более
    double mMaxAmpMomResidual ;
    // фаза моменат сопротивления эл двигателя при обесточенных обмотках, не более
    double mPhMomResidual;
    //мом инерции ротора
    double mJ0;
    // момент сухого трения
    double mMomDryFriction;

    // масса штанги
    double mMassBar;

    // НАГРУЗКА
    //мом инерции нагрузки
    double mJPayLoad;
    // коэффициент сопротивления нагрузки
    double mCx;
    //масса левой нвавески
    double mMassLeft;
    //масса правой навескт
    double mMassRight;
    //общая длина штанги
    double mBarLength;
    //высота штанги
    double mBarHight;
    // коэффициент трния нагрузки
    double mCv;
    // внешний момент, заложенный в модель
    double mModelMomOut;


    // тип управления
    bool mbIsVeloCtrl;


    // максимум напряжения Uqu во время пеходношго процесса
    double mQuiqMaxU0;


    QElectDriver mDriverModel;
    QElectMotor mElectMotorModel;
    QLoad mLoadModel;
    QLoad mLoadReal;


    QLoad *mpLoadModel;
    QLoad *mpLoadReal;

    // параметры  реального привода
    // индуктивность катушки статора
    double mRealInductL;
    // сопротивление катушки статора
    double mRealResist;

    // индуктивный поток
    double mRealPsi_f;
    // момент сопротивления эл двигателя при обесточенных обмотках, не более
    double mRealMaxAmpMomResidual;
    // фаза моменат сопротивления эл двигателя при обесточенных обмотках, не более
    double mRealPhMomResidual;
    //мом инерции ротора
    double mRealJ0;
    // масса штанги
    double mRealMassBar;
    // внешний момент, заложенный в модель
    double mRealMomOut;

    // НАГРУЗКА
    //мом инерции нагрузки
    double mRealJPayLoad;
    // коэффициент сопротивления нагрузки
    double mRealCx;
    //масса левой нвавески
    double mRealMassLeft;
    //масса правой навескт
    double mRealMassRight;
    //общая длина штанги
    double mRealBarLength;
    //высота штанги
    double mRealBarHight;
    // коэффициент тренния нагрузки
    double mRealCv;
    //

    QElectDriver mDriverReal;

    QElectMotor mElectMotorReal;

    //
    QCtrl mCtrlReal;
    //  QCtrl mCtrlModel;
    // полетное время
    double mMovingT ;

    TComp mCmpArrLamb[4];

    //начальнвая  угловая скорость вращения
    double mOmegaBegin;

    //  начальное угловое пложение
    double mTettaBegin;

    QBrls mBrls;
    QBrls mBrlsReal;

    // папка с графиками результатов
    wchar_t mwchOutPutFold[400];

    // путь к файлу с данными аэродинамических элементов планера
    wchar_t mwchInputFileName[400];
    QString mqstrInputFileName;

    bool mbtableWidget_2Init;
    bool mCtrlVelo;

    QCtrlVeloQuiq mModelCtrlVeloQuiq;
    QCtrlVelo mModelCtrlVelo;
    QCtrl mModelCtrl;
    QCtrlPos mModelCtrlPos;
    QCtrlFollow mModelCtrlFollow;    
    QCtrlDnmkCorrect_M_R mModelCtrlDnmkCorrect_M_R;
    QCtrlDnmkCorrect_W mModelCtrlDnmkCorrect_W;
    QCtrlTrad mModelCtrlTrad;
    QCtrl* mpModelCtrl;

    QCtrlVeloQuiq mRealCtrlVeloQuiq;
    QCtrlVelo mRealCtrlVelo;
    QCtrl mRealCtrl;
    QCtrlPos mRealCtrlPos;
    QCtrlFollow mRealCtrlFollow;
    QCtrlTrad mRealCtrlTrad;

    QCtrl* mpRealCtrl;

    // интервал обновления вычислений параметров по Mout и R
    double mTPeriodParamsRenew;
    // постоянная времени интеграторов
    double mTimeIntegrator;

    // вектор цели слежения привода в начале движения
    double marrObjective[2];


    // имитатор измерительнрой сиситемы
    QMeasurmentImitator mMeasurmentImitator;

    // фильтр
    QFiltr mFiltr;
    // шаг управления/фильтрации
     double mh;

    QDriveMoveImit mDriveMoveImit;


    // вектор разбросов параметров двигателя
    double  marrSpreadParams[QPARAMS];

    // тип идентификации параметров
    int mTypeOfID;

    TEnvironment mEnvironment;

    QAbstractTraj *mpTargTraj;
    QPrimitiveCirleUniformTraj mPrimitiveCirleUniformTraj;



    void inputData();

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


    void on_comboBox_3_currentIndexChanged(int index);

    void on_comboBox_3_currentIndexChanged(const QString &arg1);

    void on_tableWidget_2_currentCellChanged(int currentRow, int currentColumn, int previousRow, int previousColumn);

    void on_tableWidget_2_cellChanged(int row, int column);

    void read_tableWidget_2(double *arr); 

   void init_tableWidget_2(const int ncols, double *arr);

   void on_comboBox_currentIndexChanged(const QString &arg1);



private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
