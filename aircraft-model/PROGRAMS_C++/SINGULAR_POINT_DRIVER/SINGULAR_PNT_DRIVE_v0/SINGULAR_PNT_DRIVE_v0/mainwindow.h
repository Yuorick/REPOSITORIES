#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDoubleSpinBox>
#include "ElectDriver.h"
#include "Ctrl.h"
#include "Comp.h"
#include <QTableWidget>
//#include "CtrlVeloQuiq.h"
//#include "CtrlVelo.h"
#include "Ctrl.h"
//#include "CtrlPos.h"
#include "CtrlFollow.h"
#include "MeasurmentImitator.h"
#include "Filtr.h"
//#include "DriveMoveImit.h"

//#include "CtrlTrad.h"
#include "Environment.h"
#include "AbstractTraj.h"
#include "PrimitiveCirleUniformTraj.h"
#include "EquallyAcceleratedMotion.h"
#include "WarShip.h"
#include "WarSins.h"
#include "CtrlFollowAdapt1.h"
#include "UniformRectLinMotion.h"
#include "MobilePlatf.h"
#include "OperativeCtrl.h"
#include "MeasInstr.h"
#include "OEChannel.h"
#include "MobilePlatfCtrl.h"


class QDoubleSpinBox;
class QlectDriver;
class QCtrl;
class QBrls;


class QCtrl;

class TComp;
class QCtrlFollow;
class QMeasurmentImitator;
class QFiltr;
//class QDriveMoveImit;
class QAbstractTraj ;
class QPrimitiveCirleUniformTraj;
class QEquallyAcceleratedMotion;
class QWarShip;
class QWarSins;
class QCtrlFollowAdapt1;
class QUniformRectLinMotion;
class QDevice;
class QTargTrackingOperation;


class TEnvironment;

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
    double mModelInductL;
    // сопротивление катушки статора
    double mModelResist;

    // индуктивный поток
    double mModelPsi_f;
    // момент сопротивления эл двигателя при обесточенных обмотках, не более
    double mModelMaxAmpMomResidual ;
    // фаза моменат сопротивления эл двигателя при обесточенных обмотках, не более
    double mModelPhMomResidual;
    //мом инерции ротора
    double mModelJ0;
    // момент сухого трения
    double mMomDryFriction;

    // масса штанги
    double mMassBar;

    // НАГРУЗКА
    //мом инерции нагрузки
    double mModelJPayLoad;
    // коэффициент сопротивления нагрузки
    double mModelCx;
    //масса левой нвавески
    double mMassLeft;
    //масса правой навескт
    double mMassRight;
    //общая длина штанги
    double mBarLength;
    //высота штанги
    double mBarHight;
    // коэффициент трния нагрузки
    double mModelCv;
    // внешний момент, заложенный в модель
   // double mModelMomOut;


   // QElectDriver mDriverModel;
   // QElectMotor mElectMotorModel;
   // QLoad mLoadModel;

    QElectDriver mModelOecHorDriver;
    QElectMotor mModelOecHorElectMotor;

    QElectDriver mModelOecVertDriver;
    QElectMotor mModelOecVertElectMotor;

    QLoad mHorLoadModel;
    QLoad mVertLoadModel;

    QLoad *mpHorLoadModel;
    QLoad *mpVertLoadModel;

    QBrls mHorBrlsModel;
    QBrls mVertBrlsModel;


/////////////////////////////////////////////////////////////////


    QLoad mHorLoadReal;
    QLoad mVertLoadReal;

    QLoad *mpHorLoadReal;
    QLoad *mpVertLoadReal;

    QBrls mHorBrlsReal;
    QBrls mVertBrlsReal;



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
   // double mRealMomOut;

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

    QElectDriver mOecHorDriver;
    QElectMotor mOecHorElectMotor;

    QElectDriver mOecVertDriver;
    QElectMotor mOecVertElectMotor;


    // начальное время
    double mT0;
    // полетное время
    double mMovingT ;

    TComp mCmpArrLamb[4];

    //начальнвая  угловая скорость вращения
    double mOmegaBegin;

    //  начальное угловое пложение
    double mTettaBegin;

    QBrls mBrls;
   // QBrls mBrlsReal;



    // папка с графиками результатов
    wchar_t mwchOutPutFold[400];

    // путь к файлу с данными аэродинамических элементов планера
    wchar_t mwchInputFileName[400];
    QString mqstrInputFileName;

    bool mbtableWidget_2Init;



// управление
    QOperativeCtrl mOperCtrl;

    QMobilePlatfCtrl mOecMobilePlatfCtrl;

    QCtrl* mpHorCtrl;
    QCtrlFollow mCtrlHorFollow;
    QCtrlFollowAdapt1 mCtrlHorFollowAdapt1;

    QCtrl* mpVertCtrl;
    QCtrlFollow mCtrlVertFollow;
    QCtrlFollowAdapt1 mCtrlVertFollowAdapt1;
///
    // вектор цели слежения привода для движения по окружности
   // double marrObjective[2];
  //  double mTargR;

    // имитатор измерительнрой сиситемы

    QMeasurmentImitator mDriveMeasImitator;

    // фильтр
    QFiltr mFiltr;
    // шаг управления/фильтрации
     double mh;


    // вектор разбросов параметров двигателя
    double  marrSpreadParams[QPARAMS];

    // тип идентификации параметров
    int mTypeOfID;

    TEnvironment mEnvironment;

    QAbstractTraj *mpTargTraj;
    QPrimitiveCirleUniformTraj mPrimitiveCirleUniformTraj;
    QUniformRectLinMotion mUniformRectLinMotion;
    QEquallyAcceleratedMotion mEquallyAcceleratedMotion;

    // вектор параметров позиционироваия привода истинный
    double marrOecPosParams[6];
    // вектор параметров позиционироваия привода оценочный
    double marrOecPosParams_Zv[6];

    //корабль
    QWarShip mWarShip;
    double mvalMaxQ;
    double mvalMaxTet;
    double mvalMaxPsi;

    // СИНС
    QWarSins mSins;
    // темп выдачи информации СИНС
    double mTimeTempSins;

    QPlatf *mpPltfOec;
    QMobilePlatf mMobilePlatfOec;
   // QMeasInstr* mpMeasInstr;
    QOEChannel mOEChannel;
    QDevice *mpDevice;

    //QElectDriver  mOecHorDriver;
    //QElectDriver  mOecVertDriver;

    // алгоритм фильтрации/ если mbPSK == true, то фильтрация происходит в ПСК, если false, то в КГСК
   // bool mbPSK;

    bool mbtableWidget_Init;

    bool mbtableWidget_4Init;


    void inputData();

private slots:

    void on_pushButton_2_clicked();



    void on_pushButton_3_clicked();







    void on_comboBox_3_currentIndexChanged(int index);





    void on_tableWidget_2_cellChanged(int row, int column);

    void read_tableWidget_2(double *arr); 

   void init_tableWidget_2(const int ncols, double *arr);





   void on_comboBox_4_currentIndexChanged(int index);

   void on_tableWidget_4_cellChanged(int row, int column);

   void read_tableWidget_4(double *arr);

   void write_tableWidget_4(double *arr);

   void on_tableWidget_cellChanged(int row, int column);

   int createInputDataReport(wchar_t*FileName, const bool bHeader, QTargTrackingOperation &TargTrackingOperation);

   int createOutputDataReport(wchar_t*FileName, const bool bHeader
                              ,const double* arrHorSyst,const double* arrHorDisp
                              ,const double*  arrVertSyst,const double* arrVertDisp);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
