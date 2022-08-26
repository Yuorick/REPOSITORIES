#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include "Receiver.h"
#include "SubWaterAnt.h"
#include "OrbitalAnt.h"
#include "MtrxAnt.h"


class QReceiver;
class QSubWaterAnt;
class QOrbitalAnt;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
void  inputData();

QString mQstrOutPutFold;
wchar_t mwchOutPutFold[400];

double mr;
int mQuant;
double ma;
double mlamb;
double mScan_q0;
double mScan_e0;

QReceiver mReceiver;

QSubWaterAnt *mpSubWaterAnt;

QOrbitalAnt mOrbitalAnt;

QMtrxAnt mMtrxAnt;
double m_e_True;
double m_e_Zv;
double m_qu_True;
double m_qu_Zv;

int mIndexcomboBox;

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton_3_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_comboBox_currentIndexChanged(int index);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
