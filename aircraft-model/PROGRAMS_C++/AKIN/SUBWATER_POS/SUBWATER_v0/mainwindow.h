#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

class TTable_1D;
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

     wchar_t mpwchPrflFIle[400];

     wchar_t mwchOutPutFold[400];

private slots:
    void on_pushButton_clicked();

    void on_pushButton_5_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;
};
double f(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

double f1(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

double f2(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

double F(const double z, const double VAlza
         ,TTable_1D &tblPrfl,const double VAlCosTetta);

#endif // MAINWINDOW_H
