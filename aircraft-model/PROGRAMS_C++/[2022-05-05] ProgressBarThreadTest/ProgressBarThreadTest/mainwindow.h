#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "mythread.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_pushButtonStartTask_clicked();

    void on_pushButtonStopTask_clicked();

    void onTaskProgress(int step, double error);
    void onTaskFinished();
    void onCalcFinished(int result);

private:
    Ui::MainWindow *ui;
    MyThread *myThread;
};

#endif // MAINWINDOW_H
