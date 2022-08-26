#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "form0.h"


// https://evileg.com/ru/post/112/


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    double mx;
    double my;

private:
    Ui::MainWindow *ui;

    private slots:
        // Слоты от кнопок главного окна
        void on_pushButton_clicked();


        void on_pushButton_2_clicked();

private:

        // второе и третье окна
        Form0 *sAnotherForm2;

};

#endif // MAINWINDOW_H
