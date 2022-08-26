#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    // Инициализируем второе окно
    sAnotherForm2 = new Form1();

    // подключаем к слоту запуска главного окна по кнопке во втором окне
    connect(sAnotherForm2, &Form1::firstWindow, this, &MainWindow::show);

}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    sAnotherForm2->show();  // Показываем второе окно
   // this->close();    // Закрываем основное окно
}
