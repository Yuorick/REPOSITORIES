#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);


    // Инициализируем второе окно
    sAnotherForm2 = new Form0();

    // подключаем к слоту запуска главного окна по кнопке во втором окне
    connect(sAnotherForm2, &Form0::firstWindow, this, &MainWindow::show);

    mx = 5.;
    my = 0.;
    sAnotherForm2->mpx = &mx;
    sAnotherForm2->mpy = &my;
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_clicked()
{

    sAnotherForm2->show();  // Показываем второе окно
    this->close();    // Закрываем основное окно
}

void MainWindow::on_pushButton_2_clicked()
{
    double z =  my;
    int ii=0;
}
