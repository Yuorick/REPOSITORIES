#include "form0.h"
#include "ui_form0.h"

Form0::Form0(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::Form0)
{
    ui->setupUi(this);
}

Form0::~Form0()
{
    delete ui;
}


void Form0::on_pushButton_clicked()
{
   // this->close();      // Закрываем окно
   // emit firstWindow(); // И вызываем сигнал на открытие главного окна
}

void Form0::on_pushButton_2_clicked()
{
    //double temp = (*mpx)*10.;
    //(*mpy) = temp;
}
