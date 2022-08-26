#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    myThread(0)
{
    ui->setupUi(this);

    // очистка интерфейса
    ui->pushButtonStartTask->setEnabled(true);
    ui->pushButtonStopTask->setEnabled(false);
    ui->progressBar->setValue(0);
    ui->labelTaskInfo1->setText("");
    ui->labelTaskInfo2->setText("");
}

MainWindow::~MainWindow()
{
    delete ui;
}

// обработчик нажатия кнопки Start (запуск задачи)
void MainWindow::on_pushButtonStartTask_clicked()
{
    qDebug() << "MainWindow::on_pushButtonStartTask_clicked";

    int maxStep = 1000;

    // очистка интерфейса
    ui->pushButtonStartTask->setEnabled(false);
    ui->pushButtonStopTask->setEnabled(true);
    ui->progressBar->setValue(0);
    ui->progressBar->setRange(0, maxStep);
    ui->labelTaskInfo1->setText("In progress...");
    ui->labelTaskInfo2->setText("");

    // создаем поток и привязываем сигналы к слотам
    myThread = new MyThread;
    QObject::connect(myThread, SIGNAL(finished()), this, SLOT(onTaskFinished()));
    QObject::connect(myThread, SIGNAL(progress(int,double)), this, SLOT(onTaskProgress(int,double)));
    QObject::connect(myThread, SIGNAL(calcfinished(int)), this, SLOT(onCalcFinished(int)));

    // задаём параметры выполнения задачи
    // здесь можно передать массив измерений, профиль скорости звука и пр. параметры
    myThread->setMaxStep(maxStep);
    myThread->setErrorThreshold(0.01);

    // запускаем поток (при этом автоматически запустится функция MyThread::run())
    myThread->start();
}

// обработчик нажатия кнопки Stop (принудительная остановка задачи пользователем)
void MainWindow::on_pushButtonStopTask_clicked()
{
    qDebug() << "MainWindow::on_pushButtonStopTask_clicked";
    myThread->userBreak(); // устанавливаем в потоке флаг завершения
}

// обработчик сигнала о текущем прогрессе выполнения задачи
void MainWindow::onTaskProgress(int step, double error)
{
    qDebug() << "MainWindow::onTaskProgress" << step << error;
    // вывод в интерфейс
    ui->progressBar->setValue(step);
    ui->labelTaskInfo2->setText(QString("Step %1, Error %2").arg(step).arg(error, 0, 'f', 6));
}

// обработчик сигнала о факте завершении потока (этот сигнал поток автоматически отправляет при выходе из функции MyThread::run())
void MainWindow::onTaskFinished()
{
    qDebug() << "MainWindow::onTaskFinished";
    // интерфейс
    ui->pushButtonStartTask->setEnabled(true);
    ui->pushButtonStopTask->setEnabled(false);

    // удаляем завершившийся поток, обнуляем ссылку
    delete myThread;
    myThread = 0;
}

// обработчик сигнала о факте завершения и результате вычислений
void MainWindow::onCalcFinished(int result)
{
    switch(result)
    {
    case MYTHREAD_RESULT_OK:
        ui->labelTaskInfo1->setText("Task finished successfully");
        break;
    case MYTHREAD_RESULT_MAX_STEP:
        ui->labelTaskInfo1->setText("Task failed (maximum iterations reached)");
        break;
    case MYTHREAD_RESULT_USER_CANCEL:
        ui->labelTaskInfo1->setText("Task failed (cancelled by user)");
        break;
    default:
        ui->labelTaskInfo1->setText("Task failed (unknown reason)");
        break;
    }
}
