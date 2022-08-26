#include "mythread.h"
#include <QDebug>

MyThread::MyThread(QObject *parent) :
    QThread(parent),
  userBreakRequested(false),
  maxStep(1),
  errorThreshold(0.1)
{
}

void MyThread::setMaxStep(int value)
{
    maxStep = value;
}

int MyThread::getMaxStep() const
{
    return maxStep;
}

void MyThread::setErrorThreshold(double value)
{
    errorThreshold = value;
}

double MyThread::getErrorThreshold() const
{
    return errorThreshold;
}

void MyThread::run()
{
    qDebug() << "MyThread::run(), in";

    int result = MYTHREAD_RESULT_MAX_STEP; // результат по умолчанию
    for(int i = 0; i <= maxStep; i++ )
    {
        qDebug() << "MyThread::run()" << i;

        // !следующую строку следует заменить на реальные вычисления!
        QThread::msleep(100); // этот вызов здесь нужен для имитации длительности вычислений

        // после очередной итерации следует вычислить ошибку (невязку)
        double error = 1/double(1+i);

        // отправляем в интерфейс сигнал о ходе выполнения задачи
        emit progress(i, error);

        // проверяем различные условия досрочного выхода из цикла
        if(error < errorThreshold)
        {
            result = MYTHREAD_RESULT_OK;
            break; // если достигли нужной точности, завершаем вычисления
        }
        if(userBreakRequested)
        {
            result = MYTHREAD_RESULT_USER_CANCEL;
            break; // если пользователь потребовал, завершаем вычисления
        }
    }
    // отправляем в интерфейс результат выполнения задачи
    emit calcfinished(result);

    qDebug() << "MyThread::run(), out";
}
