#ifndef MYTHREAD_H
#define MYTHREAD_H

#include <QThread>

enum MYTHREAD_RESULT
{
    MYTHREAD_RESULT_OK,         // вычисления завершены успешно (достигли требуемой точности)
    MYTHREAD_RESULT_MAX_STEP,   // вычисления провалены (достигли максимально допустимого кол-ва итераций)
    MYTHREAD_RESULT_USER_CANCEL // вычисления провалены (прервано пользователем)
};

class MyThread : public QThread
{
    Q_OBJECT
public:
    explicit MyThread(QObject *parent = 0);

    void setMaxStep(int value);
    int getMaxStep() const;
    void setErrorThreshold(double value);
    double getErrorThreshold() const;

    void userBreak() { userBreakRequested = true; } // установка флага завершения работы потока

signals:
    void progress( int step, double error); // step - номер итерации, error - ошибка (невязка)
    void calcfinished(int result);

public slots:

private:
    void run();

private:
    bool userBreakRequested;
    int maxStep; // максимально допустимое количество итераций
    double errorThreshold; // требуемая точность

};

#endif // MYTHREAD_H
