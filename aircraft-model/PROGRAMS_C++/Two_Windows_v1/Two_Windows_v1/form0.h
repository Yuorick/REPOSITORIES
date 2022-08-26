#ifndef FORM0_H
#define FORM0_H

#include <QWidget>

namespace Ui {
class Form0;
}

class Form0 : public QWidget
{
    Q_OBJECT

public:
    explicit Form0(QWidget *parent = nullptr);
    ~Form0();

    double *mpx;
    double *mpy;

signals:
    void firstWindow();  // Сигнал для первого окна на открытие

private slots:
    // Слот-обработчик нажатия кнопки
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::Form0 *ui;
};

#endif // FORM0_H

