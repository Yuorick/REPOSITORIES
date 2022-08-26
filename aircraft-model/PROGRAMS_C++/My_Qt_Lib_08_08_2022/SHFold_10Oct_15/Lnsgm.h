#ifndef LNSGM_H
#define LNSGM_H
// класс описывает отрезок , лежащий на прямой

class QLnSgm
{
public:
    double ma;
    double mb;


    QLnSgm();

    QLnSgm (const QLnSgm &R2);

    // оператор присваивания
    QLnSgm &operator=(const QLnSgm  &R2);

    // парам констр
    QLnSgm(const double a1,const double b1);

    double calcLeng();

    static bool productOfSgms( QLnSgm sgm0,  QLnSgm sgm1, QLnSgm* sgmRez);

    bool  isValInsideSgm(const double x);

    static int SegmMinusSegm( QLnSgm sgm0, QLnSgm sgm1, QLnSgm* parrSgmRez);
};

#endif // LNSGM_H
