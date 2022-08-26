#ifndef CAVITYANTENNA_H
#define CAVITYANTENNA_H


class QCavityAntenna
{
public:
    QCavityAntenna();

    QCavityAntenna (const QCavityAntenna &R);

    QCavityAntenna &operator=(const QCavityAntenna  &R);

    QCavityAntenna (const double Lin_a, const double Lin_d
                          ,  const double DM_a,  const double DM_L,  const double DM_d);

    // сечение канала линейки (большая стенка), [м]
    double mLin_a;
    // шаг между щелями, [м]
    double mLin_d;
    // сечение канала делителя мощности, [м]
    double mDM_a;
    // длина витка делителя мощности, [м]
    double mDM_L;
    // шаг выходов в делители мощности, [м]
    double mDM_d;

    void createInpDataArray (double *arrInpData) ;

    double calcLamb_w(const double VAlLamb);

    double calcU(const double VAlLamb);

    double calcV(const double VAlLamb);

    double calcLambDM(const double VAlLamb);

    void calc_DM_params_array(const double VAlLamb, double *arrParamsDM);

    int calcDM_k(const double VAlLamb);

};

#endif // CAVITYANTENNA_H
