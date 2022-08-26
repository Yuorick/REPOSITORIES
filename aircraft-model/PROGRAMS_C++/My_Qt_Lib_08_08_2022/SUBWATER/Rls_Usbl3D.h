#ifndef RLS_USBL3D_H
#define RLS_USBL3D_H
#include "HidroRLS.h"


class QRls_Usbl3D : public QHidroRLS
{
public:
    QRls_Usbl3D();

    QRls_Usbl3D (const  QRls_Usbl3D &R);

    QRls_Usbl3D  &operator=( const QRls_Usbl3D  &R);

    QRls_Usbl3D  (const double SigR, const double DiagWidth
                             ,const double Sig_q, const double Sig_e);

    virtual void calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                  , double &Sig_q, double &Sig_e);

    // ТОчность по углу курса (betta)
        double mSig_q;
    //  точность определения угла места (eps)
        double mSig_e;

};

#endif // RLS_USBL3D_H
