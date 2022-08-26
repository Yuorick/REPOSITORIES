#ifndef RLS_USBL2D_H
#define RLS_USBL2D_H
#include "HidroRLS.h"


class QRls_Usbl2D : public QHidroRLS
{
public:
    QRls_Usbl2D();

    QRls_Usbl2D (const  QRls_Usbl2D &R);

    QRls_Usbl2D  &operator=( const QRls_Usbl2D  &R);

    QRls_Usbl2D  (const double SigR, const double DiagWidth,const double SigV);

    virtual void calc_qZv_and_eZv(const QTrueMeasParams trueMeasParams, double &qZv, double &eZv
                                  , double &Sig_q, double &Sig_e);

    // ТОчность по углу курса (betta)
        double mSig_q;




};

#endif // RLS_USBL2D_H
