#ifndef QORIENTATION_H
#define QORIENTATION_H

#define VELO_C  299792458.0
class QCavityAntenna;
class QOrientation
{
public:
    QOrientation();

    static void calcUnitVector(const double VAlV, const double VAlU, double *arrUnit);

    //static void calcEilersAngs( double *arr_e0, double *arr_e1, double *arrMu, double *matrRez);

    static void processingExperimentalData(QCavityAntenna cavityAntenna, double *arrBaseFreqData
        ,const int quantBaseFreq,double *arrWorkFreq, const int quantWorkFreq
        , const int IPredeterminationType, const double* ARrPredeterminatedAngs
        ,double *arrDataOut, double *arrOutMu, QCavityAntenna *pCavityAntennaEst, double *arrNeviaz);

    static double calcLamb_(const double f);

    static double estimate_aLin_and_dLin(const double *ArrBaseFreqData,const int quantBaseFreq
                                         , double *pval_aLin, double *pval_dLin);

    static double calcFGr_u(const double *ArrBaseFreqData,const int quantBaseFreq, const double *arrMu);

    static void calcLin_dF_and_JacF(const double *ArrBaseData,const int quantBaseFreq
                                      ,const double * arrMu,double * arrF,double * arrJacF);



    static double fncFiMu(const double lamb, const double uZv, const double *arrMu) ;




    static void calc_a_and_d_Lin_Prev(const double lamb1, const double lamb2
                                                    ,const double u1, const double u2, double* arrMu );

    static double fnc_d(const double lamb1,const double u1,const double a);

    static void calcLin_dFi_and_JacFi(const double lamb, const double u
                                                    , const double* arrMu, double*arrFi, double*arrJacFi );

    static double estimate_aDM_and_dDM_and_LDM(const double* arrBaseData,const int quantBaseFreq
               ,const int Ik, double*val_aDM,double*val_dDM,double*val_LDM, double *arrDM_params);

    static bool calc_DM_params_init(const double lamb1, const double lamb2, const double lamb3
                                              ,const double u1, const double u2, const double u3,double*  arrx );

    static double calc_f(const double ly1,const double ly2,const double  lamb1
                                       , const double lamb2, const double lamb3,const double x);

    static double calcG(const double  lam, const double mu, const double x);

    static double calc_g(const double  lam, const double x);

    static double calc_dg(const double  lam, const double x);

    static double calc_d2g(const double  lam, const double x);

    static double calc_df(const double y1,const double y2,const double  lamb1
                                       , const double lamb2, const double lamb3,const double x);

    static double calcdG(const double  lam, const double mu, const double x);

    static double calc_d2f(const double y1,const double y2,const double  lamb1
                                         , const double lamb2, const double lamb3,const double x);

    static double calcd2G(const double  lam, const double mu, const double x);

    static double calcFGr_V(const double *ArrBaseFreqData,const int quantBaseFreq, const double *arrDM_params);

    static double fncFi_DM(const double lamb, const double VZv, const double *arrDM_params);

    static void calcDM_dF_and_JacF(const double *ArrBaseData,const int quantBaseFreq
                                                 ,const double * arrDM_params,double * arrF,double * arrJacF);

    static void calcDM_dFi_and_JacFi(const double lamb, const double u
          , const double* arrDM_params, double*arrdFi, double*arrJacFi );

    static void calcRotationAngs( double *arr_e0, double *arr_e1
                                                , const int IPredeterminationType, const double* ARrPredeterminatedAngs
                                                , double *arrMu, double *matrRez);

    static double fncV_params(const double lamb, const double *arrDM_params);

    static double estimate_aDM_and_dDM_and_LDM_(const QCavityAntenna cavityAntenna0, const double* ArrBaseData, const int quantBaseFreq
                                 , double*val_aDM, double*val_dDM, double*val_LDM);

    static double calNeviazV(  QCavityAntenna cavityAntenna,const double* ArrBaseData,const int quantBaseFreq);

};

#endif // QORIENTATION_H
