#ifndef Solver2D_Angs_H
#define Solver2D_Angs_H
#include "PosSolver.h"

//class QDataExchange;
class QSolver2D_Angs : public QPosSolver
{
public:
    QSolver2D_Angs();

    QSolver2D_Angs (const  QSolver2D_Angs &R);

    QSolver2D_Angs  &operator=( const QSolver2D_Angs  &R);

    QSolver2D_Angs (const QBigMeasure *parrBigMeasures
                    , const int  QntMeas    ,const TTable_1D tblEstPrfl
                    , const double Deepth, const double Toler
                    ,const double *arrPAnt_PSK,const double *arrSBeaconXY_GSK);

    // парам конструктор 2
  /* QSolver2D_Angs (const QBigMeasure *parrBigMeasures
                   , const int  QntMeas    ,const TTable_1D tblEstPrfl
                   ,const double*arrGPSAprioriPosParams, const double Deepth
                   , const double Toler, const int ProgonNum, const int QuantIterMax
                   ,const double *arrPAnt_PSK ,const double *arrSBeaconXY_GSK);*/

    double marrPAnt_PSK[3];
    double marrSBeaconXY_GSK[2];

    QVector <double>mvct_e_psk ;

    void fill_vct_e_PSK();


    virtual void  collectGrlsMeasure(const QBigMeasure &Meas
                                                 ,double *arrMeasure,double *arrD);

    virtual bool calc_arrNeviaz_and_dArrNeviaz_po_dX(double *arrX, const int NUmBigMeasure
                           , double *arrNeviaz , double *dArrNeviaz_po_dX);



    virtual bool calc_arrNeviaz_and_dArrNeviaz_po_dX_3D( const int NUmBigMeasure
      ,const double e_est,double *arr_dV_po_dJ ,double *arrt7
      , double *arrNeviaz , double *dArrNeviaz_po_dX);

   virtual bool calc_arrNeviaz(double *arrX, const int NUmBigMeasure
                                               , double *arrNeviaz );


    void calc_arrI(double *arrZ, double *arrI);

    static double calcCorseNeviaz(const double qest, const double qzv);

  //  bool calc_2D_arrNeviaz_and_dArrNeviaz_po_dX_old(double *arrX, QBigMeasure Meas, QDataExchange DataExchange
    //                                       , double *arrNeviaz , double *dArrNeviaz_po_dX);

    void calc_arrZ(double *arrY, double *arrZ);


    static void calcMtrx_dJ_po_dX(double *arrX, double *arr_e_psk
                                                  , double *arr_dJ_po_dX);

    static void calcVectJ_(double *arrX,  double *arr_e_psk,double *arrJ);

    void calc_arr_e_PSK(QBigMeasure Meas, double *arr_e_psk);

    virtual bool calc_arrNeviaz_3D( const int NUmBigMeasure
                                                    ,const double e_est, double *arrNeviaz );

    //
    void calc_dAngs_po_dSgps(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                             , double *arr_dXpos_po_dSgps, QVector<double> *pvct_dAngs_po_dSgps);

    bool calc_dFGrA_po_dSgps(double *arrX0,double * parrMtrx_dFGrA_po_dSgps);

   // bool calc_dFGrA_po_dX0pos(double *arrAngX0,double * parrMtrx_dFGrA_po_dX0pos);

    bool  calc_arrdfiA_po_dXpos(double *arrX
                           ,const int NUmBigMeasure,double *arrdfi_po_dXpos);    

     void calc_dAngs_po_dHBeacon(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                 , double *arr_dX_po_dHBeacon,QVector<double> &vct_dAngs_po_dH);

    bool calc_dFGrA_po_dH(double *arrAngX0,double * parrMtrx_dFGrA_po_dSgps);

    bool  calc_arrdfiA_po_dH(double *arrX,const int NUmBigMeasure,double *arrdfi_po_dH);

    bool  calc_K_qe(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv,QVector<double> &vct_K_dqe);

    bool  calc_arrfi_and_arrdfi_po_dqe(double *arrAngX0
                                                ,const int NUmBigMeasure,double *arr_df_po_dqe);

    bool  calc_dAngs_po_dAngSins(double* arrAngX0,double* arrMtrx_dFGr_po_dXang_Inv
           , double *arr_dXPos_po_AngSins,QVector<double> &vct_dXa_po_AngSins);


    //bool  calc_dFGrA_po_dSins(double* arrAngX0,double*  parrMtrx_dFGrA_po_dSins);

   // bool  calc_arrdfiA_po_dAngSins(double *arrX
            //               ,const int NUmBigMeasure,double *arrdfi_po_dSins);

    bool  calc_dAngs_po_dProfile(double* arrAngX0,double* arrMtrx_dFGr_po_dXang_Inv
                                 ,double* arr_dXpos_po_dProfile,QVector<double> &vct_dAngs_po_dProfile);




    bool  calc_arrdfiA_po_dProfile(double *arrX
                           ,const int NUmBigMeasure,double *arrdfi_po_dProfile);

    bool  calc_AngKgps(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                       ,double *arrMtrx_dFGrLBL_po_dX_Inv,QVector<double> &vct__Kgps_per_1m);


    bool  calc_AngKsinsQ(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                         ,double *arrMtrx_dFGrLBL_po_dX_Inv ,QVector<double> &vct_AngK_SinsQ);

    bool  calc_dFusbl_po_dXlbl(double *arrAngX0,double *arrMtrx_dFusbl_po_dXlbl);

    bool  calc_AngKsins_Psi_Tetta(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                          ,double *arrMtrx_dFGrLBL_po_dX_Inv,QVector<double> &vct_AngKgps_Psi_Tetta_per_1_rad);


    bool  calc_AngKt(double *arrAngX0,double *arrMtrx_dFGr_po_dXang_Inv
                                     ,double *arrMtrx_dFGrLBL_po_dX_Inv ,QVector<double> &vct_AngKt_per_1sec);

};
#endif
