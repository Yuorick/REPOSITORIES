#ifndef SUBWATERBEAM_H
#define SUBWATERBEAM_H

class TTable_1D;
class QDataExchange;
class QBigMeasure;

enum enumTypeOf_000{ VAR0,VAR1};

class QSubWaterBeam
{
public:
    QSubWaterBeam();

    static int calcColCountFrom000(wchar_t*NameOf000File);

    static int ReadDataFrom000(wchar_t*NameOfDataFile, int *piNumRows, double *parrData);

    static int ReadDataFrom000_11_03_2022(wchar_t*NameOfDataFile, int *piNumRows, double *parrData);



    static bool replace(wchar_t*str);

    static bool replace(char*str) ;

    static void sortProfile(double *arrProfile, const int quantRows);

    static void adjustProfile(double *parrData, const int iNumRows, int &iNumRowsNew);

    static bool createProfileTbl(wchar_t*NameOfDataFile, TTable_1D *ptbl, enumTypeOf_000 TypeOf_000);


};

   double calcIntegral_I(const double VAlC0,const double VAlZn1
,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

   double fncSq1(const double VAlC0,const double VAlCn,const double VAlCosTetta);

   double calcIntegral_I1(const double VAlC0,const double VAlZn1
    ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

   double calcIntegral_I2(const double VAlC0,const double VAlZn1
    ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);

   double calcIntegral_I3(const double VAlC0,const double VAlZn1
    ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta);


   // вычисление горизоньтального расстояния луча от антенны
   double calcXHoriz(const double VAlza,const double VAlzm, TTable_1D &tblPrfl
                                        ,const double VAlCosTetta);

int fncCmp1( const void *a, const void *b);

double calcSumIntegral_I( const double VAlza,const double VAlzm
    , TTable_1D &tblPrfl,const double VAlCosTetta
    ,double (*f)(const double VAlC0,const double VAlZn1
    ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta));

// вычисление максимального значения скорости звука c(z)
// на отрезке [za;zm]
double calcMaxZ( const double VAlza,const double VAlzm, TTable_1D &tblPrfl);


//вычисление критического начального угла скольжения Tetta
//при котором горизонтальная координата луча на заданной
// глубине будет максимальна
double calcTettaCrit( const double VAlza,const double VAlzm, TTable_1D &tblPrfl);


// вычисление критического радиуса распространения звуковой волны
//на заданной глубине
double calcXMCrit( const double VAlza,const double VAlzm, TTable_1D &tblPrfl);

//вычисление Tetta
//Tetta- угол скольжения луча, исходящего из антенны
bool calcTetta( const double VAlza,const double VAlzm,const double VAlxm
               , TTable_1D &tblPrfl, double &valTetta);

double calc_dfi_po_dtetta( const double VAlza,const double VAlzm,const double tetta0
               , TTable_1D &tblPrfl);

void calc_gradTetta_po_dZ( const double VAlza,const double VAlzm
      ,const double tetta0, TTable_1D &tblPrfl, double *arrgrad);

double calc_t( const double VAlza,const double VAlzm
                                   ,const double tetta0, TTable_1D &tblPrfl);

//void calc_grad_t_po_dZ( const double VAlza,const double VAlzm
  //    ,const double tetta0, TTable_1D &tblPrfl,const double *ARrGradTetta, double *arrgrad);

double calc_CurveLength( const double VAlza,const double VAlzm
                                   ,const double tetta0, TTable_1D &tblPrfl);

double calc_CurTetta( const double VAlza,const double VAlzm
                                   ,const double tetta0, TTable_1D &tblPrfl);

bool calcTetta_InpGSK(TTable_1D &tblPrfl, const double *arrAntY,const double *arrTrueBeaconPos
               ,  double &valTetta);

bool calc_t_InpGSK(TTable_1D &tblPrfl, const double *arrAntY,const double *arrTrueBeaconPos
               ,  double &val_t);

void calc_dt_po_db( const double VAlza,const double VAlzm
      ,const double tetta0, TTable_1D &tblPrfl, double *arr_dt_po_db);



void calc_dZ_po_dY(const double *arrY,double * arr_dTetta_po_dZ,double *arr_dZ_po_dY);

void calc_dZ_po_dYAnt(const double *arrYAnt,const double *arrTrueBeaconPos, double *arr_dZ_po_dYAnt);

 bool calc_t_and_dt_po_dYAnt_InpGSK(TTable_1D &tblPrfl, const double *arrAntY
                ,const double *arrTrueBeaconPos,  double &val_t, double* arrdt_po_dYAnt);

bool calc_t_and_dt_po_dZ( const double VAlza,const double VAlzm,const double VAlxm
                          , TTable_1D &tblPrfl, double &valTetta, double *arr_dTetta_po_dZ
                          , double &val_t, double *arr_dt_po_dZ);

bool calc_t_and_dt_po_dX(TTable_1D &tblPrfl,  double *arrEilers0,  double *arrSZvAntPSK
                , double *arrSZvVessGSK,const double *arrSZvBeaconGSK
                         ,double &val_t, double* arrdt_po_dX/*, QDataExchange &DataExchange*/);

void calc_dZ_3Vars_po_dY_6Vars(const double *arrYAnt,const double *arrTrueBeaconPos, double *arr_dZ_po_dY);

void calc_dz_po_dY(const double *arrY,double *arr_dz_po_dY);

bool transfMeasure_to_GSK(TTable_1D &tblPrfl,QBigMeasure &BigMeasureInp
                          , double* arrAntPosParams,  double *arrSZv_GSK, double*valTZv);

//---------------------------------
// пересчет замера 3D в АСПК
bool transfMeasure_to_ASPK(TTable_1D &ptblEstPrfl,QBigMeasure &BigMeasureInp, double* arrAntPosParams,  double *arrSZv_ASPK, double*valTZv);

void calc_dZ_po_dSBeacon(const double *arrSAnt
                         ,const double *arrSBeacon, double *arr_dZ_po_dSBeacon);

void calc_dZ_po_dSBeacon(const double *arrSAnt
                         ,const double *arrSBeacon, double *arr_dZ_po_dSBeacon);

bool calc_TZapr_and_dTZapr_po_dTOtv(TTable_1D &tblPrfl, double* arrSAnt_GSK,double* arrSAntWave_GSK
               ,const double val_q_gsk, const double val_tetta
          ,const double val_tOtv, double *arrGskBeacon,double&val_tZapr,double& val_dtZapr_po_dtOtv);

bool calcBeamPosGSK_from_t_and_dBeamPosGSK_po_dt(TTable_1D &tblPrfl, double* arrSAnt_GSK
                                                 ,const double val_q_gsk, const double val_tetta
                                                      ,const double val_t,double *arrBeamGskPos ,double *arr_dBeamGskPos_po_dt);

bool calcDeepth_and_dDeepth_po_dt(TTable_1D &tblPrfl, const double  VAlZa, const double val_tetta
                                 ,const double val_t, double &val_z, double &val_dz_po_dt);

bool calc_t_and_dt_po_dSBeacon_InpGSK(TTable_1D &tblPrfl, const double *arrSAnt
                ,const double *arrSBeacon,  double &val_t, double* arrdt_po_dSBeacon);

bool transf_GSK_XYZ_to_USBL3D(TTable_1D &tblPrfl, double *arrGSK_XYZ, double *arrSVessGSK, double *arrMu
  ,double* arrAntPosParams,  double *pval_q, double *pval_e,  double *pval_t);

bool calc_tZapr(TTable_1D &tblPrfl, const double *arrSAnt
                ,const double *arrBeaconPos,const double *arrBeaconVelo,  double &val_t);

bool transfMeasure_to_GSK(TTable_1D &tblPrfl,QBigMeasure &BigMeasureInp
                          ,double* arrAntPosParams,const double *arrBeaconVelo
                          ,  double *arrSZv_GSK, double*valTZv);

bool calc_TZapr_and_dTZapr_po_dTOtv_(TTable_1D &tblPrfl, double* arrSAnt_GSK
                                     ,double* arrSAntWave_GSK,const double *arrBeaconVelo,const double VAlTobr
                                     , const double val_q_gsk, const double val_tetta
                                     ,const double val_tOtv, double *arrGskBeacon,double&val_tZapr,double& val_dtZapr_po_dtOtv);

bool calc_t_for_Vess_gsk(TTable_1D &tblPrfl,  double *arrEilers0,  double *arrSAntPSK
                , double *arrSVessGSK,const double *arrSBeaconGSK,double &val_t);

void calcSphericalAnglesKGSK(TTable_1D &tblPrfl, const double VAl_ASSK_q, const double VAl_ASSK_e, double *arrEilersMu
       ,double* arrAntPosParams, double *pval_KGSK_q, double *pval_KGSK_e);
#endif // SUBWATERBEAM_H
