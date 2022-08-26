#ifndef CONTOBJ6_H
#define CONTOBJ6_H

// QUANT_GDG_PRMS - к-во переменных в задаче по КО вариантнт 1
// последовательность
// 0. XRelative
// 1. YRelative
// 2. ZRelative
// 3. BetRelative
// 4. Eps RLK
// 5. Eps IU
// 6. Alf RLK
// 7. Alf IU



#define  QUANT_GDG_PRMS 6
//#define LENX_V2 3

class QContObj6
{
public:
    QContObj6();

    //1 вектор Эйлеровых углов ЦКП, курсовой угол, угол килевой и бортовой качек
    double marrEilerCntrKP[3];

    //1 вектор измеренных Эйлеровых углов ЦКП, курсовой угол, угол килевой и бортовой качек
    double marrEilerCntrKP_Zv[3];

    //2 вектор дисперсий  ошибок оценивания Эйлеровых углов ЦКП
    double marrDispEilerCntrKP[3];

    // 3 вектор истинных параметров пространственной ориентации грани РЛК
    //  в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXTrue_RLK_PSK[QUANT_GDG_PRMS];

    //4 вектор истинных параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXTrue_IU_PSK[QUANT_GDG_PRMS];

    //5 вектор первичных оценок параметров пространственной ориентации грани РЛК
    // в ПСК -вектор параллакса, угол Betta, угол Eps
    double marrXZv_RLK_PSK[QUANT_GDG_PRMS];

    //6 вектор первичных оценок параметров пространственной ориентации основания ИУ
    // в ПСК - вектор параллакса, угол Betta, угол Eps
    double marrXZv_IU_PSK[QUANT_GDG_PRMS];

    //7 вектор истинных координат КО в КГСК
    double marrS_KGSK[3];

    //7 вектор измеренных координат КО в КГСК
    double marrSZv_KGSK[3];
    //8 Коррел матрица ошибок первичных измерений РЛК
    double marrCorMtrxRLK[9];
    //9 Коррел матрица ошибок первичных измерений КО в КГСК
  //  double marrCorMtrx_CO_KGSK[9];


    //10 СКЗ ошибки по углу ОЭК
    double mAngleSigIU;

    //11 СКЗ ошибки дальномера
    double mDistSigIU;

    //  СКЗ ошибки позиционирования корабля-носителя

    double mVessSKZ;
    //  СКЗ ошибки позиционирования КО

    double mCntrlObjSKZ;





//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

    QContObj6 (const QContObj6 &R);

    QContObj6 &operator=(const QContObj6  &R);

    // парм консьтруктор 1
    QContObj6(const double *arrEilerCntrKP, const double *arrDispEilerCntrKP
                ,const double * arrXTrue_RLK_PSK,const double * arrXTrue_IU_PSK
                ,const double *arrXZv_RLK_PSK,const double *arrXZv_IU_PSK
                ,const double * arrS_KGSK,const double * arrCorMtrx_CO_KGSK
                , const double AngleSigIU, const double DistSigIU
                              ,const double VessSKZ,const double  CntrlObjSKZ );
//--------------------------------------------------------------------------

void imitateMeasures(const bool BNoise,double *arrVZvRLK
                     , double *arrVZvIU);



static void calc_f(const bool BDalnomer,double *arrX_RLK,double *arrX_IU
                                , double *arrVZv_RLK, double *arrf);

//static void calc_Viu(double *arrXRelative_RLK,double *arrX_IU
//                     , double *arrVZv_RLK, double *arrViu);


static void calc_df_po_dXj(const bool BDalnomer,double *arrX_RLK,double *arrX_IU
                    , double *arrV_RLK,const int NUmj, double *arr_df_po_dXj);

//static void recalcPositionFromKGSK_to_SphericalSK(double *arrEilerCntrKP,double *arrS_KGSK,double *arrXGdg,double *arrVGdg);

static double estimateParams_Var_1(QContObj6 *pContObjectArr,const int QuantObj
                                       ,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst, double *arrMtrxCor);

static double imitate_and_estimateX_Var1(QContObj6 *pContObjectArr,const int QuantObj
     , const bool BNoise, double *arrXEst, double *arrMtrxCor);

static double calcSumNeviazka_Var1(double *arrXCur,QContObj6 *pContObjectArr
                 ,const int QuantObj ,double *parrVZv_RLK
                 ,double *parrVZv_IU);



static void calc_dViu_po_dXj(double *arrXRelative_RLK,double *arrX_IU
                             , double *arrV_RLK,const int NUmj, double *arr_dViu_po_dXj);

static void calc_dViu_po_dX_Razn(double *arrXRelative_RLK,double *arrX_IU
                                              , double *arrV_RLK, double *arr_dViu_po_dX);

double calcNeviazka_Var1(double *arrXCur,double *arrVZv_RLK, double *parrVZv_IU);

static double imitate_and_estimateAngles_Var1_(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise , double *arrXEst, double *arrMtrxCor);

static double estimateAngles_Var_1(QContObj6 *pContObjectArr,const int QuantObj
           ,double *parrVZv_RLK
           ,double *parrVZv_IU, double *arrXEst, double *arrMtrxCor);


static double doOneIterationAngles_Var1(double *arrX0,QContObj6 *pContObjectArr
                        ,const int QuantObj,double *parrVZv_RLK
                        ,double *parrVZv_IU,double *arrXRez, double *arrMtrx_dFGr_po_dX_Inv);

static bool calcCorMtrx_Var1(QContObj6 *pContObjectArr,const int QuantObj
                             , double *arrMtrx_dFGr_po_dX_Inv,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst, double *arrMtrxCor);

static  void  calc_dViu_po_dX_Var1_(double *arrXRelative_RLK,double *arrX_IU
                                                , double *arrVZv_RLK,double *arr_dViu_po_dX);

static bool doOneIteration_Var1_New(double *arrX0,QContObj6 *pContObjectArr
            ,const int QuantObj,double *parrVZv_RLK
            ,double *parrVZv_IU, double *arrXRez, double *arrMtrx_dF_po_dX_Inv);

static void  calc_FGr_and_dFGr_po_dX_Var1(double *arrX0,QContObj6 *pContObjectArr
                    ,const int QuantObj,double *parrVZv_RLK
                    ,double *parrVZv_IU,double *parrF,double *parrMtrx_dF_po_dX);

void  calc_fi_and_dfi_po_dx_Var1(double *arrX0, double *arrVZv_RLK
          ,double *arrVZv_IU,double *arrf,double *arrMtrx_df_po_dX);

static void calc_S_IUSPK(double *arrXRelative_RLK,double *arrX_IU
                         , double *arrVZv_RLK, double *arrS_IUSPK);

void calcMatr_dg_po_dXj(double *arrX0, double *arrVZv_RLK
                        ,const int NUmj,  double *arr_dg_po_dXj);

static void calc_RelativeEilers(double *arrEilers_RLK,double *arrEilers_IU
                                             ,double *arrEilers_Relative);

static void calc_RelativeParalacs(double *arrX6_RLK,double *arrX6_IU, double *arrX6_Relative);

static void calc_RelativePosition(double *arrX6_RLK,double *arrX6_IU
                                               , double *arrX6_Relative);

static void  calc_FGrAngs_and_dFGr_po_dX_Var1(double *arrX0,QContObj6 *pContObjectArr
                    ,const int QuantObj ,double *parrVZv_RLK
                    ,double *parrVZv_IU,double *parrF,double *parrMtrx_dF_po_dX);

void  calc_fiAng_and_dfiAng_po_dx_Var1(double *arrX0,double *arrVZv_RLK
          ,double *arrVZv_IU,double *arrfi,double *arrMtrx_dfi_po_dX);

static double imitate_and_estimateAngs_Var2(QContObj6 *pContObjectArr,const int QuantObj
                        , const bool BNoise, double *arrXEst, double *arrMtrxCor);

static double estimateAngs_Var_2(QContObj6 *pContObjectArr,const int QuantObj
                                                ,double *parrVZv_RLK, double *arrXRLKEst, double *arrMtrxCor);

static double calcSumNeviazka_Var2(double *arrXCur,QContObj6 *pContObjectArr
                                                ,const int QuantObj,double *parrVZv_RLK);

double calcNeviazka_Var2(double *arrX, double *arrVZv_RLK);

static bool doOneIterationAngs_Var2(double *arrX0,QContObj6 *pContObjectArr
                                             ,const int QuantObj,double *parrVZv_RLK
                                            , double *arrDel, double *arrMtrx_dFGr_po_dX_Inv);
static void  calc_FGrAngs_and_dFGrAngs_po_dX_Var2(double *arrX0,QContObj6 *pContObjectArr
               ,const int QuantObj ,double *parrVZv_RLK
               ,double *parrF,double *parrMtrx_dF_po_dX);

void  calc_fiAngs_and_dfiAngs_po_dx_Var2(double *arrX0, double *arrVZv_RLK
          ,double *arrfi,double *arrMtrx_dfi_po_dX);

static void calc_VGadget(double *arrEilerCntrKP,double *arrX0
                      , double *arrSZv_KGSK, double *arrf);

//void  calc_dVGdg_po_dX_Var2( double *arrX0
    //                      , double *arrS_KGSK,double *arr_g);

 static void  calc_dV_po_dX_Var2(double *arrX0, double *arrEilerCntrKP,double *arrS_KGSK,double *arr_g);

//void calcMatr_dg_po_dXj_Var2(double *arrX0, const bool BDalnomer ,double *arrSZv_KGSK
                //   ,const int NUmj,  double *arr_dg_po_dXj);


static bool calcCorMtrx_Var2(QContObj6 *pContObjectArr,const int QuantObj
             , double *arrXEst, double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor);

static double imitate_and_estimateParams_Var2(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor);

static double estimateParams_Var_2(QContObj6 *pContObjectArr,const int QuantObj
                                                ,double *parrVZv_IU, double *arrXIUEst, double *arrMtrxCor);

static bool doOneIteration_Var2(double *arrX0,QContObj6 *pContObjectArr
         ,const int QuantObj,double *parrVZv_IU, double *arrDel, double *arrMtrx_dFGr_po_dX_Inv);

static void calc_FGr_and_dFGr_po_dX_Var2(double *arrX0,QContObj6 *pContObjectArr
            ,const int QuantObj ,double *parrVZv_RLK
            ,double *parrF,double *parrMtrx_dF_po_dX);

void  calc_fi_and_dfi_po_dx_Var2(double *arrX0, double *arrVZv_RLK
          ,double *arrfi,double *arrMtrx_dfi_po_dX);

static double imitate_and_estimateStarProcessing(QContObj6 *pContObjectArr,const int QuantObj
       , const bool BNoise, double *arrXEst, double *arrMtrxCor);

void imitateStarMeasure(const bool BNoise, double *arrVZvIU);

static void recombine(double *arrXCombinedInp,double *arrXInp_IU, double *arrXOut_RLK,double *arrXOut_IU);

static void combine(double *arrXInp_RLK,double *arrXInp_IU,double *arrXOut);

static void recalcPositionFromKGSK_to_SphericalSK(double *arrEilerCntrKP0
         , double *arrS_KGSK0,double *arrXGdg0,double *arrVGdg0);

static void calc_Viu(double *arrXRLKCur0,double *arrX_IU0
                                 , double *arrVZv_RLK0, double *arrViu0);

 void calcMtrx_QKQT(double *parrVZv_IU, double *arrX,double *arrQKQT);

void calcMtrx_gTKg(double *parrVZv_RLK
              ,double *parrVZv_IU, double *arrX, double *arrgTgKggT);

static bool calcCorMtrxAngs_Var1(QContObj6 *pContObjectArr,const int QuantObj
                                            ,double *parrVZv_RLK,double *parrVZv_IU, double *arrXEst
                                            , double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor);

void calcMtrxAngs_QKQT(double *arrVZv_RLK, double *arrX,double *arrQKQT);

void calcMtrxAngs_gTKg(double *arrVZv_RLK
             ,double *arrVZv_IU, double *arrX, double *arrgTgKg);

void calcMtrx_JKJT(double *arrXEst, double *arrS_KGSK,double *arrJKJT_Cur);

void calcMtrx_CKCT(double *arrXEst,double *arrS_KGSK,double *arrCKCT_Cur);

void calcMtrx_g2TKg2(double *arrXEst,double *arrS_KGSK,double *arrg2TKg2_Cur);

static bool calcCorMtrxAngs_Var2(QContObj6 *pContObjectArr,const int QuantObj
                                            , double *arrXEst, double *arrMtrx_dFGr_po_dX_Inv, double *arrMtrxCor);

void calcMtrxAngs_JKJT(double *arrX0, double *arrS_KGSK,double *arrJKJT_Cur);

void calcMtrxAngs_CKCT(double *arrX0,double *arrS_KGSK,double *arrCKCT_Cur);

void calcMtrxAngs_g2TKg2(double *arrX0,double *arrS_KGSK,double *arrg2TKg2_Cur);

static void calc_dVGadget_po_dMu_Razn(double *arrX, double *arrMu
                                      , double *arrS_KGSK, double *arr_dViu_po_dMu);

static void calc_dVGadget_po_dMuj_Razn(double *arrX, double *arrMu, double *arrS_KGSK,const int NUmj, double *arr_dViu_po_dMuj);

static void calc_dVGadget_po_dMu(double *arrX, double *arrMu
                                            , double *arrS_KGSK, double *arr_dViu_po_dMu);

static void calc_IdealRelativeEilers_Var1(double *arrEilersTrue_RLK,double *arrEilersTrue_IU
                                                     ,double *arrEilersEst_IU ,double *arrEilersRlkIdeal);


};

void setDblArray(long double *arrInp, const int LEN, double *arrOut);
void setLongDblArray( double *arrInp, const int LEN,long double *arrOut);

#endif // CONTOBJ6_H
