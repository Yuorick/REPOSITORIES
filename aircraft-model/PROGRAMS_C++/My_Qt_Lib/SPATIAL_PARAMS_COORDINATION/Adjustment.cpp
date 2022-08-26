#include <string.h>
#include "Adjustment.h"
#include "CalcCorMatrx.h"
#include "MatrixProccess.h"


QAdjustment::QAdjustment()
{

}

//Вычисление координат цели в сферической сиситеме координат прибора
//INPUT:
//arrS_KGSK [3] - координаты цели в КГСК
//arrSins[3] - углы СИНС- КУ, КК, БК
//arrX[NUM_GADG_PARAMS] - вектор параметров позиуционирования прибора в ПСК
//OUTPUT:
//arrVIU[2] -углы гориз и вертикальногот наведения ПУГН и ПУВН
void QAdjustment::calc_VIU(double * arrS_KGSK ,double * arrSins
            ,double *arrXGdg,  double &valR, double &valBet, double &valEps)
{

    // 1. истинное положение в ПСК
    double matrPereh_PSK_V_KGSK[9] = {0.};
    calcMatr_PSK_v_KGSK(arrSins,matrPereh_PSK_V_KGSK) ;
    double arrSTrue_PSK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_PSK_V_KGSK,3, 3, arrS_KGSK,1, arrSTrue_PSK) ;
    ///
    // 2. истинное положение в ПСК-УСТРОЙСТВО
    double arrSTrue_PSK_Gdg[3] = {0.};
    MtrxMinusMatrx(arrSTrue_PSK, arrXGdg,1, 3, arrSTrue_PSK_Gdg);
    ///

    // 3. истинное положение  в АСПК
    double matrPereh_ASK_V_PSK[9]= {0.};

    calcMtrx3_ASPK_v_PSK(&(arrXGdg[3]),matrPereh_ASK_V_PSK) ;


     ///

    double arrSTrue_ASPK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrSTrue_PSK_Gdg,1, arrSTrue_ASPK) ;
    ///

    // 4. положение в АСфСАК
    double temp = 0.;
    recalcCoord_INTO_Spherical(arrSTrue_ASPK, valR, valBet, valEps);
    ///

}

//Вычисление Якобиана сферических координат цели
// по параметрам позиционирования разностным методом
//INPUT:
//arrTargS_KGSK [3] - координаты цели в КГСК
//arrSins[3] - углы СИНС- КУ, КК, БК
//arrX[NUM_GADG_PARAMS] - вектор параметров позиционирования прибора в ПСК
//OUTPUT:
//arrVIU[2] -углы гориз и вертикальногот наведения ПУГН и ПУВН
// arr_dV_po_dX[2 * NUM_GADG_PARAMS] - якобиан
void QAdjustment::calc_VIU_and_dVIU_po_dX(double * arrS_KGSK ,double * arrSins
            ,double *arrXGdg, double &valR, double &valBet, double &valEps, double *arr_dVIU_po_dX)
{
    calc_VIU(arrS_KGSK , arrSins, arrXGdg, valR, valBet, valEps);
    for(int j =0; j < NUM_GADG_PARAMS; ++j)
    {
        double arrT  [2] = {0.};

        calc_dVIU_po_dXj(arrS_KGSK , arrSins,arrXGdg, j, arrT ) ;
        changeCol(arr_dVIU_po_dX, 2 ,  NUM_GADG_PARAMS,j ,arrT);
    }

}

//-----------------------------------------

//Вычисление координат цели в сферической сиситеме координат прибора
//INPUT:
//arrTargS_KGSK [3] - координаты цели в КГСК
//arrSins[3] - углы СИНС- КУ, КК, БК
//arrX[NUM_GADG_PARAMS] - вектор параметров позиуционирования прибора в ПСК
//OUTPUT:
//arrVIU[2] -углы гориз и вертикальногот наведения ПУГН и ПУВН
void QAdjustment::recalc_VIU_from_RLK_to_IU(double * arrVZv_RLK
            ,double *arrX_RLK,double *arrX_IU, double &valBet, double &valEps)
{

    ///
    // 1.  положение в АСПК
    double arrS_ASPK[3] = {0.};
    recalcSphericalCoord_INTO_Rectangular(arrVZv_RLK[1],arrVZv_RLK[0]
                                               ,arrVZv_RLK[2], arrS_ASPK);
    ///

    // 2.  положение  в ПСК-РЛК
    double matrPereh_ASK_V_PSK[9]= {0.};    
    calcMtrx3_ASPK_v_PSK(&(arrX_RLK[3]),matrPereh_ASK_V_PSK) ;


    double arrS_PSK_RLK[3] = {0.};
    MtrxMultMatrx(matrPereh_ASK_V_PSK,3, 3, arrS_ASPK,1, arrS_PSK_RLK) ;
    ///

    // 3. положение в ПСК
    double arrS_PSK[3] = {0.};
    MtrxSumMatrx(arrS_PSK_RLK, arrX_RLK,1, 3, arrS_PSK) ;

    // 4. положение в ПСК-ИУ
    double arrS_PSK_IU[3] = {0.};
    MtrxMinusMatrx(arrS_PSK, arrX_IU,1, 3, arrS_PSK_IU) ;
    ///

    // 5. положение в ИУСПК
    double matrPereh_ASK_V_PSK1[9]= {0.};
    calcMtrx3_ASPK_v_PSK(&(arrX_IU[3]),matrPereh_ASK_V_PSK1) ;


    double arrS_IUSPK[3] = {0.};
    MtrxTranspMultMatrx(matrPereh_ASK_V_PSK1,3, 3, arrS_PSK_IU,1, arrS_IUSPK) ;
    ///

    // 6.  положение в СфСК-ИУ

    double temp = 0;
    recalcCoord_INTO_Spherical(arrS_IUSPK, temp, valBet, valEps);
    ///

}

//Вычисление Якобиана сферических координат цели
// по параметрам позиционирования
//INPUT:
//arrTargS_KGSK [3] - координаты цели в КГСК
//arrSins[3] - углы СИНС- КУ, КК, БК
//arrX[NUM_GADG_PARAMS] - вектор параметров позиционирования прибора в ПСК
//OUTPUT:
//arrVIU[2] -углы гориз и вертикальногот наведения ПУГН и ПУВН
// arr_dVIU_po_dX[2 * NUM_GADG_PARAMS] - якобиан
void QAdjustment::recalc_VIU_from_RLK_to_IU_and_dVIU_po_dX(double * arrVZv_RLK
            ,double *arrX_RLK,double *arrX_IU,  double &valBet, double &valEps, double *arr_dVIU_po_dX)
{
    recalc_VIU_from_RLK_to_IU(arrVZv_RLK,arrX_RLK,arrX_IU, valBet, valEps);
    for(int j =0; j < LENX_V1; ++j)
    {
        double arrT  [2] = {0.};

        calc_dVIU_po_dXj_V1(arrVZv_RLK ,arrX_RLK,arrX_IU, j, arrT ) ;
        changeCol(arr_dVIU_po_dX, 2 ,  LENX_V1,j ,arrT);
    }


}

//-----------------------------------------
// вычисление разностной производной функции f по переменной
// с номером j
//функция f зависит от 7 переменных
//от переменных вектора позиционирования РЛК и ИУ  в ПСК - arrX

//OUTPUT:
//arr_df_po_dViu[2 + BDalnomer] - вектор частных производных ыектор функции Viu
//по переменной с номером NUmj (нумерация с 0)
void QAdjustment::calc_dVIU_po_dXj(double *arrS_KGSK ,double * arrSins,  double *arrXGdg
                         ,const int NUmj, double *arr_dVIU_po_dXj)
{


    double arrXGdg1[NUM_GADG_PARAMS] = {0.};
    memcpy(arrXGdg1, arrXGdg, NUM_GADG_PARAMS * sizeof(double));
    double delta = 0.;
    if (NUmj < 3)
    {
     delta = 0.1;

    }
    else
    {
        delta = 0.0001;


    }
    arrXGdg1[NUmj]+= delta;
    double arrVIU[2] = {0.},arrVIU1[2] = {0.}, arrT[2] = {0.},arrDelta[2] = {0.};
    double temp;
    calc_VIU( arrS_KGSK ,arrSins ,arrXGdg,temp, arrVIU[0], arrVIU[1]);
    calc_VIU( arrS_KGSK ,arrSins ,arrXGdg1, temp, arrVIU1[0], arrVIU1[1]);

    MtrxMinusMatrx(arrVIU1, arrVIU,1, 2, arrT);
    MatrxMultScalar(arrT, 1, 2, 1./delta,arrDelta);

    for (int i =0; i < 2;++i)
    {
      arr_dVIU_po_dXj[i] =  arrDelta[i] ;

    }
}

void QAdjustment::calc_dVIU_po_dXj_V1(double *arrVZv_RLK,double *arrX_RLK
                ,double *arrX_IU, const int NUmj, double *arr_dVIU_po_dXj)
{


    double arrX_RLK1[NUM_GADG_PARAMS] = {0.};
    memcpy(arrX_RLK1, arrX_RLK, NUM_GADG_PARAMS * sizeof(double));


    double arrX_IU1[NUM_GADG_PARAMS] = {0.};
    memcpy(arrX_IU1, arrX_IU, NUM_GADG_PARAMS * sizeof(double));
    double delta = 0.;
    if (NUmj <3)
    {
     delta = 0.1;
     arrX_RLK1[NUmj]+= delta;
    }
    else
    {

    delta = 0.0001;
    switch(NUmj)
    {
    case 3:// betta RLK
      arrX_RLK1[3]+= delta;
        break;
    case 4: // betta IU
        arrX_IU1[3]+= delta;
        break;
    case 5:  // Eps RLK
        arrX_RLK1[4]+= delta;
        break;
    case 6: // Eps IU
        arrX_IU1[4]+= delta;

        break;
    case 7:  //Alf RLK
        arrX_RLK1[5]+= delta;
        break;

    case 8:  //Alf IU
        arrX_IU1[5]+= delta;
        break;

    default:
        break;
    }
    }


    double arrVIU[2] = {0.},arrVIU1[2] = {0.}, arrT[2] = {0.},arrDelta[2] = {0.};
    double temp;
    recalc_VIU_from_RLK_to_IU(arrVZv_RLK,arrX_RLK,arrX_IU, arrVIU[0], arrVIU[1]);
    recalc_VIU_from_RLK_to_IU(arrVZv_RLK,arrX_RLK1,arrX_IU1, arrVIU1[0], arrVIU1[1]);

    MtrxMinusMatrx(arrVIU1, arrVIU,1, 2, arrT);
    MatrxMultScalar(arrT, 1, 2, 1./delta,arrDelta);

    for (int i =0; i < 2;++i)
    {
      arr_dVIU_po_dXj[i] =  arrDelta[i] ;

    }

}

