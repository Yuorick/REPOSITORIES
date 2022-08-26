//---------------------------------------------------------------------------


#pragma hdrstop
#include  <string.h>
#include <math.h>
#include "Kirn_Vers1.h"
#include "Comp.h"
#include "MatrixProccess.h"

const double CONST_DIST_MODULE = 32.8;

const double    CONST_LAMBDA= 3.15;

const double CONST_DGR_WIDTH = 0.043589958111;

const double TET0707__ = 1.39115573;
/// Решение задачи ММП методом Ньютона
// Коэффициенты отражеия находятся путем решения сиситемы линейных уравнений
// После этого задача сводиться к решению сиситемы 2-х нелинейных уравнений
// относительно  palfTrg  и palfAnp
// Матрица Якоби рассчитывается разностным методом, для этогно по каждой переменной
// palfTrg  и palfAnp даются приращения и частные производные вычисляются разногстным методом
 //INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // palfTrg  и palfAnp - начальные пнриближения угла цели и антипода
 // OUTPUT:
 // pZTarg - комплексный коэфф  отражения цели
 // pZAnt -  комплексный коэфф  отражения антипода
 // palfTrg - угол цели
 // palfAnp -  угол антипода
 int   solvNewtonMeth_KirnStage1(const double valSigNoise,TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr )
 {

  int breturn = -3;
  int i = 0;
  double arr_FGreek[2] ={0.},arr_dFGreek[4] ={0.},arr_dFGreekInv[4] ={0.} ;
  double arrX[2] ={0.}, arrXT[2] ={0.}; // вектор с решениями
  arrX[0] = *palfTrg ;
  arrX[1] = *palfAnp;
  double del = -2.;

  for (i = 0; i < 10; i++)
  {

	double arr_F00[2] ={0.},arr_dF00[4] ={0.}, arrDelX[2] ={0.}, arrDelX1[2] ={0.} ;

	calc_vectG_and_mtrxH_NewtonMeth_Razn (pcmpSZv
		,  arrX[0], arrX[1], arr_FGreek,  arr_dFGreek, pZTarg, pZAnt   ) ;
		 //	 проверка положит определ
		double parrRez[4] ={0.}, parrRez1[4] ={0.};
		MatrxMultScalar(arr_dFGreek, 2, 2, 0.000001,parrRez);
		double valPos = arr_dFGreek[0] *arr_dFGreek[3] - arr_dFGreek[1] * arr_dFGreek[2];
		if (!((arr_dFGreek[0] > 0.)&&(valPos > 0.)))
		{
		 breturn = -3;
		 break;
		}
	///


	if(!InverseMtrx2(arr_dFGreek, arr_dFGreekInv))
	{
	 *palfTrg = arrX[0]   ;
	 *palfAnp = arrX[1]  ;
	 breturn = -2;
	 break;
	}

	MtrxMultMatrx(arr_dFGreekInv ,2, 2, arr_FGreek,1, arrDelX) ;
	del = NormVect(arrDelX, 2);
	double arrT[2] ={0.};
	MatrxMultScalar(arrDelX, 2, 1, 0.2,arrDelX1);
	memcpy( arrDelX, arrDelX1, 2 * sizeof(double));
	MtrxMinusMatrx(arrX, arrDelX,1, 2, arrXT) ;
	memcpy( arrX, arrXT, 2 * sizeof(double));
	if (del< 0.00001 )
	{

	  *palfTrg = arrX[0]   ;
	  *palfAnp = arrX[1]  ;

	 // calcMtrxCorrel_Mistake(valSigNoise, *palfTrg, *palfAnp, *pZTarg,  *pZAnt, arr_dFGreekInv, arrMtrxCorr);
	  breturn= 0;
	  break;
	}
 }

  return breturn;

}

//--------------------------------------------------------------------
// Вычисление вектор функции  левой части сиситемы уравнений F =0
// и матрицы частных производных dF/dx
// матрица dF/dx вычисляется разностным способом
// INPUT:
// pcmpSZv[4] - массив измерений дианрамм, комплексный
// alfTrg  и alfAnp -  угола цели и угол антипода
// OUTPUT:
// pcmpZTarg - комплексный коэфф  отражения цели
// pcmpZAnp -  комплексный коэфф  отражения антипода
//  arr_FGreek[2] - левая часть сиситемы уравнений
//  arr_dFGreek[4] - матрица частных производных левой части
int   calc_vectG_and_mtrxH_NewtonMeth_Razn (TComp *pcmpSZv, double alfTrg, double alfAnp
		,double* arr_FGreek, double*  arr_dFGreek , TComp*  pcmpZTarg,TComp* pcmpZAnp  )
 {

	 calc_vectF_from_Alfa__( pcmpSZv,  alfTrg, alfAnp,  arr_FGreek, pcmpZTarg, pcmpZAnp  );
	 double delt = 0.0000001;
	 double arr_FGreekTemp[2] ={0.}, arrT0[2] ={0.};
	 TComp cmpZTarg0(0.,0.), cmpZAnp0(0.,0.);
	 calc_vectF_from_Alfa__( pcmpSZv,  alfTrg + delt, alfAnp,  arr_FGreekTemp ,&cmpZTarg0, &cmpZAnp0  );
	 MtrxMinusMatrx(arr_FGreekTemp , arr_FGreek, 1, 2, arrT0);
	 MatrxDivideScalar(arrT0, 1, 2, delt,arr_dFGreek);

	 calc_vectF_from_Alfa__( pcmpSZv,  alfTrg , alfAnp + delt,  arr_FGreekTemp ,&cmpZTarg0, &cmpZAnp0  );
	 MtrxMinusMatrx(arr_FGreekTemp , arr_FGreek, 1, 2, arrT0);
	 MatrxDivideScalar(arrT0, 1, 2, delt,&arr_dFGreek[2]);
	 arr_dFGreek [2] =arr_dFGreek [1] ;
	 return 0;
 }

// НАхождение вектора  левой части сиситемы уравнений при фиксированных углах цели и антипода
// pcmpSZv[4] - массив измерений дианрамм, комплексный
// alfTrg  и alfAnp -  угола цели и угол антипода
// OUTPUT:
// pcmpZTarg - комплексный коэфф  отражения цели
// pcmpZAnp -  комплексный коэфф  отражения антипода
// arr_FGreek[2] - вектор левой части сиситемы уравнений
void  calc_vectF_from_Alfa__( TComp *pcmpSZv, double alfTrg, double alfAnp
  ,double* arr_FGreek, TComp*  pcmpZTarg,TComp* pcmpZAnp  )
{

  find_ZTarg_and_ZAnt__( pcmpSZv, pcmpZTarg, pcmpZAnp,  alfTrg,  alfAnp );

  double arr_dFGreek[4] ={0.}; // нигде не исполтьзуется
  calc_F__ ( pcmpSZv, *pcmpZTarg, *pcmpZAnp,  alfTrg,  alfAnp
	  ,  arr_FGreek );
}

// Нахождение оптимальных коэффициетов отражения  цели и антипода при фиксированных углах
// положения цели и антипода
// для вертикального веера
// INPUT:
// pcmpSZv[4]   - массив измерений по диаграммам (комплексные числа)
// alfTrg- угол места цели  в радианах
// alfAnp - угол места антипода в радианах

// OUTPUT:
// *pZTarg -  коэффициент отражения цели
// *pZAnt   - коэффициент отражения  антипода
// возвращает:
// -2 - если матрица якоби вырождена
// -3 - Если метод не сошелся
bool  find_ZTarg_and_ZAnt__( TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnp, double alfTrg, double alfAnp )
{
	// нахождение коэффициентов отражения при фиксированных углах
	// 2 задач квадрат программирования без ограничений 2-го порядка
	double arrZ[4] ={0.};
	double arrAATSum[16] ={0.}, arrAATCur[16] ={0.}, arrATSSum[4]={0.}, arrATSCur[4]={0.}, arrT0[16] ={0.}
	,arrT1[4]={0.}, arrAATSumInv[16]={0.};
	int numDiagr = 4;
	for (int i =0; i < numDiagr; i++)
	{

	 calc_ATA_andATS__(i, pcmpSZv[i], alfTrg,  alfAnp, arrAATCur, arrATSCur );
	 MtrxSumMatrx(arrAATSum, arrAATCur,4, 4, arrT0) ;
	 memcpy(arrAATSum,arrT0, 16 * sizeof(double));
	 MtrxSumMatrx(arrATSSum, arrATSCur,4, 1, arrT1) ;
	 memcpy(arrATSSum,arrT1, 4 * sizeof(double));
	}

	if(! InverseMtrx4(arrAATSum, arrAATSumInv)) return false ;
	MtrxMultMatrx(arrAATSumInv,4, 4, arrATSSum ,1, arrZ) ;
	(*pZTarg).m_Re = arrZ[0] ;
	(*pZTarg).m_Im = arrZ[1] ;
	(*pZAnp).m_Re = arrZ[2] ;
	(*pZAnp).m_Im = arrZ[3] ;
	return true;
}

//
void  calc_ATA_andATS__(int iNumDiagr, TComp cmpSZv, double alfTrg, double alfAnp , double*arrAATCur, double*arrATSCur )
{
   double arrA[8] = {0.}, arrA0[8] = {0.}, arrS[2] = {0.};
   TComp cmp0 = fncF__ (iNumDiagr, alfTrg);
   TComp cmp1 = fncF__ (iNumDiagr, alfAnp);
   arrA[0] =  cmp0.m_Re;
   arrA[1] =  -cmp0.m_Im;
   arrA[2] =  cmp1.m_Re;
   arrA[3] =  -cmp1.m_Im;
   arrA[4] =  cmp0.m_Im;
   arrA[5] =  cmp0.m_Re;
   arrA[6] =  cmp1.m_Im;
   arrA[7] =  cmp1.m_Re;
   memcpy(arrA0, arrA, 8 * sizeof(double));

   arrS[0] =  cmpSZv.m_Re;
   arrS[1] =  cmpSZv.m_Im;

   MtrxTranspMultMatrx(arrA,2, 4, arrA0,4, arrAATCur) ;
   MtrxTranspMultMatrx(arrA,2, 4, arrS,1, arrATSCur) ;

}


void  calc_F__ ( TComp *pcmpSZv, TComp ZTarg, TComp ZAnt
, double alfTrg, double alfAnp, double* arr_FGreek )
{
	double arr_Part_F[2] ={0.},arr_Part_dF[4] ={0.}, arrT[2]={0.}, arrT1[4]={0.};
	memset(arr_FGreek, 0, 2 * sizeof(double));

	int NumDiagr = 4 ;
	for (int i =0; i < NumDiagr; i++)
	{
	calcPartial_F__(i,pcmpSZv[i], ZTarg,  ZAnt
		, alfTrg, alfAnp , arr_Part_F ) ;
	 MtrxSumMatrx(arr_FGreek, arr_Part_F,1, 2, arrT) ;
	 memcpy(arr_FGreek,arrT, 2 * sizeof(double));


	}
}

//--------------------------------------------------------
void  calcPartial_F__(int iNumDiagr, TComp cmpSZv,TComp ZTarg, TComp ZAnp
		,double alfTrg,double alfAnp ,double* arr_Part_F )
{
  TComp cmpMesTarg = fncF__ (iNumDiagr,alfTrg);;
  TComp cmpMesAnp = fncF__ (iNumDiagr, alfAnp);;
  TComp cmpA = ZTarg * cmpMesTarg + ZAnp * cmpMesAnp  -cmpSZv;

  TComp cmp_dF1_po_dx = dF_po_dTet__ (iNumDiagr, alfTrg);
  TComp cmp_dF2_po_dx = dF_po_dTet__ (iNumDiagr, alfAnp);



  double val_dA1_po_dx1 =  ZTarg.m_Re * cmp_dF1_po_dx.m_Re -  ZTarg.m_Im * cmp_dF1_po_dx.m_Im ;
  double val_dA1_po_dx2 =  ZAnp.m_Re * cmp_dF2_po_dx.m_Re -  ZAnp.m_Im * cmp_dF2_po_dx.m_Im ;

  double val_dA2_po_dx1 =  ZTarg.m_Im * cmp_dF1_po_dx.m_Re +  ZTarg.m_Re * cmp_dF1_po_dx.m_Im ;
  double val_dA2_po_dx2 =  ZAnp.m_Im * cmp_dF2_po_dx.m_Re +  ZAnp.m_Re * cmp_dF2_po_dx.m_Im ;
  // вектор функция
  arr_Part_F[0] = (ZTarg.m_Re * cmpA.m_Re + ZTarg.m_Im * cmpA.m_Im) * cmp_dF1_po_dx.m_Re
		   +(-ZTarg.m_Im * cmpA.m_Re + ZTarg.m_Re * cmpA.m_Im ) * cmp_dF1_po_dx.m_Im ;

  arr_Part_F[1] = (ZAnp.m_Re * cmpA.m_Re + ZAnp.m_Im * cmpA.m_Im) * cmp_dF2_po_dx.m_Re
		   +(-ZAnp.m_Im * cmpA.m_Re + ZAnp.m_Re * cmpA.m_Im ) * cmp_dF2_po_dx.m_Im ;

}

TComp fncF__ (int iNumDiagr, const double tet)
{
  TComp cmpRez(0.,0.) ;
  TComp delFas(0.,  CONST_DIST_MODULE * sin(tet)/ 2./ CONST_LAMBDA * 2. * M_PI);

  TComp delFas0(0.,  -delFas.m_Im - delFas.m_Im * 2.);
  TComp cmpDia(fncDiagrSimple__(CONST_DGR_WIDTH, tet),0.);
  TComp cmpFasCur ( 0.,((double)iNumDiagr) * 2. * delFas.m_Im+delFas0.m_Im);
  TComp cmpTemp = exp_(cmpFasCur);
  cmpRez =  cmpDia* cmpTemp;
  return  cmpRez;
}

// диаграмма sin(ax)/ax
double fncDiagrSimple__(double valWidthDgr, double valTetRad)
{
double coeff =  TET0707__ / valWidthDgr;
return fncDiagrSinx_div_x__(valTetRad * coeff);
}

double fncDiagrSinx_div_x__(double tet)
{
if (fabs(tet)< 0.0000001) return 1.;

return sin(tet)/tet;
}


TComp dF_po_dTet__(int iNumDiagr, const double tet)
{
	double  delFas =  CONST_DIST_MODULE * sin(tet)/ 2./ CONST_LAMBDA * 2. * M_PI;
	double  d_delFas =  CONST_DIST_MODULE * cos(tet)/ 2./ CONST_LAMBDA * 2. * M_PI;
	double temp = ((double)( -  4 + 1 + 2 * iNumDiagr));
	TComp cmpTemp = exp_(delFas * temp );
	TComp cmp025(0.25,0);
	cmpTemp *= cmp025;

	double valS =  fncDiagrSimple__(CONST_DGR_WIDTH, tet);
	double valdS = fncDerivDiagrSimple__(CONST_DGR_WIDTH, tet);
	TComp cmpTemp1(valdS ,  temp * d_delFas * valS);
	TComp cmpTemp2 = cmpTemp * cmpTemp1;

return  cmpTemp2 ;
}


// производная диаграмма sin(ax)/ax
double fncDerivDiagrSimple__(double valWidthDgr, double valTetRad)
{
double coeff =  TET0707__ / valWidthDgr;
return coeff* fncDerivDiagrSinx_div_x__(coeff *valTetRad);
}

double fncDerivDiagrSinx_div_x__(double tet)
{
if (fabs(tet)< 0.0000001) return 0.;
return (tet * cos(tet) -sin(tet))/ tet/tet;
}
#pragma package(smart_init)

