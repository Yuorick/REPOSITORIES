//---------------------------------------------------------------------------


#pragma hdrstop

#include "Article_Rad_El_2016.h"
#include "Far.h"
#include "Faceta.h"
#include <vcl.h>
#include <float.h>
#include <math.h>
#include "Comp.h"
#include <stdio.h>
#include <stdlib.h>
#include "YrWriteShapeFile.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include "Equations.h"
#include "MatrixProccess.h"
#include "Gauss.h"
#include "Diagr.h"
#include "URPolygon.h"
#include "YrWrite.h"


// для статьи Р и Э 2016
// построенгие графиков метода статистических испытаний для фиксированных значений
// угла цели и антипода
// сравнение квадратичного метода с методом миаксимального правдоподобия
// задача минимизации изметода максимального правдоподобия решается перебором на сетке
//wchFoldName  - путь к папаеке с отчетом
// NIsp  - число испытаний
// alfUMTrg, cmpKTarg - ум цели, коеффиц отражения цели
// alfUMAntp, cmpKAntp - ум антип, коеффиц отражения антип
void makeMonteCarloGraphs(TFar Far, wchar_t *wchFoldName1, const int NIsp, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp)

{
 if (Far.m_N != 4)
  {
   ShowMessage(L"ERROR: Far.m_N != 4");
   return;
  }
 int numPoints = 0;
  TComp *pcmpSZv = (TComp *)malloc( Far.m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( Far.m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[9] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);

  int iNCols = 20;
  double * parrBuff = (double *) malloc(iNCols * NIsp * sizeof(double));
  TComp z1(0.,0.),z2(0.,0.);
  TURPointXY *ppntarrZ1Izm = new TURPointXY[ NIsp];
  TURPointXY *ppntarrZ2Izm = new TURPointXY[ NIsp];
  double arrEstAng[3] = {0.};
  TComp cmparrK[3];
  TComp cmparrZ[3];


  const int iQUANT_STEPS = 1000;
  for (int i = 0; i < NIsp; i++)
	 {
		Far.ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);
		double valXiSquare= 10000000.;
		Far.fncPolyn2(pcmpSZv,  arrEstAng, cmparrK, arrMtrxCorr , cmparrZ, &valXiSquare);
		z1 = cmparrZ[0];
		z2 = cmparrZ[1];
		valEstAngTarg =  arrEstAng[0];
		valEstAngAntp =  arrEstAng[1];
		cmpEstKTarg =  cmparrK[0];
		cmpEstKAntp  =  cmparrK[1];
	   double valXi2  =0., valEstRSM =0., valRSMDisp0 =0., valRSMDisp =0.;



	  TComp cmpZTargMMP, cmpZAntMMP;
	  double valalfTrgMMP, valalfAnpMMP, valMinObj;
	  double arrMtrxCorrMMP[4];
	  double *parrObj = new double  [iQUANT_STEPS * iQUANT_STEPS];
	  solveMinMMP_Perebor_Na_Setke(Far ,  iQUANT_STEPS, pcmpSZv, &cmpZTargMMP, &cmpZAntMMP
	   , &valalfTrgMMP, &valalfAnpMMP, arrMtrxCorrMMP, parrObj , &valMinObj ) ;
	   if ((arrMtrxCorrMMP[0] < 0.) && (arrMtrxCorrMMP[3] < 0.))
	   {
		  ShowMessage(L"ERROR (arrMtrxCorrMMP[0] < 0.) && (arrMtrxCorrMMP[3] < 0.)");
		  return;
	   }
	  parrBuff[ i * iNCols + 1] = alfUMTrg * 1000.;
	  parrBuff[ i * iNCols + 2] =  alfUMAntp * 1000.;

	  parrBuff[ i * iNCols + 3] = valEstAngTarg * 1000.;
	  parrBuff[ i * iNCols + 4] = valEstAngAntp * 1000.;

	  parrBuff[ i * iNCols + 5] = (valEstAngTarg  - alfUMTrg )* 1000.;
	  parrBuff[ i * iNCols + 6] = (valEstAngAntp - alfUMAntp)* 1000.;
	  parrBuff[ i * iNCols + 7] = sqrt(arrMtrxCorr[0])* 1000. *3.;
	  parrBuff[ i * iNCols + 8] = -parrBuff[ i * iNCols + 7] ;
	  parrBuff[ i * iNCols + 9] =  sqrt(arrMtrxCorr[3])* 1000. *3.;
	  parrBuff[ i * iNCols +10] = -parrBuff[ i * iNCols + 9];

	  parrBuff[ i * iNCols] = (double)(i);
	  parrBuff[ i * iNCols + 11] = valXi2;

	  //для ММП
	  parrBuff[ i * iNCols + 12] = valalfTrgMMP * 1000.;
	  parrBuff[ i * iNCols + 13] = valalfAnpMMP * 1000.;

	  parrBuff[ i * iNCols + 14] = (valalfTrgMMP  - alfUMTrg )* 1000.;
	  parrBuff[ i * iNCols + 15] = (valalfAnpMMP - alfUMAntp)* 1000.;
	  parrBuff[ i * iNCols + 16] = sqrt(arrMtrxCorrMMP[0])* 1000. *3.;
	  parrBuff[ i * iNCols + 17] = -parrBuff[ i * iNCols + 16] ;
	  parrBuff[ i * iNCols + 18] =  sqrt(arrMtrxCorrMMP[3])* 1000. *3.;
	  parrBuff[ i * iNCols + 19] = -parrBuff[ i * iNCols + 18];

	  ppntarrZ1Izm[i].X =  z1.m_Re;
	  ppntarrZ1Izm[i].Y =  z1.m_Im;
	  ppntarrZ2Izm[i].X =  z2.m_Re;
	  ppntarrZ2Izm[i].Y =  z2.m_Im;

	  delete parrObj;


	 }
	 free(pcmpSZv);
	 free(pcmpS);

	 wchar_t *wcharrFileNames = new wchar_t[iNCols * 30];
	 memset(wcharrFileNames, 0,iNCols * 30* sizeof (wchar_t));
	 wcscpy(&wcharrFileNames[0],L"n");
	 wcscpy(&wcharrFileNames[30],L"plnUMTrgTrue");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnUMAntpTrue");
	 wcscpy(&wcharrFileNames[30 * 3],L"plnUMTrgIzm");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnUMAntpIzm");
	 wcscpy(&wcharrFileNames[30 * 5],L"plnUMTrgErr");
	 wcscpy(&wcharrFileNames[30 * 6],L"plnUMAntpErr");
	 wcscpy(&wcharrFileNames[30 * 7],L"plnUMTrg_3Sig");
	 wcscpy(&wcharrFileNames[30 * 8],L"plnUMTrg_Minuis3Sig");
	 wcscpy(&wcharrFileNames[30 * 9],L"plnUMAntp_3Sig");
	 wcscpy(&wcharrFileNames[30 * 10],L"plnUMAntp_Minuis3Sig");
	 wcscpy(&wcharrFileNames[30 * 11],L"Xi2");

	 wcscpy(&wcharrFileNames[30 * 12],L"plnUMTrgIzm_MMP");
	 wcscpy(&wcharrFileNames[30 * 13],L"plnUMAntpIzm_MMP");
	 wcscpy(&wcharrFileNames[30 * 14],L"plnUMTrgErr_MMP");
	 wcscpy(&wcharrFileNames[30 * 15],L"plnUMAntpErr_MMP");
	 wcscpy(&wcharrFileNames[30 * 16],L"plnUMTrg_3Sig_MMP");
	 wcscpy(&wcharrFileNames[30 * 17],L"plnUMTrg_Minuis3Sig_MMP");
	 wcscpy(&wcharrFileNames[30 * 18],L"plnUMAntp_3Sig_MMP");
	 wcscpy(&wcharrFileNames[30 * 19],L"plnUMAntp_Minuis3Sig_MMP");



	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i < iNCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,iNCols // - к-во переменных о корорых накоплена информация в буфере
								  ,NIsp //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , scalex  //  масштаб по оси X
								  , scaley  // масштаб по оси Y
								   ) ;
	 }

free(parrBuff);

	///


	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., double (NIsp) + 3.
	 ,-30.,30., 2.) ;
	// графики z1 и z2 измерений
  wchar_t wchZFileName[300] ={0};
  wcscpy(  wchZFileName,  wchFoldName);
  wcscat(wchZFileName, L"ZIzm1.shp");
  TURPointXY::WriteSetSHPFiles(wchZFileName,ppntarrZ1Izm, NIsp ) ;
  wcscpy(  wchZFileName,  wchFoldName);
  wcscat(wchZFileName, L"ZIzm2.shp");
  TURPointXY::WriteSetSHPFiles(wchZFileName,ppntarrZ2Izm, NIsp ) ;




delete []ppntarrZ1Izm;
delete []ppntarrZ2Izm;
 delete wcharrFileNames;
const double VAL_WAVE_CONST = 2. * M_PI * Far.m_D / Far.mFaceta.mLambda;
// построение  окружности с центром в 0
TURPointXY pointZero(0.,0.);
TURPolygon plgCircle0 = TURPolygon::fncCreateCircle(pointZero, 1., 500);
wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"Circle0.shp");
plgCircle0.WriteSetSHPFiles(wchAxesFileName,&plgCircle0, 1) ;
	 ///
	 Far.ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

		Far.fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	  , &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr, &z1, &z2 );
// посмтроение окружности  с центром в точке z1
TURPointXY pointZ1(z1.m_Re,z1.m_Im);
TURPolygon plgCircle1 = TURPolygon::fncCreateCircle(pointZ1, 3. * sqrt(arrMtrxCorr[0])*VAL_WAVE_CONST, 500);
wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"Circle1.shp");
plgCircle1.WriteSetSHPFiles(wchAxesFileName,&plgCircle1, 1) ;

wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"z1.shp");
pointZ1.WriteSetSHPFiles(wchAxesFileName,&pointZ1, 1) ;
///


// посмтроение окружности  с центром в точке z2
TURPointXY pointZ2(z2.m_Re,z2.m_Im);
TURPolygon plgCircle2 = TURPolygon::fncCreateCircle(pointZ2, 3. * sqrt(arrMtrxCorr[3])*VAL_WAVE_CONST, 500);
wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"Circle2.shp");
plgCircle2.WriteSetSHPFiles(wchAxesFileName,&plgCircle2, 1) ;

wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"z2.shp");
pointZ2.WriteSetSHPFiles(wchAxesFileName,&pointZ2, 1) ;
///
	///



	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-1.5, 1.5
	 ,-1.5, 1.5, 0.2) ;



}


/// Решение задачи ММП методом  перебора  на сетке iQUANT_STEPSxiQUANT_STEPS по углам цели и антипода
// для 4-х строковых диаграмм
 //INPUT:
 // iQUANT_STEPS - число разбиений по каждой переменной
 // pcmpSZv - массив измерений дианрамм, комплексный
 // valSigNoise - скз шума
 // OUTPUT:
 // pZTarg - комплексный коэфф  отражения цели
 // pZAnt -  комплексный коэфф  отражения антипода
 // palfTrg - угол цели
 // palfAnp -  угол антипода
 // arrMtrxCorr[4] - кореляц матрица ошибок оценивания углов цели и антипода
 // parrObj[iQUANT_STEPS *iQUANT_STEPS] - массив значений целевой функции
 // *pvalMinObj - минимум функции правдоподобия
 int   solveMinMMP_Perebor_Na_Setke(const TFar Far
  , const int iQUANT_STEPS, TComp *pcmpSZv, TComp *pZTarg, TComp *pZAnt
	   , double *palfTrg, double *palfAnp, double *arrMtrxCorr, double *parrObj , double *pvalMinObj )
 {
  if(Far.m_N != 4) return -3;

  const double valSigNoise =  sqrt( Far.mparrDisp[0]*4. );
  // полная диграмма АМ в одну стророну (ноль диаграммы)
  const double VAL_Tet0 = Far.mLambda /Far.m_N/ Far.m_D;

  const double VALNu =  2.* M_PI* Far.m_D / Far.mLambda ; // констатнта
  const double VALDiapazon = VALNu * sin(VAL_Tet0);  // граница изменения обобщ угла mu
  const double VALStep = 2. * VALDiapazon / ((double)iQUANT_STEPS - 1.);// шаг по обобщенным углам
  double valMu0 = -VALDiapazon; // обоб угол цели
  double valMu1 = -VALDiapazon; // обоб угол антипода

  TComp cmparrV[8], cmparrVH[8], cmparrTemp[4],cmparrTempInv[4], cmparrA[2], cmpZTargTemp(0.,0.), cmpZAntTemp(0.,0.);

  double alfTrgTemp = 0., alfAnpTemp = 0., arrMtrxCorrTemp[4] ={0.};
  double objMin = 1000000000.;

  for ( int i = 0; i < iQUANT_STEPS; i++)
  for ( int j = 0; j < iQUANT_STEPS; j++)
  {
	valMu0 = -VALDiapazon + ((double)i) * VALStep;
	valMu1 = -VALDiapazon - VALStep/ 100. + ((double)j) * VALStep;
	fncCreateMtrxV(valMu0, valMu1, cmparrV);

	HermiteSoprMatr(cmparrV, 4, 2, cmparrVH) ;
	MtrxMultMatrx(cmparrVH,2, 4, cmparrV,2, cmparrTemp);
	bool brez = InverseMtrx2( cmparrTemp, cmparrTempInv) ;
	if (!brez)
	{
	   TComp cmpVectV0[4], cmpVectV1[4], cmparrX[2], cmparrB[2]
	, cmparrProjection[4],cmparrPerp[4], cmparrGram[4], cmparrGramInv[4];
	fncCreateVectV(valMu0,  cmpVectV0) ;
	TComp cmpNorm2V = TComp::scalProd(cmpVectV0, cmpVectV0, 4);
	TComp cmpVHS = TComp::scalProd(cmpVectV0, pcmpSZv, 4);
	TComp cmpCoef = cmpVHS/ cmpNorm2V;
	TComp cmparrProj[4], cmparrPer[4] ;
	MatrxMultScalar(cmpVectV0, 4, 1, cmpCoef,cmparrProj);
	MtrxMinusMatrx(pcmpSZv , cmparrProj,4, 1, cmparrPer);
	TComp cmpRez =  TComp::scalProd(cmparrPer, cmparrPer, 4);
	parrObj[ i * iQUANT_STEPS + j] =  cmpRez.m_Re;
	continue;
	}
	TComp cmparrT0[2],cmparrT1[4],cmparrT2[4],cmparrT3[4];
	MtrxMultMatrx(cmparrVH,2, 4, pcmpSZv, 1, cmparrT0);
	MtrxMultMatrx(cmparrTempInv,2, 2, cmparrT0, 1, cmparrA);

	MtrxMultMatrx(cmparrV,4, 2, cmparrA, 1, cmparrT1);
	MtrxMinusMatrx(pcmpSZv, cmparrT1,4, 1, cmparrT2);
	HermiteSoprMatr(cmparrT2, 4, 1, cmparrT3) ;
	TComp cmpF(0.,0.);
	MtrxMultMatrx(cmparrT3,1, 4, cmparrT2, 1, &cmpF);
	parrObj[ i * iQUANT_STEPS + j] =  cmpF.m_Re;

	// альтернативное решение
  /*	TComp cmpVectV0[4], cmpVectV1[4], cmparrX[2], cmparrB[2]
	, cmparrProjection[4],cmparrPerp[4], cmparrGram[4], cmparrGramInv[4];
	fncCreateVectV(valMu0,  cmpVectV0) ;
	fncCreateVectV(valMu1,  cmpVectV1) ;
	cmparrGram[0] = TComp::scalProd(cmpVectV0, cmpVectV0, 4);
	cmparrGram[1] = TComp::scalProd(cmpVectV1, cmpVectV0, 4);
	cmparrGram[2] = TComp::scalProd(cmpVectV0, cmpVectV1, 4);
	cmparrGram[3] = TComp::scalProd(cmpVectV1, cmpVectV1, 4);

	InverseCmpMtrx2__( cmparrGram, cmparrGramInv);

	cmparrB[0] =  TComp::scalProd(pcmpSZv, cmpVectV0, 4);
	cmparrB[1] =  TComp::scalProd(pcmpSZv, cmpVectV1, 4);

	MtrxMultMatrx(cmparrGramInv,2, 2, cmparrB,1, cmparrX);

	MtrxMultMatrx(cmparrV,4, 2, cmparrX,1, cmparrProjection);

	MtrxMinusMatrx(pcmpSZv, cmparrProjection,4, 1, cmparrPerp);

	TComp cmpObjTemp = TComp::scalProd(cmparrPerp, cmparrPerp, 4);  */

	int iiiiii=0;

  }
 // double  parrZ[100000]={0.};
  const int nrows = iQUANT_STEPS;
  const int ncols = iQUANT_STEPS;
  const double xllcorner = -VALDiapazon;
  const double yllcorner = -VALDiapazon;
  const double cellsize = VALStep;
  const double NODATA_value = DBL_MAX;
 // for (int i = 0; i < nrows/ 2; i++)
 // for (int j = 0; j < ncols; j++)
 // {
 //  parrZ [i * ncols + j] = 1.;
 // }
   TYrWrite::WriteMassiveInFltFile(L"E:\\СТАТЬЯ_Р_и_Э_2016_v2\\New\\rastr1.flt", parrObj,  nrows,
							  ncols, xllcorner , yllcorner,
							 cellsize, NODATA_value ) ;
  int inum = -1;
  *pvalMinObj = dblArrMin( parrObj,iQUANT_STEPS * iQUANT_STEPS, &inum) ;
  int i0 =  inum/iQUANT_STEPS;
  int j0 =  inum % iQUANT_STEPS;
  valMu0 = -VALDiapazon + ((double)i0) * VALStep;
  valMu1 = -VALDiapazon + ((double)j0) * VALStep;
  double temp0 = valMu0 /VALNu;
	if (temp0 < -0.99999999999) temp0 = -0.99999999999;
	if (temp0 > 0.99999999999) temp0 = 0.99999999999;
	*palfTrg = asin( temp0);

	double temp1 = valMu1 /VALNu;
	if (temp1 < -0.99999999999) temp1= -0.99999999999;
	if (temp1 > 0.99999999999) temp1 = 0.99999999999;
	*palfAnp  = asin( temp1);

  // find_ZTarg_and_ZAnt(pcmpSZv, pZTarg, pZAnt ,*palfTrg,  *palfAnp );

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
// solvNewtonMeth_Razn(valSigNoise,pcmpSZv, pZTarg, pZAnt
 //  	   , palfTrg, palfAnp, arrMtrxCorr ) ;

  return 0;

}

double  dblArrMin( double *arr,const int LENArr, int *pinum)
{
	double valreturn =  DBL_MAX;

	*pinum = -1;
	for (int i = 0; i < LENArr; i++)
	{
	  if (arr[i] < valreturn )
	  {
		valreturn =  arr[i];
		*pinum = i;

	  }
	}
	return valreturn;


}


void fncCreateMtrxV(const double valMu0, const double valMu1, TComp *cmparrV)
{
  for (int i =0; i < 4; i++)
  {
	cmparrV [ 2 * i] = exp_ (((double)i) * valMu0);
	cmparrV [ 2 * i + 1] = exp_ (((double)i) * valMu1);
  }
}

void fncCreateVectV(const double valMu,  TComp *cmpVectV)
{
  for (int i =0; i < 4; i++)
  {
	cmpVectV [  i] = exp_ (((double)i) * valMu);

  }
}





#pragma package(smart_init)
