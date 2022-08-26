//---------------------------------------------------------------------------
#pragma hdrstop
#include <math.h>
#include <wchar.h>
#include "ProcessExperDiagrs.h"
#include "UrPointXY.h"
#include "URPolyLine.h"
#include "DiagrSinX.h"
#include "YrWriteShapeFile.h"
#include "MeasStand.h"

extern const double TET0707;

//------------------------------------------------------------------------------
void TProcessExperDiagrs::extractDiagrArray(double *parrData, int iNumRows, int iNumCol, double *parrOut)
{ for(int i = 0; i < iNumRows; i++)
	{ parrOut[i] =  parrData[ i * 7 +  iNumCol];
	}
}
//------------------------------------------------------------------------------
TURPolyLine  TProcessExperDiagrs::extractLevelSubLine( TURPolyLine plnInp, double valLevel)
{
    TURPolyLine plnInp1 =  plnInp;
	int iarg = -1;
	int iNumPoints = 0;
	int iFirstPoint = -1;
	int iEndPoint = -1;

	for (int i =0; i < plnInp1.NumPoints; i++)
	{ if(plnInp1.Points[i].Y >= valLevel )
		{ if(iNumPoints ==0)
			{ iFirstPoint = i;
				iNumPoints++;
				continue;
			}
			if (iEndPoint == -1)  iNumPoints++;
		}
		if (plnInp1.Points[i].Y < valLevel )
		{ if (iFirstPoint == -1) continue;
			else  { iEndPoint = i-1; break;}
		}
	}
	TURPolyLine plnOut( &(plnInp1.Points[iFirstPoint]), iNumPoints);

	return plnOut;
}
//------------------------------------------------------------------------------
void  TProcessExperDiagrs::normDiagr( TURPolyLine &plnInp, int &iarg)
{ double valMax = -1000000.;

	iarg = -1;
	for(int i =0; i < plnInp.NumPoints; i++)
	{ if( plnInp.Points[i].Y > valMax)
		{ valMax =plnInp.Points[i].Y;
		 	iarg =i;
		}
	}

	double temp =  plnInp.Points[iarg].X ;
	for (int i =0; i < plnInp.NumPoints; i++)
	{ plnInp.Points[i].Y = plnInp.Points[i].Y / valMax ;
	  plnInp.Points[i].X =  plnInp.Points[i].X - temp ;
	 // plnInp.Points[i].X *= M_PI/ 3000.  ;
	}
}
//------------------------------------------------------------------------------
double  TProcessExperDiagrs::approximateDiagr( TURPolyLine &plnInp
	, double eps,double *pvalCoeff, double *pvalX)
{
	bool brez = false;
	double valCoeffTemp = TET0707/50. * 1000.;
	double valXTemp =  0.;
	for (int i =0 ; i < 100; i++)
	{
		double valCoeffTemp1 = findOptCoeff ( plnInp,  valXTemp ) ;
		double valXTemp1 = findOptX ( plnInp,   valCoeffTemp1 ) ;
		if ((fabs(valCoeffTemp1 - valCoeffTemp) < eps) && (fabs(valXTemp1 - valXTemp) < eps) )
		{
			*pvalCoeff = valCoeffTemp1;
			*pvalX =  valXTemp1;
			brez = true;
			break;
		}
		valCoeffTemp = valCoeffTemp1;
		valXTemp =  valXTemp1;
	}
	return sqrt(calcNeviazka(plnInp, *pvalCoeff, *pvalX)/plnInp.NumPoints);
}
//------------------------------------------------------------------------------
double TProcessExperDiagrs::findOptCoeff ( TURPolyLine &plnInp, double valXTemp )
{ // перебор по диаграммам с tet07 от 1 мрад до 51 мрад с шагом 0,05мрад
	double min = 10000000000.;
	double optCoeff= -1.;
	double step = 0.05;
	int iNc = 50. / step;
	for (int i =0; i < iNc; i++)
	{ double valTet07 = (1. + ((double)i) * step)/ 1000.;
		double valCoeff = TET0707/valTet07;
		double valNeviazka = calcNeviazka(plnInp, valCoeff, valXTemp);
		if (valNeviazka < min)
		{
			min = valNeviazka;
			optCoeff = valCoeff;
		}
	}
	return optCoeff;
}
//------------------------------------------------------------------------------
double TProcessExperDiagrs::findOptX ( TURPolyLine &plnInp, double valCoeff )
{ // перебор по диаграммам со сдвигом от -25 мрад до 25 мрад  с шагом 0,05 мрад
	double min  = 10000000000.;
	double optX0= -10000000.;
	double step = 0.05;
	int iNc = 50. / step;
	for (int i=0;i<iNc;i++)
	{
		double valX = (-25. + ((double)i) * step)/ 1000.;
		double valNeviazka = calcNeviazka(plnInp, valCoeff, valX);
		if (valNeviazka < min)
		{
		    min 	= valNeviazka;
			optX0= valX;
		}
	}
	return optX0;
}
//------------------------------------------------------------------------------
double TProcessExperDiagrs::calcNeviazka(TURPolyLine &plnInp, double valCoeff, double valX)
{ double sum = 0.;
	double a;
	for (int i = 0; i < plnInp.NumPoints; i++)
	{ a=  plnInp.Points[i].Y - fncDiagrSinx_div_x(valCoeff * (plnInp.Points[i].X - valX));
		a=  a/fncDiagrSinx_div_x(valCoeff * (plnInp.Points[i].X - valX));
		sum += a*a;
		//sum +=(plnInp.Points[i].Y - fncDiagrSinx_div_x(valCoeff * (plnInp.Points[i].X - valX)))*
		//(plnInp.Points[i].Y - fncDiagrSinx_div_x(valCoeff * (plnInp.Points[i].X - valX)));
	}
	return sum;
}
//------------------------------------------------------------------------------
 // построние графика парциальной диаграммы смещенной на  valSmeshenieGoriz и valSMeshenieVert
void TProcessExperDiagrs::createGraphDiagrAndApproxDiagr(wchar_t *Fold,TURPolyLine &plnInp,  double valCoeff, double valX)
{   TURPolyLine plnSinx = plnInp ;
	TURPolyLine plnTemp = plnInp ;

	for (int i=0 ; i < plnInp.NumPoints; i++)
	{   plnTemp.Points[i].X =  valCoeff * (plnInp.Points[i].X - valX);
		plnSinx.Points[i].X =  valCoeff * (plnInp.Points[i].X - valX);
		plnSinx.Points[i].Y = fncDiagrSinx_div_x(plnSinx.Points[i].X);
	}
	wchar_t wchFileName1[300] ={0};
	wcscpy(  wchFileName1,  Fold);
	wcscat(wchFileName1, L"\\DiagrOrig.shp");
	plnSinx.WriteSetSHPFiles(wchFileName1, &plnTemp, 1);


	wchar_t wchFileName2[300] ={0};
	wcscpy( wchFileName2,  Fold);
	wcscat(wchFileName2, L"\\SinXDidX.shp");
	plnSinx.WriteSetSHPFiles(wchFileName2, &plnSinx, 1);

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  Fold);
	wcscat(wchAxesFileName, L"\\AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-200., 200.
	 ,-13.,13.5,  0.05) ;
}
//------------------------------------------------------------------------------

int TProcessExperDiagrs::extractFreqNum_FromFileName(wchar_t *FileName)
{
  wchar_t *pstr = wcsrchr(FileName, L'\\') ;
  for (int i =0; i < 100; i++)
  {
	if ((pstr[0] >= 48) && (pstr[0] <=  57) )
	{
	  break;
	}
	pstr++;
  }
  wchar_t arrstr[2] ={0};
  arrstr[0] =  pstr[0];
  int nFreq= _wtoi(arrstr);

	return nFreq;
}


// извлечение массива измерений тройки диаграмм с номером  iNumTriple
// из буфера parrBuff, полученного из текстового файла Каурова
// NUmRows - к-во строк буфкера
// 1 + 2 * NUmDiagr - к-во столбцов буфера
// NUmDiagr  - к-во измереннных диаграмм
// parrMeas - массив измерений
// формат столбцов буфера:
// parrBuff[0] - угол наклона центральной оси антенны относительно горизонта
// центральная ось проходит через середину угла, образованного осями 3 и 4 диаграмм
// parrBuff[1] по parrBuff[NUmDiagr]  - амплитуды диаграмм
// parrBuff[NUmDiagr + 1] по parrBuff[2 *NUmDiagr]  - фазы  диаграмм  в градусах
void TProcessExperDiagrs::extractMeasuresArrFromKaurovStandFile(double *parrBuff, const int NUmRows
   , int NUmDiagr, int iNumTriple, TMeasStand *parrMeas)
{
	for (int i = 0; i < NUmRows; i++)
	{
	double *parr = &parrBuff[ i * ( 2 * NUmDiagr + 1)];

	 parrMeas[i].mRadAxeAnt = parr[0] * M_PI / 3. / 1000.;
	 for (int j = 0; j < 3; j++)
	 {

	   double alf = parr[  NUmDiagr + iNumTriple + j] * M_PI / 180.;
	   parrMeas[i].mpcmparrMeas[j].m_Re = parr[ 1 +  iNumTriple + j] * cos (alf);
	   parrMeas[i].mpcmparrMeas[j].m_Im = parr[ 1 +  iNumTriple + j] * sin (alf);


	 }
	}
}


  // ПЕРЕГРУЖЕННАЯ!!!
// извлечение массива измерений ансамбля  диаграмм с номером  iNumAns
// из буфера parrBuff, полученного из текстового файла Каурова
// NUmRows - к-во строк буфкера
// 1 + 2 * NUmDiagr - к-во столбцов буфера
// NUmDiagr  - к-во измереннных диаграмм
// parrMeas - массив измерений
// lenAns - длина ансамбля
// формат столбцов буфера:
// parrBuff[0] - угол наклона центральной оси антенны относительно горизонта
// центральная ось проходит через середину угла, образованного осями 3 и 4 диаграмм
// parrBuff[1] по parrBuff[NUmDiagr]  - амплитуды диаграмм
// parrBuff[NUmDiagr + 1] по parrBuff[2 *NUmDiagr]  - фазы  диаграмм  в градусах
void TProcessExperDiagrs::extractMeasuresArrFromKaurovStandFile(double *parrBuff, const int NUmRows
   , int NUmDiagr, int iNumAns, int lenAns,  TMeasStand *parrMeas)
{
	for (int i = 0; i < NUmRows; i++)
	{
	double *parr = &parrBuff[ i * ( 2 * NUmDiagr + 1)];

	 parrMeas[i].mRadAxeAnt = parr[0] * M_PI / 3. / 1000.;
	 for (int j = 0; j < lenAns; j++)
	 {

	   double alf = parr[  NUmDiagr + iNumAns + j] * M_PI / 180.;
	   parrMeas[i].mpcmparrMeas[j].m_Re = parr[ 1 +  iNumAns + j] * cos (alf);
	   parrMeas[i].mpcmparrMeas[j].m_Im = parr[ 1 +  iNumAns + j] * sin (alf);


	 }
	}
}

void TProcessExperDiagrs::drawAplitudeGraphs(wchar_t *wchFileName,double *parrBuff, int numDiagrBegin, int countDiagr
   , const int NUmBuffRows, const int NUmBuffCols, TURPointXY pntSdvig, double scalex, double  scaley)
{
  TURPolyLine *ppln = new TURPolyLine[countDiagr];//( const int iNumParts, const int iNumPoints) ;
   for (int i =0; i < countDiagr; i++)
   {
	 ppln[i] =  TURPolyLine(1, NUmBuffRows);
   }
   for (int i = 0; i < NUmBuffRows; i++)
   {
	  double *pCur = &parrBuff[ NUmBuffCols *i];
	  for (int j =0; j < countDiagr; j++)
	  {
		ppln[j].Points[i].X =  pCur[0] * M_PI / 3. / 1000. * scalex  + pntSdvig.X;
		ppln[j].Points[i].Y =  pCur[ 1 + numDiagrBegin + j] * scaley + pntSdvig.Y;
	  }
   }

  TURPolyLine::WriteSetSHPFiles(wchFileName, ppln, countDiagr);
  delete [] ppln;
}

void   TProcessExperDiagrs::drawApproximatedAplitudeGraphs(wchar_t *wchFileName1, double *arrDiagrSdvigX
	,double *arrDiagrMultCoeffY, double *arrDiagrMultCoeffX )
{
  int countDiagr = 6;
 TURPolyLine *ppln = new TURPolyLine[countDiagr];//( const int iNumParts, const int iNumPoints) ;
   for (int i =0; i < countDiagr; i++)
   {
	 TURPointXY pntSdvig0( arrDiagrSdvigX[i], -50.);
	 ppln[i] =  TURPolyLine::create_SinX_Div_X_Line( pntSdvig0, arrDiagrMultCoeffX[i], arrDiagrMultCoeffY[i], 2);
   }


  TURPolyLine::WriteSetSHPFiles(wchFileName1, ppln, countDiagr);


  delete [] ppln;
}

#pragma package(smart_init)
