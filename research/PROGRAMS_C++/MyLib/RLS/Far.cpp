//---------------------------------------------------------------------------


#pragma hdrstop
#include "Far.h"

 #include "Faceta.h"


#include <vcl.h>


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
#include "AM_2D.h"
#include "SincDgr.h"
#include "Diagrams.h"
#include "DiagrSinX.h"

extern const double TET0707;
//#include "DiagrSet.h"

TFar::~TFar()
{
	if(mparrDisp)
	{
		free( mparrDisp) ;
		mparrDisp = NULL ;
		}

	if(mparrAmplFactDisp) free( mparrAmplFactDisp) ;
	mparrAmplFactDisp = NULL ;

}
TFar::TFar()
{
 // количество модулей (фасет)
	m_N = 0;
 // расстояние между модулями
	 m_D = 0;
 // длина волны
	 mLambda = 0;
 // массив фасет
  mFaceta = TFaceta();
  //
	mparrDisp = NULL;
  mparrAmplFactDisp = NULL;

}
// Конструктор копирования
TFar::TFar (const TFar &R2)
 {
	// количество модулей (фасет)
	m_N = R2.m_N;
	// расстояние между модулями
	m_D = R2.m_D;
	// длина волны
	mLambda = R2.mLambda;
	mFaceta = R2.mFaceta;



	if(R2.mparrDisp != NULL)
	{
	mparrDisp = (double*)malloc(m_N * sizeof(double));

	if(mparrDisp == NULL)
	{
	ShowMessage(L"Not memory for Parts") ;
	Abort() ;
	}

	memcpy( mparrDisp,R2.mparrDisp, R2.m_N  * sizeof(double));
	}
	///


	if(R2.mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));

	if(mparrAmplFactDisp == NULL)
	{
	ShowMessage(L"Not memory for Parts") ;
	Abort() ;
	}

	memcpy( mparrAmplFactDisp,R2.mparrAmplFactDisp, R2.m_N  * sizeof(double));
	}

 }

// парам констр  1
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta)
 {
	 m_N = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;
	 ///
	if (mparrDisp != NULL)
	{
	mparrDisp = NULL;
	}
	mparrDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}

	if (mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = NULL;
	}
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrAmplFactDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}
	for (int i = 0; i < m_N; i++)
	{
	 mparrDisp[i] =  Faceta.calcNoiseSKZ() * Faceta.calcNoiseSKZ();
	 mparrAmplFactDisp[i] = Faceta.calcAmplFactSKZ() * Faceta.calcAmplFactSKZ() ;
	}


 }

 // парам констр  2
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double *parrDisp)
 {
	 m_N = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;
	 if (mparrDisp != NULL)
	 {
	 mparrDisp = NULL;
	 }
	 if(parrDisp != NULL)
	 {
	   mparrDisp = (double*)malloc(m_N * sizeof(double));
	   if (mparrDisp == NULL)
	   {
		   ShowMessage(L"Error 000");
	   }
	   memcpy(mparrDisp, parrDisp,m_N * sizeof(double));
	 }
	 ///
	if (mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = NULL;
	}
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrAmplFactDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}
	for (int i = 0; i < m_N; i++)
	{
	 mparrAmplFactDisp[i] = Faceta.calcAmplFactSKZ() * Faceta.calcAmplFactSKZ() ;
	}

 }


 // парам констр  2
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double *parrDisp, double *parrAmplFactDisp)
 {
	 m_N = N ;
	 m_D= D;
	 mLambda = Lambda ;
	 mFaceta = Faceta;
	 if (mparrDisp != NULL)
	 {
	 mparrDisp = NULL;
	 }
	 if(parrDisp != NULL)
	 {
	   mparrDisp = (double*)malloc(m_N * sizeof(double));
	   if (mparrDisp == NULL)
	   {
		   ShowMessage(L"Error 000");
	   }
	   memcpy(mparrDisp, parrDisp,m_N * sizeof(double));
	 }
	 ///
	if (mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = NULL;
	}
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrAmplFactDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}
	 memcpy(mparrAmplFactDisp, parrAmplFactDisp,m_N * sizeof(double));
	 double sum0 =0., sum1 = 0.;
	 for (int i = 0; i < N; i++)
	 {
	   sum0 +=  parrDisp[i];
	   sum1 +=  parrAmplFactDisp[i];
	 }

	 mFaceta.mSigEmitAmplFact = sqrt(sum1 / ((double)m_N) *  ((double)mFaceta.m_n));
	 mFaceta.mSigEmitNoise =  sqrt(sum0)/ sqrt((double)(m_N *  mFaceta.m_n));

 }

 // парам констр  3
 __fastcall TFar::TFar(const int N,const double D,const double Lambda
   ,TFaceta Faceta, double valFacetaDisp)
 {
	m_N = N ;
	m_D= D;
	mLambda = Lambda ;
	mFaceta = Faceta;
	mFaceta.mSigEmitNoise =  sqrt(valFacetaDisp/ mFaceta.m_n);
	mFaceta.mSigEmitAmplFact = 0.;
	if (mparrDisp != NULL)
	{
	mparrDisp = NULL;
	}
	mparrDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}

	if (mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = NULL;
	}
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrAmplFactDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}
	for (int i = 0; i < m_N; i++)
	{
	 mparrDisp[i] =  valFacetaDisp;
	 mparrAmplFactDisp[i] = 0. ;
	}


 }

 // парам констр   4
 __fastcall TFar::TFar(const int N)
 {

 // количество модулей (фасет)
  m_N = N;
 // расстояние между модулями
   m_D = 32.8;
 // длина волны
   mLambda = 3.15;;
 // массив фасет
  mFaceta = TFaceta();
  	 ///
	if (mparrDisp != NULL)
	{
	mparrDisp = NULL;
	}
	mparrDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}

	if (mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = NULL;
	}
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
	if (mparrAmplFactDisp == NULL)
	{
	   ShowMessage(L"Error 000");
	}
	for (int i = 0; i < m_N; i++)
	{
	 mparrDisp[i] =  mFaceta.calcNoiseSKZ() * mFaceta.calcNoiseSKZ();
	 mparrAmplFactDisp[i] = mFaceta.calcAmplFactSKZ() * mFaceta.calcAmplFactSKZ() ;
	}


 }

 // парам конструктор   создание ФАР с заданным числом строк  NumRows
 // из ФАР Far
 //  строки куонструируесмой ФАР являются суммами соответствующих строк  Far
 //
 //
 //
  __fastcall TFar::TFar(const TFar Far, const int NumRows)
{
  if (Far.m_N % NumRows != 0)
  {
	return ;
  }
 //
 mLambda =  Far.mLambda;
 int qRow = Far.m_N / NumRows;  // сколько исходных  строк в строке
 double *parrDisp = new double [NumRows];
 double *parrAmpFactDisp = new double [NumRows];
 memset(parrDisp, 0, NumRows * sizeof(double));
 memset(parrAmpFactDisp, 0, NumRows * sizeof(double));
 for (int i = 0; i < NumRows; i++)
 {

   for (int j =0; j < qRow; j++)
   {

	parrDisp[i] += Far.mparrDisp[ i * qRow + j];
	parrAmpFactDisp [i] += Far.mparrAmplFactDisp[i * qRow + j];
   }
   parrAmpFactDisp [i] = parrAmpFactDisp [i]/((double)qRow);  // ????????????
 }
 ///

 double valD = qRow * Far.m_D;


 TFaceta FacetaUn(Far.mFaceta.m_n * qRow  , Far.mFaceta.m_d, Far.mFaceta.mLambda, parrDisp[0], parrAmpFactDisp [0]);

 TFar FarUn(NumRows,valD , mLambda ,FacetaUn,parrDisp, parrAmpFactDisp );
 delete []parrDisp ;
 delete []parrAmpFactDisp;
 *this = FarUn;

}
 // оператор присваивания
  TFar &TFar::operator=(const TFar  &R2)
{
// количество модулей (фасет)
	m_N = R2.m_N;
	// расстояние между модулями
	m_D = R2.m_D;
	// длина волны
	mLambda = R2.mLambda;
	mFaceta = R2.mFaceta;

	mparrDisp = NULL;

	if(R2.mparrDisp != NULL)
	{
	mparrDisp = (double*)malloc(m_N * sizeof(double));

	if(mparrDisp == NULL)
	{
	ShowMessage(L"Not memory for Parts") ;
	Abort() ;
	}

	memcpy( mparrDisp,R2.mparrDisp, R2.m_N  * sizeof(double));
	}
///
	mparrAmplFactDisp = NULL;

	if(R2.mparrAmplFactDisp != NULL)
	{
	mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));

	if(mparrAmplFactDisp == NULL)
	{
	ShowMessage(L"Not memory for Parts") ;
	Abort() ;
	}

	memcpy( mparrAmplFactDisp,R2.mparrAmplFactDisp, R2.m_N  * sizeof(double));
	}

  return *this ;
}



// создание одномерной ФАР из 2-х мерной
// BVertical = true если по вертикали (углу места)
// BVertical = false если по горизонтали(курсовому углу)
 __fastcall TFar::TFar(TFar_2D Far_2D, const bool BVertical)
 {

	 if (BVertical)
	 {
		// количество модулей (фасет)
		m_N = Far_2D.mNumAMRows;
		// расстояние между модулями
		m_D = ((double)(Far_2D.mpAm2D[0].mNumEmitRows)) * Far_2D.mpAm2D[0].mdRow;
		// длина волны
		mLambda = Far_2D.mLambda;
		// фасета
		// поиск персвого рабочего АМ
		int iWorkingAM = -1;
		for (int i = 0; i < Far_2D.mNumAMCols * Far_2D.mNumAMRows; i++)
		{
		if (Far_2D.mpbarrAM [i])
		{
		iWorkingAM = i;
		break;
		}
		}
		const double SigEmitNoise = sqrt((double)(Far_2D.mpAm2D[iWorkingAM ].mNumEmitCols * Far_2D.mNumAMCols )  )* Far_2D.mpAm2D[iWorkingAM ].mSigEmitNoise;
		const double SigEmitAmplFact =  Far_2D.mpAm2D[iWorkingAM ].mSigEmitAmplFact/sqrt((double)(Far_2D.mpAm2D[iWorkingAM ].mNumEmitCols * Far_2D.mNumAMCols));
		mFaceta = TFaceta(Far_2D.mpAm2D[0].mNumEmitRows,Far_2D.mpAm2D[0].mdRow,mLambda
		, SigEmitNoise , SigEmitAmplFact  );
		///
		if (mparrDisp != NULL)
		{
		mparrDisp = NULL;
		}
		mparrDisp = (double*)malloc(m_N * sizeof(double));
		if (mparrDisp == NULL)
		{
		ShowMessage(L"Error 000");
		}
		memset(mparrDisp, 0, m_N * sizeof(double));
		if (mparrAmplFactDisp != NULL)
		{
		mparrAmplFactDisp = NULL;
		}
		mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
		if (mparrAmplFactDisp == NULL)
		{
		ShowMessage(L"Error 000");
		}
		memset(mparrAmplFactDisp, 0, m_N * sizeof(double));
		for (int i = 0; i < m_N; i++)
		{
		for (int j = 0; j < Far_2D.mNumAMCols; j++)
		{
		if (Far_2D.mpbarrAM[i * (Far_2D.mNumAMCols) + j])
		{
		mparrDisp[i] += Far_2D.mpAm2D[i * (Far_2D.mNumAMCols)+ j].calcSumDisp();
		mparrAmplFactDisp[i] += Far_2D.mpAm2D[i *(Far_2D.mNumAMCols) + j].calcSumAmplFactDisp();
		}
		}

	}


	 }
	 else
	 {
// количество модулей (фасет)
		m_N = Far_2D.mNumAMCols;
		// расстояние между модулями
		m_D = ((double)Far_2D.mpAm2D[0].mNumEmitCols) * Far_2D.mpAm2D[0].mdCol;
		// длина волны
		mLambda = Far_2D.mLambda;
		// фасета
		// поиск персвого рабочего АМ
		int iWorkingAM = -1;
		for (int i = 0; i < Far_2D.mNumAMCols * Far_2D.mNumAMCols; i++)
		{
		if (Far_2D.mpbarrAM [i])
		{
		iWorkingAM = i;
		break;
		}
		}
		const double SigEmitNoise = sqrt((double)(Far_2D.mpAm2D[iWorkingAM ].mNumEmitRows * Far_2D.mNumAMRows )  )* Far_2D.mpAm2D[iWorkingAM ].mSigEmitNoise;
		const double SigEmitAmplFact = sqrt((double)(Far_2D.mpAm2D[iWorkingAM ].mNumEmitRows * Far_2D.mNumAMRows )  )* Far_2D.mpAm2D[iWorkingAM ].mSigEmitAmplFact;
		mFaceta = TFaceta(Far_2D.mpAm2D[0].mNumEmitCols,Far_2D.mpAm2D[0].mdCol,mLambda
		, SigEmitNoise , SigEmitAmplFact  );
		///
		if (mparrDisp != NULL)
		{
		mparrDisp = NULL;
		}
		mparrDisp = (double*)malloc(m_N * sizeof(double));
		if (mparrDisp == NULL)
		{
		ShowMessage(L"Error 000");
		}
		memset(mparrDisp, 0, m_N * sizeof(double));
		if (mparrAmplFactDisp != NULL)
		{
		mparrAmplFactDisp = NULL;
		}
		mparrAmplFactDisp = (double*)malloc(m_N * sizeof(double));
		if (mparrAmplFactDisp == NULL)
		{
		ShowMessage(L"Error 000");
		}
		memset(mparrAmplFactDisp, 0, m_N * sizeof(double));
		for (int i = 0; i < m_N; i++)
		{
			for (int j = 0; j < Far_2D.mNumAMRows; j++)
			{
			if (Far_2D.mpbarrAM[i  +  j * Far_2D.mNumAMCols ])
			{
			mparrDisp[i] += Far_2D.mpAm2D[i  +  j * Far_2D.mNumAMCols ].calcSumDisp();
			mparrAmplFactDisp[i] += Far_2D.mpAm2D[i  +  j * Far_2D.mNumAMCols ].calcSumAmplFactDisp();
			}
			}
	   }

 }
}

// СКЗ суммарного шума в диаграмме ФАР
double TFar::calcNoiseDisp()
{
 double valrez = 0.;
 for (int i = 0; i <  m_N; i++)
 {
   valrez += mparrDisp[i];
 }
return valrez;
}

// СКЗ суммарного шума КУ ( 1 + delta) ФАР
double TFar::calcAmplFactDisp()
{

 double valrez = 0.;
 for (int i = 0; i <  m_N; i++)
 {
   valrez += mparrAmplFactDisp[i];
 }
 return valrez / sqrt((double)m_N);
}


// пеленгационная функция СРМ
// INPUT:
// pcmpSZv - массив замеров по строкам, комплексный
//  OUTPUT -  массив  частных производных пеленгац ф-ции по строковым диаграммам
// pcarr_dPel_po_dSi
 double TFar::fncPeleng(TComp *pcmpSZv, TComp *pcarr_dPel_po_dSi)
 {

  TComp cmpUp(0.,0.);//
  TComp cmpDown (0.,0.);//
  for (int i = 0; i < m_N/2; i++)
  {
	cmpUp += pcmpSZv[i];
	cmpDown += pcmpSZv[m_N/2 + i];
  }
  TComp cmpRazn = cmpUp -   cmpDown;
  TComp cmpSum  = cmpUp +   cmpDown;
  TComp cmpPel =  cmpRazn/cmpSum;

   if(pcarr_dPel_po_dSi)
   {
	TComp cmp2(2.,0.);
	TComp cmpMin2(-2.,0.);
	for (int i = 0; i < m_N/2; i++)
	{
	TComp cmpt0 = cmp2 *cmpDown/ TComp(1000.,0.);
	TComp cmpt1 = (cmpSum*cmpSum)/ TComp(1000.,0.);
	pcarr_dPel_po_dSi[           i] = cmpt0 / cmpt1 ;

	TComp cmpt2 = cmpMin2 *cmpUp/ TComp(1000.,0.);
	pcarr_dPel_po_dSi[m_N/2 +i] = cmpt2 / cmpt1 ;
	}
  }

  return cmpPel.m_Im;
 }

double TFar::calcDispPelengFnc(TComp *pcmpSZv)
{
  TComp *pcarr_dPel_po_dSi = (TComp*)malloc(m_N * sizeof(TComp));
  fncPeleng(pcmpSZv, pcarr_dPel_po_dSi);
  double rez = 0.;
  for (int i = 0; i < m_N; i++)
  {
	rez +=  pcarr_dPel_po_dSi[i].modul() * pcarr_dPel_po_dSi[i].modul() *  mparrDisp[i]/2. ;
  }

  free(pcarr_dPel_po_dSi);
  return rez;

}

// вычисление оценки угла РСМ и расчет дисперсии ошибки для пеленг
// функции tg(y/2)
//  INPUT:
 // pcmpSZv - массив замеров по строкам, комплексный
 // OUTPUT:
 // pDisp - дисперсия ошибки оценивания
double TFar::fncEstimateRsmTetta(TComp *pcmpSZv,  double *pDisp )
{
	double valPelFnc = fncPeleng(pcmpSZv, NULL);
	double valDispPelFnc = calcDispPelengFnc(pcmpSZv);
	double valCoef =  mFaceta.mLambda/( M_PI * m_D * ((double)m_N)/ 2.);
	double valTet = -atan(valPelFnc)* valCoef;
	double valDeriv = 1./ (1. + valPelFnc * valPelFnc );
	*pDisp =  valDispPelFnc * valDeriv * valDeriv * valCoef* valCoef;
	return valTet;
}


// INPUT:
// pcmpSZv - скорректированный массив измерений строковых диаграмм
// OUTPUT:
// valXi2 - хи квадрат
// valEstRSM - измерение угла в мдрд
// valRSMDisp - дисперсия ошибки измерения угла
//

void TFar::fncMeasureRSMProcessing(TComp *pcmpSZv, double *valXi2, double *valEstRSM, double *valRSMDisp  )
{
	*valEstRSM = fncEstimateRsmTetta (pcmpSZv,  valRSMDisp );

	double valMuZv =  sin(*valEstRSM) * m_D * 2. * M_PI/ mLambda;
	*valXi2 = fncXi2_5P10_UM_GeneralCase( pcmpSZv, valMuZv) ;


 }

 // вычисление функции хи квадрат для изделия 5П10 для случая работы в штатном режиме
 // - то есть, когда все модули в рабочем состоянии
//по углу места
// INPUT:
// valNoiseSKZ - шум в суммарной диаграмме
// pcmpSZv - замеры строковых диаграмм, скомпенсированные по модулю
//valMuZv - величина Mu = D* sin(tet)*2* M_PI / lambda
//
double  TFar::fncXi2_5P10_UM(double valNoiseSKZ, TComp *pcmpSZv, double valMuZv)
{
	TComp cmparrQZv[4];
	for (int i=0; i < 4; i++)
	{
	  TComp cmpRotate = exp_(-((double)i) * valMuZv);
	  cmparrQZv[i] = pcmpSZv[i] * cmpRotate;
	}
	TComp cmp0_75(3./4.,0);
	TComp cmpB =  cmparrQZv[0] +  cmparrQZv[3];
	cmpB *= cmp0_75;
	cmpB  += cmparrQZv[1];
	cmpB  += cmparrQZv[2];

	double temp0 = 21. *(cmparrQZv[0].modul()*cmparrQZv[0].modul() + cmparrQZv[3].modul()*cmparrQZv[3].modul())/8. ;
	double temp1 = 7. *(cmparrQZv[1].modul()*cmparrQZv[1].modul() + cmparrQZv[2].modul()*cmparrQZv[2].modul())/2. ;

	double temp2 = temp0 + temp1 -   cmpB.modul()*cmpB.modul();
	double temp3 = temp2/valNoiseSKZ/valNoiseSKZ ;

	// альтернативный рассчет
	double sig0Sq = valNoiseSKZ*valNoiseSKZ * 8./ 21.;
	double sig1Sq = valNoiseSKZ*valNoiseSKZ * 2./ 7.;
	double sig0Sqmin1 = 1./ sig0Sq;
	TComp cmpsig0Sqmin1(sig0Sqmin1,0.);
	double sig1Sqmin1 = 1./ sig1Sq;
	TComp cmpsig1Sqmin1(sig1Sqmin1,0.);

	TComp cmpsigSqSum ( 1./(sig0Sqmin1 +  sig1Sqmin1),0.);
	TComp cmpT1 = cmparrQZv[0] +  cmparrQZv[3];
	TComp cmpT2 = cmpT1 *cmpsig0Sqmin1;

	TComp cmpT3 =  cmparrQZv[1] +  cmparrQZv[2];
	TComp cmpT4 =  cmpT3 *cmpsig1Sqmin1;

	TComp cmpT5 = cmpT2 + cmpT4;
	return temp3 ;
}

 /*
// вычисление функции хи квадрат для изделия 5П10 для случая рабты в штатном режиме
 // - то есть, когда все модули в рабочем состоянии
//по углу места
// INPUT:
// valNoiseSKZ - шум в суммарной диаграмме
// pcmpSZv - замеры строковых диаграмм, скомпенсированные по модулю
//valMuZv - величина Mu = D* sin(tet)*2* M_PI / lambda
//
double  TFar::fncXi2_5P10_UM_GeneralCase(TComp *pcmpSZv, double valMuZv)
{
	TComp cmparrQZv[4];
	TComp cmpQZvSum(0.,0.), cmpSumSigMin2(0.,0.);
	for (int i=0; i < 4; i++)
	{
	  TComp cmpRotate = exp_(-((double)i) * valMuZv);
	  cmparrQZv[i] = pcmpSZv[i] * cmpRotate;
	  TComp cmpDispMinus2(1./mparrDisp[i], 0.);
	  cmpSumSigMin2 += cmpDispMinus2;
	  cmpQZvSum +=  cmparrQZv[i]* cmpDispMinus2;
	  }

	TComp cmpA0 =  cmpQZvSum/  cmpSumSigMin2;
	double valF0 = 0.;
	for (int i =0; i < 4; i++)
	{
	 TComp cmpTemp = cmparrQZv[i] -cmpA0;
	 valF0 +=  cmpTemp.modul() * cmpTemp.modul()/mparrDisp[i];
	}
	valF0 = valF0;
   return valF0;
}
*/

// вычисление функции хи квадрат для изделия 5П10 для случая рабты в штатном режиме
 // - то есть, когда все модули в рабочем состоянии
//по углу места
// INPUT:
// valNoiseSKZ - шум в суммарной диаграмме
// pcmpSZv - замеры строковых диаграмм, скомпенсированные по модулю
//valMuZv - величина Mu = D* sin(tet)*2* M_PI / lambda
//
double  TFar::fncXi2_5P10_UM_GeneralCase(TComp *pcmpSZv, double valMuZv)
{
	TComp *cmparrQZv = new TComp[m_N];
	TComp cmpQZvSum(0.,0.), cmpSumSigMin2(0.,0.);
	for (int i=0; i < m_N; i++)
	{
	  TComp cmpRotate = exp_(-((double)i) * valMuZv);
	  cmparrQZv[i] = pcmpSZv[i] * cmpRotate;
	  TComp cmpDispMinus2(1./mparrDisp[i], 0.);
	  cmpSumSigMin2 += cmpDispMinus2;
	  cmpQZvSum +=  cmparrQZv[i]* cmpDispMinus2;
	  }

	TComp cmpA0 =  cmpQZvSum/  cmpSumSigMin2;
	double valF0 = 0.;
	for (int i =0; i < m_N; i++)
	{
	 TComp cmpTemp = cmparrQZv[i] -cmpA0;
	 valF0 +=  cmpTemp.modul() * cmpTemp.modul()/mparrDisp[i];
	}
	valF0 = valF0;
	delete []cmparrQZv ;
   return valF0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////// ФУНКЦИИ МЕТОДА СТРОКОВЫХ ДИАГРАММ  /////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // OUTPUT:
 // * z1- цель
 // *z2 -  антипод
 // arrMtrxCorr_fi - корреляционная матирица
 int   TFar::solvQuadrEqMeth(TComp *pcmpSZv, TComp * z1, TComp *z2 , double *arrMtrxCorr_fi0 )
{
  TComp cmpa(1., 0.);

  TComp cmpb = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0] * pcmpSZv[3] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  TComp cmpc = ( pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2] * pcmpSZv[2] )/ ( pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1] * pcmpSZv[1] ) ;
  SolvEq2( cmpa, cmpb,  cmpc, *z1, *z2);

  // вычисление корреляц матрицы ошибок аргументов z1, z2

  memset( arrMtrxCorr_fi0, 0, 4 * sizeof(double));
  TComp y1(0.,0.), y2(0.,0.), cmparrHInv[16];
  y1 = ( pcmpSZv[0] * (*z2) - pcmpSZv[1])/ ((*z2) - (*z1));
  y2 = ( pcmpSZv[1] - pcmpSZv[0] * (*z1))/ ((*z2) - (*z1));


  TComp cmp3(3.,0.), cmpMin1(-1.,0.), cmpMin3(-3.,0.), cmp2(2.,0.)
	,cmpMin2(-2.,0.), cmp6(6.,0.), cmpMin6(-6.,0.), cmpSq1(1. + sqrt(3.),0.)
	, cmpSq2(1. - sqrt(3.),0.);
  cmparrHInv[0]  = y1 * y2 *(*z2)*(*z2) * ((*z2) - (cmp3 * (*z1)));
  cmparrHInv[1]  = cmp6 *y1 * y2 *(*z1)*(*z2);
  cmparrHInv[2]  = cmpMin3 * y1 * y2 *((*z2) + (*z1));
  cmparrHInv[3]  = cmp2 * y1 * y2 ;

  cmparrHInv[4]  = cmpMin1 * (*z1)* (*z2) * (*z2) * y2  * ((*z2) -  (*z1));
  cmparrHInv[5]  = (*z2) * y2 *  ((*z2) - (*z1))*  ((*z2) + cmp2 *(*z1));
  cmparrHInv[6]  = cmpMin1 * y2  *((*z2) - (*z1))*  ((*z1) + cmp2 *(*z2));
  cmparrHInv[7]  =  y2 * ((*z2) -  (*z1));

  cmparrHInv[8]   = cmpMin1 * y1 * y2 *(*z1)*(*z1) * ((*z1) - cmp3 * (*z2));
  cmparrHInv[9]   = cmpMin6 *y1 * y2 *(*z1)*(*z2);
  cmparrHInv[10]  = cmp3 * y1 * y2 *((*z2) + (*z1));
  cmparrHInv[11]  = cmpMin2 * y1 * y2 ;

  cmparrHInv[12]  = cmpMin1 * (*z1)* (*z1) * (*z2) * y1  * ((*z2) -  (*z1));
  cmparrHInv[13]  = (*z1) * y1 *  ((*z2) - (*z1))*  ((*z1) + cmp2 *(*z2));
  cmparrHInv[14]  = cmpMin1 * y1 * ((*z2) - (*z1))*  ((*z2) + cmp2 *(*z1));
  cmparrHInv[15]  =  y1 * ((*z2) -  (*z1));

  TComp cmpTemp = y1 * y2* ((*z2) -  (*z1)) * ((*z2) -  (*z1)) *((*z2) -  (*z1));
  for (int i = 0; i <16; i++)
  {
   cmparrHInv[i] = cmparrHInv[i] / cmpTemp;
  }

  for (int i = 0; i < m_N; i++)
  {
	arrMtrxCorr_fi0[0] += cmparrHInv[4 + i].modul() * cmparrHInv[4 + i].modul() *mparrDisp[i];
	arrMtrxCorr_fi0[3] += cmparrHInv[12 + i].modul() * cmparrHInv[12 + i].modul() *mparrDisp[i];
	TComp z1Sopr = (*z1).Sopr() ;
	TComp cmpbSopr = cmparrHInv[12 + i].Sopr();

	TComp cmpt =  cmparrHInv[4 + i] * z1Sopr * cmpbSopr * (*z2);
	arrMtrxCorr_fi0[1] += cmpt.m_Re *mparrDisp[i];
  }
  arrMtrxCorr_fi0[0] = arrMtrxCorr_fi0[0]/ (*z1).modul()/ (*z1).modul();
  arrMtrxCorr_fi0[3] = arrMtrxCorr_fi0[3]/ (*z2).modul()/ (*z2).modul();
  arrMtrxCorr_fi0[1]  = arrMtrxCorr_fi0[1]/ (*z1).modul()/ (*z1).modul()
	   / (*z2).modul()/ (*z2).modul();
  arrMtrxCorr_fi0[2]  = arrMtrxCorr_fi0[1]  ;


  return 0;
}


// вычисление корреляц матрицы ошибок z1, z2
void TFar::calcMtrxCorr_z1_z2(TComp *pcmpSZv, TComp z1, TComp z2,  double *arrMtrxCorr )
{
  TComp cmparr_k1[4], cmparr_k2[4];
  calcArr_k(pcmpSZv,  z1, cmparr_k1);
  calcArr_k(pcmpSZv,  z2, cmparr_k2);
  arrMtrxCorr[0] = calcPhaseDisp( pcmpSZv,  z1, cmparr_k1) ;
  arrMtrxCorr[3] = calcPhaseDisp( pcmpSZv,  z2, cmparr_k2);
  arrMtrxCorr[1] = calcSmeshMoment(pcmpSZv,   z1, cmparr_k1, z2, cmparr_k2);
  arrMtrxCorr[2] = arrMtrxCorr[1];
}

// дисперсия фазового угла
double TFar::calcPhaseDisp( TComp *pcmpSZv, TComp z, TComp *cmparr_k)
{

	TComp cmpFig2(2., 0.);
	TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
	TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
	TComp cmpC = (pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2]* pcmpSZv[2]);
	TComp cmpFDeriv_Po_z = cmpFig2 * z * cmpA  + cmpB;
	double sum = 0.;
	for (int i = 0; i < 4; i++)
	{
	  sum +=  cmparr_k[i].modul() * cmparr_k[i].modul()*mparrDisp[i];

	}

	return sum / (cmpFDeriv_Po_z.modul() * cmpFDeriv_Po_z.modul() * z.modul()* z.modul());
}

TComp TFar::calcFDeriv_po_z( TComp *pcmpSZv, TComp z)
{
	TComp cmpFig2(2., 0.);
	TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
	TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
	TComp cmpC = (pcmpSZv[1] * pcmpSZv[3] - pcmpSZv[2]* pcmpSZv[2]);
	return ( cmpFig2 * z * cmpA  + cmpB);
}

// Вычисление массива коэффициентов k для корня квадратного уравнения
//INPUT:
//pcmpSZv[4] - массив измерений
//z - корень
// OUTPUT:
// cmparr_k[4] - массив коэффициентов
void   TFar::calcArr_k(TComp *pcmpSZv, TComp z, TComp *cmparr_k)
{
 TComp cmpA = (pcmpSZv[0] * pcmpSZv[2] - pcmpSZv[1]* pcmpSZv[1]);
 TComp cmpB = (pcmpSZv[1] * pcmpSZv[2] - pcmpSZv[0]* pcmpSZv[3]);
 TComp cmpFig1(1.,0.);
 TComp cmpFig2(2.,0.);
 TComp cmpMinFig2(-2.,0.);

 cmparr_k[0] =   pcmpSZv[2] * z* z - (pcmpSZv[3] * z);

 cmparr_k[1] = cmpMinFig2 * pcmpSZv[1] * z*z + pcmpSZv[2] * z +  pcmpSZv[3] ;

 cmparr_k[2] = pcmpSZv[0] * z* z + pcmpSZv[1] * z   - (cmpFig2 * pcmpSZv[2]);

 cmparr_k[3] =  pcmpSZv[1]  -  pcmpSZv[0] * z;

}

// сиешанный момоент ошибок фазовых углов
double TFar::calcSmeshMoment( TComp *pcmpSZv,  TComp z1,TComp* cmparr_k1,TComp z2, TComp* cmparr_k2)

{
	double arr_k1[4] = {0.},arr_k2[4] = {0.}, arrsum[4] ={0.}, arrT0[4] = {0.}, arrT1[4] = {0.};
	for (int i = 0; i < 4; i++)
	{
	 doArr(cmparr_k1[i], arr_k1);
	 doArr(cmparr_k2[i], arr_k2);
	 MtrxMultMatrxTransp(arr_k1,2, 2, arr_k2, 2, arrT0) ;
	 MtrxSumMatrx(arrT0, arrsum,2, 2, arrT1) ;
	 memcpy( arrsum, arrT1, 4 * sizeof(double));
	}

	double arrz1[2] ={0.},arrz2[2] ={0.}, arrFDeriv1[4] = {0.}, arrFDeriv2[4] = {0.}, arrFDerivT1[4] = {0.}, arrFDerivT2[4] = {0.};
	TComp cmpDeriv1 = calcFDeriv_po_z( pcmpSZv, z1);
	TComp cmpDeriv2 = calcFDeriv_po_z( pcmpSZv, z2);
	TComp cmpMult = z1 * z2 * cmpDeriv1 * cmpDeriv2;
	doArr(cmpDeriv1, arrFDerivT1);
	MatrxDivideScalar(arrFDerivT1, 2, 2, cmpMult.modul()  ,arrFDeriv1);
	doArr(cmpDeriv2, arrFDerivT2);
	MatrxDivideScalar(arrFDerivT2, 2, 2, cmpMult.modul()  ,arrFDeriv2);

	arrz1[0] = -z1.m_Im;
	arrz1[1] =  z1.m_Re;
	arrz2[0] = -z2.m_Im;
	arrz2[1] =  z2.m_Re;
	double arrT2[2] ={0.}, arrT3[2] ={0.}, arrT4[2] ={0.};
	MtrxMultMatrxTransp(arrz1, 1, 2, arrFDeriv1,2, arrT2) ;
	MtrxMultMatrx(arrT2, 1, 2, arrsum,2, arrT3) ;
	MtrxMultMatrx(arrT3, 1, 2, arrFDeriv2,2, arrT4) ;
	double rez =0.;
	MtrxMultMatrx(arrT4, 1, 2, arrz2,2, &rez) ;


	return rez;
}

void TFar::doArr(TComp cmpInp, double *arrRez)
{
  arrRez[0] =  cmpInp.m_Re;
  arrRez[1] = -cmpInp.m_Im ;
  arrRez[2] =  cmpInp.m_Im ;
  arrRez[3] =  cmpInp.m_Re;
}


//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // arrDisp[4] -  массив диспенрсий шума в строковых диаграммах
 //        arrDisp[i] - Дисперосия шума в строке с номером i
 // OUTPUT:
 // * valEstAngTarg - оценка угла цели
 // *valEstAngAntp -  оценка угла антипода
 // *cmpKTarg - оценка коеффиц отражения цели
 // *cmpKAntp - оценка коеффиц отражения артипода
 // arrMtrxCorr - корреляционная матирица ошибок оценивания углов
 int   TFar::fncEstimateMsd(TComp *pcmpSZv,  double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr )
{
  TComp z1(0.,0.), z2(0.,0.);
  double arrMtrxCorr_fi[4]= {0.};
  solvQuadrEqMeth(pcmpSZv, &z1, &z2 , arrMtrxCorr_fi );



	double ph1 = z1.phase();
	double ph2 = z2.phase();
	if (ph2 > ph1)
	{
		swap(&ph1, &ph2);
		swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
		swap(&z1,&z2);
	}
	double valCoef = arrMtrxCorr_fi[0]/ arrMtrxCorr_fi[3];


	if (valCoef > 4.  )
	{
	  swap(&ph1, &ph2);
	  swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
	  swap(&z1,&z2);
	}
	else
	{
	  if( valCoef > 1./4.)
	  {
	   //	if (ph1 > M_PI/ 4. * 1.2 )//&& (fabs(ph2)< M_PI/ 4.  ) )
		if (ph1 > M_PI/ 2.  )
		{
		ph1 -= M_PI * 2.;
		}

		//if (ph2 > M_PI/ 4. * 1.2 )// (fabs(ph1)< M_PI/ 4. ) )
		if (ph2 > M_PI/ 2. )
		{
		ph2 -= M_PI * 2.;
		}
		if( ph1< ph2)
		{
		swap(&ph1, &ph2);
		swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
		swap(&z1,&z2);

		}
	  }
    }

	// оценки углов:
	const double VAL_WAVE_CONST = 2. * M_PI * m_D / mFaceta.mLambda;

	*valEstAngTarg = asin( ph1/VAL_WAVE_CONST);
	*valEstAngAntp  = asin(ph2/VAL_WAVE_CONST);

	TComp cmpY1 = ( pcmpSZv[ 3] - pcmpSZv[2]* z2)/z1/z1/(z1-z2);
	TComp cmpY2 = (pcmpSZv[2] * z1 - pcmpSZv[3])/z2/z2/(z1-z2);

	TComp cmpFTarg ( mFaceta.fncFFaceta (*valEstAngTarg)/4., 0.);
	TComp cmpFAntp ( mFaceta.fncFFaceta (*valEstAngAntp )/4., 0.);
	*cmpKTarg= cmpY1 / cmpFTarg;
	*cmpKAntp = cmpY2 / cmpFAntp ;

	//
	// пересчет корррел матрицы
	const double constdTarg = sqrt(1. - ph1 * ph1/ VAL_WAVE_CONST/VAL_WAVE_CONST );
	const double constdAntp = sqrt(1. - ph2 * ph2/ VAL_WAVE_CONST/VAL_WAVE_CONST );

	arrMtrxCorr[0] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdTarg * arrMtrxCorr_fi[0]  ;
	arrMtrxCorr[3] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdAntp/ constdAntp * arrMtrxCorr_fi[3] ;
	arrMtrxCorr[1] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdAntp * arrMtrxCorr_fi[1];
	arrMtrxCorr[2] = arrMtrxCorr[1];

  return 0;
}


// ПЕРЕГРУЖЕННАЯ
//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // arrDisp[4] -  массив диспенрсий шума в строковых диаграммах
 //        arrDisp[i] - Дисперосия шума в строке с номером i
 // OUTPUT:
 // * valEstAngTarg - оценка угла цели
 // *valEstAngAntp -  оценка угла антипода
 // *cmpKTarg - оценка коеффиц отражения цели
 // *cmpKAntp - оценка коеффиц отражения артипода
 // arrMtrxCorr - корреляционная матирица ошибок оценивания углов
 // cmpZ1, cmpZ2 - комплексные углы цели и антипода
 int   TFar::fncEstimateMsd(TComp *pcmpSZv,  double *valEstAngTarg, double *valEstAngAntp
	  , TComp *cmpKTarg , TComp *cmpKAntp , double *arrMtrxCorr , TComp *cmpZ1 , TComp *cmpZ2)
{
  TComp z1(0.,0.), z2(0.,0.);
  double arrMtrxCorr_fi[4]= {0.};
  solvQuadrEqMeth(pcmpSZv, &z1, &z2 , arrMtrxCorr_fi );



	double ph1 = z1.phase();
	double ph2 = z2.phase();
	if (ph2 > ph1)
	{
		swap(&ph1, &ph2);
		swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
		swap(&z1,&z2);
	}



 /*	if (valCoef > 4.5  )
	{
	  swap(&ph1, &ph2);
	  swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
	  swap(&z1,&z2);
	}
	else
	{
	 if( valCoef > 1./4.)
	  {
	   //	if (ph1 > M_PI/ 4. * 1.2 )//&& (fabs(ph2)< M_PI/ 4.  ) )
		if (ph1 > M_PI/ 2.  )
		{
		ph1 -= M_PI * 2.;
		}

		//if (ph2 > M_PI/ 4. * 1.2 )// (fabs(ph1)< M_PI/ 4. ) )
		if (ph2 > M_PI/ 2. )
		{
		ph2 -= M_PI * 2.;
		}
		if( ph1< ph2)
		{
		swap(&ph1, &ph2);
		swap(&arrMtrxCorr_fi[0], &arrMtrxCorr_fi[3]);
		swap(&z1,&z2);

		}
	  }
	}  */

	// оценки углов:
	const double VAL_WAVE_CONST = 2. * M_PI * m_D / mFaceta.mLambda;

	*valEstAngTarg = asin( ph1/VAL_WAVE_CONST);
	*valEstAngAntp  = asin(ph2/VAL_WAVE_CONST);

	TComp cmpY1 = ( pcmpSZv[ 3] - pcmpSZv[2]* z2)/z1/z1/(z1-z2);
	TComp cmpY2 = (pcmpSZv[2] * z1 - pcmpSZv[3])/z2/z2/(z1-z2);

	TComp cmpFTarg ( mFaceta.fncFFaceta (*valEstAngTarg)/4., 0.);
	TComp cmpFAntp ( mFaceta.fncFFaceta (*valEstAngAntp )/4., 0.);
	*cmpKTarg= cmpY1 / cmpFTarg;
	*cmpKAntp = cmpY2 / cmpFAntp ;

	//
	// пересчет корррел матрицы
	const double constdTarg = sqrt(1. - ph1 * ph1/ VAL_WAVE_CONST/VAL_WAVE_CONST );
	const double constdAntp = sqrt(1. - ph2 * ph2/ VAL_WAVE_CONST/VAL_WAVE_CONST );

	arrMtrxCorr[0] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdTarg * arrMtrxCorr_fi[0]  ;
	arrMtrxCorr[3] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdAntp/ constdAntp * arrMtrxCorr_fi[3] ;
	arrMtrxCorr[1] = 1./ VAL_WAVE_CONST/ VAL_WAVE_CONST / constdTarg/ constdAntp * arrMtrxCorr_fi[1];
	arrMtrxCorr[2] = arrMtrxCorr[1];

  *cmpZ1 = z1;
  *cmpZ2 = z2;
  return 0;
}


//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный

 // OUTPUT:
 // arrEstAng[3] - массив оценок углов цели
 // cmparrK [3] - массив оценок коеффиц отражения целей
 // arrMtrxCorr[9] - корреляционная матирица ошибок оценивания углов
 // cmparrZ- массив комплексных углов целей
 int   TFar::fncPolyn3(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ)
{
  if (m_N%6 != 0)
  {
	return -1;
  }
 // формирование замеров 6 строковых диаграмм
 TComp cmparrSZv[6];
 int qRow = m_N/6;
 for (int i = 0; i < 6; i++)
 {
   cmparrSZv[i] = TComp(0.,0.);
   for (int j =0; j < qRow; j++)
   {
	cmparrSZv[i] += pcmpSZv[ i * qRow + j];
   }
 }
 ///

 double valD = qRow * m_D;
 const double VAL_WAVE_CONST = 2. * M_PI * valD / mFaceta.mLambda;
 TFaceta FacetaUn(mFaceta.m_n * qRow  , mFaceta.m_d, mFaceta.mLambda);

 // формирование матрицы Ганкеля
 TComp cmparrA[9];
 createHankel( 6, cmparrSZv, cmparrA );
 ///

 // обращение матрицы Ганкеля
  TComp cmparrAInv[9];
  InverseMtrx3(cmparrA, cmparrAInv);
  ///

  // нахождение вектора коэффициентов полинома 3 степени
   TComp cmparr_at[3],cmparr_a[3];
   MtrxMultMatrx(cmparrAInv,3, 3, &cmparrSZv[3],1, cmparr_at);
   MatrxMultScalar(cmparr_at, 3, 1, TComp(-1.,0.), cmparr_a);
   ///

   // решение уравнения
	SolvEq3(cmparrZ,TComp(1.,0.),cmparr_a[2],cmparr_a[1],cmparr_a[0]);

   double arrPh[3]= {0.};
   for (int i =0; i < 3; i++)
   {
	arrPh[i] =  cmparrZ[i].phase();
	arrEstAng[i] = 	 asin( arrPh[i]/VAL_WAVE_CONST);
   }
   ///

   //перестановка в порядке убывания
   for (int n =0;n < 10; n++)
   {
	 for (int i = 0; i < 2; i++)
	 {
		if (arrEstAng[i+1] > arrEstAng[i])
		{
		  double temp = arrEstAng[i];
		  arrEstAng[i] = arrEstAng[i+1];
		  arrEstAng[i+1] = temp;
		  TComp cmpt = cmparrZ[i];
		  cmparrZ[i] = cmparrZ[i+1] ;
		  cmparrZ[i+1] = cmpt ;
		}
	}
   }
   // нахождение вектора y
	 // формирование матрицы MInv
	 TComp cmparrMInv1[9], cmparrMInv[9];
	 cmparrMInv1[0] = cmparrZ[1] *  cmparrZ[2] * ( cmparrZ[2]-  cmparrZ[1]);
	 cmparrMInv1[1] = ( cmparrZ[2]-  cmparrZ[1])*( cmparrZ[2]+  cmparrZ[1]) * TComp(-1.,0.);
	 cmparrMInv1[2] = ( cmparrZ[2]-  cmparrZ[1]);

	 cmparrMInv1[3] = TComp(-1.,0.) * cmparrZ[0] *  cmparrZ[2] * ( cmparrZ[2]-  cmparrZ[0]);
	 cmparrMInv1[4] = ( cmparrZ[2]-  cmparrZ[0])*( cmparrZ[2]+  cmparrZ[0]) ;
	 cmparrMInv1[5] = ( cmparrZ[2]-  cmparrZ[0]) * TComp(-1.,0.);

	 cmparrMInv1[6] = cmparrZ[0] *  cmparrZ[1] * ( cmparrZ[1]-  cmparrZ[0]);
	 cmparrMInv1[7] = ( cmparrZ[1]-  cmparrZ[0])*( cmparrZ[1]+  cmparrZ[0]) * TComp(-1.,0.);
	 cmparrMInv1[8] = ( cmparrZ[1]-  cmparrZ[0]);

	 TComp cmpd = ( cmparrZ[1]-  cmparrZ[0])*( cmparrZ[2]-  cmparrZ[1]) *( cmparrZ[2]-  cmparrZ[0]);
	 MatrxDivideScalar(cmparrMInv1, 3, 3, cmpd,cmparrMInv);
	 ///

	 TComp cmparr_y[3];
	 MtrxMultMatrx(cmparrMInv,3, 3, cmparrSZv,1, cmparr_y);
	 /// end

	 // нахождение коэффициентов отражения целей

	 for (int i=0; i < 3; i++)
	 {
	  cmparrK[i] =  cmparr_y[i] / TComp(FacetaUn.fncFFacetaApprox(arrEstAng[i]),0.) * TComp(6.,0.);

	 }
	 ///
	 // нахождение корреляционной матрицы ошибок измерения
		//  формирование вектора C
		TComp cmparrC[3] ;
		cmparrC[0] = TComp(2.,0.)*  cmparrSZv[2]*  cmparrSZv[3]*  cmparrSZv[4] - cmparrSZv[3]*  cmparrSZv[3]*  cmparrSZv[3]
		  -cmparrSZv[1]*  cmparrSZv[4]*  cmparrSZv[4] +  cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[5]
		  - cmparrSZv[2]*  cmparrSZv[2]*  cmparrSZv[5];

		cmparrC[1] =   cmparrSZv[2]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[4]
		  +cmparrSZv[0]*  cmparrSZv[4]*  cmparrSZv[4] -  cmparrSZv[4]*  cmparrSZv[2]*  cmparrSZv[2]
		  + cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[5] - cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[5];

		cmparrC[2] =   cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[3]*  cmparrSZv[2]*  cmparrSZv[2]
		  +cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[4] -  cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[4]
		  + cmparrSZv[0]*  cmparrSZv[2]*  cmparrSZv[5] - cmparrSZv[5]*  cmparrSZv[1]*  cmparrSZv[1];
		///



		// формирование дискриминанта маьтрицы Ганкеля
		TComp cmpHanDiskr =  cmparrSZv[0]*  cmparrSZv[2]*  cmparrSZv[4] + cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[3] * TComp(2.,0.)
			   - cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[4]*  cmparrSZv[1]*  cmparrSZv[1]
			   - cmparrSZv[2]*  cmparrSZv[2]*  cmparrSZv[2];
		///

		// формирование матрицы частных проихводных dC po dS
		TComp cmparr_dC_po_dS[18];
		cmparr_dC_po_dS[0] = TComp(0.,0.);
		cmparr_dC_po_dS[1] = cmparrSZv[3]*  cmparrSZv[5] -cmparrSZv[4]*  cmparrSZv[4] ;
		cmparr_dC_po_dS[2] = TComp(2.,0.) * cmparrSZv[3]*  cmparrSZv[4] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[3] = TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[4] - TComp(3.,0.) * cmparrSZv[3]*  cmparrSZv[3]
							+ cmparrSZv[1]*  cmparrSZv[5];
		cmparr_dC_po_dS[4] = TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[3] - TComp(2.,0.) * cmparrSZv[1]*  cmparrSZv[4];
		cmparr_dC_po_dS[5] =  cmparrSZv[1]*  cmparrSZv[3] -  cmparrSZv[2]*  cmparrSZv[2];

		cmparr_dC_po_dS[6] = cmparrSZv[4]*  cmparrSZv[4] -cmparrSZv[3]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[7] = cmparrSZv[2]*  cmparrSZv[5] -cmparrSZv[3]*  cmparrSZv[4] ;
		cmparr_dC_po_dS[8] = cmparrSZv[3]*  cmparrSZv[3] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[4] +cmparrSZv[1]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[9] = TComp(2.,0.) *cmparrSZv[2]*  cmparrSZv[3] -  cmparrSZv[1]*  cmparrSZv[4] -cmparrSZv[0]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[10] = TComp(2.,0.) *cmparrSZv[0]*  cmparrSZv[4] -  cmparrSZv[1]*  cmparrSZv[3] -cmparrSZv[2]*  cmparrSZv[2] ;
		cmparr_dC_po_dS[11] = cmparrSZv[1]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[3]  ;

		cmparr_dC_po_dS[12] = cmparrSZv[2]*  cmparrSZv[5] -  cmparrSZv[3]*  cmparrSZv[4]  ;
		cmparr_dC_po_dS[13] = cmparrSZv[3]*  cmparrSZv[3] +  cmparrSZv[2]*  cmparrSZv[4] - TComp(2.,0.) * cmparrSZv[1]*  cmparrSZv[5];
		cmparr_dC_po_dS[14] = cmparrSZv[1]*  cmparrSZv[4] +  cmparrSZv[0]*  cmparrSZv[5] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[3];
		cmparr_dC_po_dS[15] = TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[3] -  cmparrSZv[2]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[4];
		cmparr_dC_po_dS[16] = cmparrSZv[1]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[3] ;
		cmparr_dC_po_dS[17] = cmparrSZv[0]*  cmparrSZv[2] -  cmparrSZv[1]*  cmparrSZv[1] ;
		///

		// формирование вектора градиента дискриминанта матр ГАнкеля
		TComp cmparr_dHank_po_dS[6];
		cmparr_dHank_po_dS[0] = cmparrSZv[2]*  cmparrSZv[4] -  cmparrSZv[3]*  cmparrSZv[3] ;
		cmparr_dHank_po_dS[1] = TComp(2.,0.) *cmparrSZv[2]*  cmparrSZv[3] -  TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[4] ;
		cmparr_dHank_po_dS[2] = cmparrSZv[0]*  cmparrSZv[4] + TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[3]
							  -  TComp(3.,0.) *cmparrSZv[2]*  cmparrSZv[2] ;
		cmparr_dHank_po_dS[3] = TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[2] -  TComp(2.,0.) *cmparrSZv[0]*  cmparrSZv[3] ;
		cmparr_dHank_po_dS[4] = cmparrSZv[0]*  cmparrSZv[2] -  cmparrSZv[1]*  cmparrSZv[1] ;
		cmparr_dHank_po_dS[5] = TComp(0.,0.);
		///
         	// нахождение матрицы da po dS
		TComp cmparr_da_po_dS[18], cmparrT0[18], cmparrT1[18], cmparrT2[18];
		MatrxDivideScalar(cmparr_dC_po_dS, 3, 6, cmpHanDiskr,cmparrT0);


		MtrxMultMatrx(cmparrC,3, 1, cmparr_dHank_po_dS,6, cmparrT1);
		MatrxDivideScalar(cmparrT1, 3, 6, cmpHanDiskr * cmpHanDiskr,cmparrT2);
		MtrxMinusMatrx(cmparrT0, cmparrT2,3, 6, cmparr_da_po_dS);

		// формирование матрицы Вандермонда
		TComp cmparrVander[9];
		cmparrVander[0] = TComp(1.,0.);
		cmparrVander[1] = cmparrZ[0];
		cmparrVander[2] = cmparrZ[0] * cmparrZ[0];
		cmparrVander[3] = TComp(1.,0.);
		cmparrVander[4] = cmparrZ[1];
		cmparrVander[5] = cmparrZ[1] * cmparrZ[1];
		cmparrVander[6] = TComp(1.,0.);
		cmparrVander[7] = cmparrZ[2];
		cmparrVander[8] = cmparrZ[2] * cmparrZ[2];
		///

		// формировнаие матрицы diag(1/dp(zj))
		TComp cmparrDiag[9];
		memset( cmparrDiag, 0, 9 * sizeof(TComp));
		TComp cmpMod (cmparrZ[0].modul(),0.);
		cmparrDiag [0] = TComp(1.,0.)/ (cmparrZ[0]- cmparrZ[1])/ (cmparrZ[0]- cmparrZ[2])* cmparrZ[0].Sopr()/cmpMod/cmpMod ;
		cmpMod = TComp (cmparrZ[1].modul(),0.);
		cmparrDiag [4] = TComp(1.,0.)/ (cmparrZ[1]- cmparrZ[0])/ (cmparrZ[1]- cmparrZ[2])*cmparrZ[1].Sopr()/cmpMod/cmpMod ;
		cmpMod = TComp (cmparrZ[2].modul(),0.);
		cmparrDiag [8] = TComp(1.,0.)/ (cmparrZ[2]- cmparrZ[0])/ (cmparrZ[2]- cmparrZ[1])*cmparrZ[2].Sopr()/cmpMod/cmpMod ;
		///

		// формирование матрицы C (другой)
		TComp cmparrMatrC[18], cmparrT3[9], cmparrT4[18];
		TComp cmparrMatrB[18];
		MtrxMultMatrx(cmparrDiag,3, 3, cmparrVander,3, cmparrT3);
		MtrxMultMatrx( cmparrT3,3, 3, cmparr_da_po_dS,6, cmparrT4);
		MatrxMultScalar(cmparrT4, 3, 6, TComp(-1.,0.),cmparrMatrB);
		///

		// формирование матрицы B
	   //	TComp cmparrMatrB[18];
	  //	for (int i =0 ; i < 18; i++)
	   //	{
	// TComp cmpMod ( 1./ cmparrMatrC[i].modul()/ cmparrMatrC[i].modul(), 0.);
		// cmparrMatrB[i] = cmparrMatrC[i] *  cmpMod;
	   //	}
		///

		// формирование действительной и мнимой частей матрицы B
		double arrReB[18] = {0.}, arrImB[18] = {0.};
		for (int i = 0; i < 18; i++)
		{
		 arrReB[i] = cmparrMatrB[i].m_Re;
		 arrImB[i] = cmparrMatrB[i].m_Im;
		}
		///
		//формирование диагональной матрицы дисперсий шума
		double arrMtrxDisp [36] = {0.};
		for (int i = 0; i < 6; i++)
		{
		  for (int j = 0; j < qRow; j++)
		  {
			arrMtrxDisp[ i * 6 + i] += mparrDisp[i * qRow + j];
		  }

		}
	   ///

	   // вычисление матрицы arrKmu
	   double arrKmu[9] = {0.}, arrT10[18] = {0.}, arrT11[9] = {0.}, arrT12[18] = {0.}, arrT14[9] = {0.};
	   MtrxMultMatrx(arrReB,3, 6, arrMtrxDisp,6, arrT10);
	   MtrxMultMatrxTransp(arrT10,3, 6, arrReB,3, arrT11) ;

	   MtrxMultMatrx(arrImB,3, 6, arrMtrxDisp,6, arrT12);
	   MtrxMultMatrxTransp(arrT12,3, 6, arrImB,3, arrT14) ;

	   MtrxSumMatrx(arrT11, arrT14, 3, 3, arrKmu) ;
	   ///

	   // вычисление матрицы ошибок
	   double valCoeff = mLambda / ( 2. * M_PI * valD );
	   MatrxMultScalar(arrKmu, 3, 3, valCoeff * valCoeff,arrMtrxCorr);


  return 0;
}

// ПЕРГРУЖЕННАЯ. в ВЫХОДНЫЕ ДАННЫЕ ВКЛЮЧЕНА ФУНКЦИЯ ХИ КВАДРАТ
//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный

 // OUTPUT:
 // arrEstAng[3] - массив оценок углов цели
 // cmparrK [3] - массив оценок коеффиц отражения целей
 // arrMtrxCorr[9] - корреляционная матирица ошибок оценивания углов
 // cmparrZ- массив комплексных углов целей
 int   TFar::fncPolyn3(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ, double *pvalXiSquare)
{
  if (m_N%6 != 0)
  {
	return -1;
  }
 // формирование замеров 6 строковых диаграмм  и их дисперсий шума
 TComp cmparrSZv[6];
 double arrDispUn [6] ={0.};
 int qRow = m_N/6;
 for (int i = 0; i < 6; i++)
 {
   cmparrSZv[i] = TComp(0.,0.);
   for (int j =0; j < qRow; j++)
   {
	cmparrSZv[i] += pcmpSZv[ i * qRow + j];
	arrDispUn[i] += mparrDisp[ i * qRow + j];
   }
 }
 ///

 const double valD = qRow * m_D;
 const double VAL_WAVE_CONST = 2. * M_PI * valD / mFaceta.mLambda;
 TFaceta FacetaUn(mFaceta.m_n * qRow  , mFaceta.m_d, mFaceta.mLambda);

 // формирование матрицы Ганкеля
 TComp cmparrA[9];
 createHankel( 6, cmparrSZv, cmparrA );
 ///

 // обращение матрицы Ганкеля
  TComp cmparrAInv[9];
  InverseMtrx3(cmparrA, cmparrAInv);
  ///

  // нахождение вектора коэффициентов полинома 3 степени
   TComp cmparr_at[3],cmparr_a[3];
   MtrxMultMatrx(cmparrAInv,3, 3, &cmparrSZv[3],1, cmparr_at);
   MatrxMultScalar(cmparr_at, 3, 1, TComp(-1.,0.), cmparr_a);
   ///

   // решение уравнения
	SolvEq3(cmparrZ,TComp(1.,0.),cmparr_a[2],cmparr_a[1],cmparr_a[0]);

   double arrPh[3]= {0.};
   for (int i =0; i < 3; i++)
   {
	arrPh[i] =  cmparrZ[i].phase();
	arrEstAng[i] = 	 asin( arrPh[i]/VAL_WAVE_CONST);
   }
   ///

   //перестановка в порядке убывания
   for (int n =0;n < 10; n++)
   {
	 for (int i = 0; i < 2; i++)
	 {
		if (arrEstAng[i+1] > arrEstAng[i])
		{
		  double temp = arrEstAng[i];
		  arrEstAng[i] = arrEstAng[i+1];
		  arrEstAng[i+1] = temp;
		  TComp cmpt = cmparrZ[i];
		  cmparrZ[i] = cmparrZ[i+1] ;
		  cmparrZ[i+1] = cmpt ;
		}
	}
   }
   // нахождение вектора y
	 // формирование матрицы MInv
	 TComp cmparrMInv1[9], cmparrMInv[9];
	 cmparrMInv1[0] = cmparrZ[1] *  cmparrZ[2] * ( cmparrZ[2]-  cmparrZ[1]);
	 cmparrMInv1[1] = ( cmparrZ[2]-  cmparrZ[1])*( cmparrZ[2]+  cmparrZ[1]) * TComp(-1.,0.);
	 cmparrMInv1[2] = ( cmparrZ[2]-  cmparrZ[1]);

	 cmparrMInv1[3] = TComp(-1.,0.) * cmparrZ[0] *  cmparrZ[2] * ( cmparrZ[2]-  cmparrZ[0]);
	 cmparrMInv1[4] = ( cmparrZ[2]-  cmparrZ[0])*( cmparrZ[2]+  cmparrZ[0]) ;
	 cmparrMInv1[5] = ( cmparrZ[2]-  cmparrZ[0]) * TComp(-1.,0.);

	 cmparrMInv1[6] = cmparrZ[0] *  cmparrZ[1] * ( cmparrZ[1]-  cmparrZ[0]);
	 cmparrMInv1[7] = ( cmparrZ[1]-  cmparrZ[0])*( cmparrZ[1]+  cmparrZ[0]) * TComp(-1.,0.);
	 cmparrMInv1[8] = ( cmparrZ[1]-  cmparrZ[0]);

	 TComp cmpd = ( cmparrZ[1]-  cmparrZ[0])*( cmparrZ[2]-  cmparrZ[1]) *( cmparrZ[2]-  cmparrZ[0]);
	 MatrxDivideScalar(cmparrMInv1, 3, 3, cmpd,cmparrMInv);
	 ///

	 TComp cmparr_y[3];
	 MtrxMultMatrx(cmparrMInv,3, 3, cmparrSZv,1, cmparr_y);
	 /// end

	 // проверка правильности нахождения  cmparr_y
  /*	 TComp cmparrM[9];
	 cmparrM[0] = TComp(1.,0.);
	 cmparrM[1] = TComp(1.,0.);
	 cmparrM[2] = TComp(1.,0.);

	 cmparrM[3] = cmparrZ[0];
	 cmparrM[4] = cmparrZ[1];
	 cmparrM[5] = cmparrZ[2];

	 cmparrM[6] = cmparrZ[0]*cmparrZ[0];
	 cmparrM[7] = cmparrZ[1]*cmparrZ[1];
	 cmparrM[8] = cmparrZ[2]*cmparrZ[2];
        //

	 /// конец проверки
 */
	 // нахождение коэффициентов отражения целей

	   // разность расстояний между фазовым центром  первой фасеты ФАР  FacetaUn
	   // и фазовым центром первой фасеты исходной антенны
		double delD =  (valD - m_D) /2.;
		TComp cmparrUnK[3];
	 for (int i=0; i < 3; i++)
	 {
	  cmparrUnK[i] = cmparr_y[i] / TComp(FacetaUn.fncFFacetaApprox(arrEstAng[i]),0.) * TComp(6.,0.);
	  //угловая разность хода
	  TComp cmpTemp = exp_(-delD * sin(arrEstAng[i]) / mLambda * 2. * M_PI);
	  cmparrK[i] =  cmparrUnK[i] * cmpTemp;

	 }
	 ///
	 // нахождение корреляционной матрицы ошибок измерения
		//  формирование вектора C
		TComp cmparrC[3] ;
		cmparrC[0] = TComp(2.,0.)*  cmparrSZv[2]*  cmparrSZv[3]*  cmparrSZv[4] - cmparrSZv[3]*  cmparrSZv[3]*  cmparrSZv[3]
		  -cmparrSZv[1]*  cmparrSZv[4]*  cmparrSZv[4] +  cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[5]
		  - cmparrSZv[2]*  cmparrSZv[2]*  cmparrSZv[5];

		cmparrC[1] =   cmparrSZv[2]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[4]
		  +cmparrSZv[0]*  cmparrSZv[4]*  cmparrSZv[4] -  cmparrSZv[4]*  cmparrSZv[2]*  cmparrSZv[2]
		  + cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[5] - cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[5];

		cmparrC[2] =   cmparrSZv[1]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[3]*  cmparrSZv[2]*  cmparrSZv[2]
		  +cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[4] -  cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[4]
		  + cmparrSZv[0]*  cmparrSZv[2]*  cmparrSZv[5] - cmparrSZv[5]*  cmparrSZv[1]*  cmparrSZv[1];
		///



		// формирование дискриминанта маьтрицы Ганкеля
		TComp cmpHanDiskr =  cmparrSZv[0]*  cmparrSZv[2]*  cmparrSZv[4] + cmparrSZv[1]*  cmparrSZv[2]*  cmparrSZv[3] * TComp(2.,0.)
			   - cmparrSZv[0]*  cmparrSZv[3]*  cmparrSZv[3] - cmparrSZv[4]*  cmparrSZv[1]*  cmparrSZv[1]
			   - cmparrSZv[2]*  cmparrSZv[2]*  cmparrSZv[2];
		///

		// формирование матрицы частных проихводных dC po dS
		TComp cmparr_dC_po_dS[18];
		cmparr_dC_po_dS[0] = TComp(0.,0.);
		cmparr_dC_po_dS[1] = cmparrSZv[3]*  cmparrSZv[5] -cmparrSZv[4]*  cmparrSZv[4] ;
		cmparr_dC_po_dS[2] = TComp(2.,0.) * cmparrSZv[3]*  cmparrSZv[4] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[3] = TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[4] - TComp(3.,0.) * cmparrSZv[3]*  cmparrSZv[3]
							+ cmparrSZv[1]*  cmparrSZv[5];
		cmparr_dC_po_dS[4] = TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[3] - TComp(2.,0.) * cmparrSZv[1]*  cmparrSZv[4];
		cmparr_dC_po_dS[5] =  cmparrSZv[1]*  cmparrSZv[3] -  cmparrSZv[2]*  cmparrSZv[2];

		cmparr_dC_po_dS[6] = cmparrSZv[4]*  cmparrSZv[4] -cmparrSZv[3]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[7] = cmparrSZv[2]*  cmparrSZv[5] -cmparrSZv[3]*  cmparrSZv[4] ;
		cmparr_dC_po_dS[8] = cmparrSZv[3]*  cmparrSZv[3] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[4] +cmparrSZv[1]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[9] = TComp(2.,0.) *cmparrSZv[2]*  cmparrSZv[3] -  cmparrSZv[1]*  cmparrSZv[4] -cmparrSZv[0]*  cmparrSZv[5] ;
		cmparr_dC_po_dS[10] = TComp(2.,0.) *cmparrSZv[0]*  cmparrSZv[4] -  cmparrSZv[1]*  cmparrSZv[3] -cmparrSZv[2]*  cmparrSZv[2] ;
		cmparr_dC_po_dS[11] = cmparrSZv[1]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[3]  ;

		cmparr_dC_po_dS[12] = cmparrSZv[2]*  cmparrSZv[5] -  cmparrSZv[3]*  cmparrSZv[4]  ;
		cmparr_dC_po_dS[13] = cmparrSZv[3]*  cmparrSZv[3] +  cmparrSZv[2]*  cmparrSZv[4] - TComp(2.,0.) * cmparrSZv[1]*  cmparrSZv[5];
		cmparr_dC_po_dS[14] = cmparrSZv[1]*  cmparrSZv[4] +  cmparrSZv[0]*  cmparrSZv[5] - TComp(2.,0.) * cmparrSZv[2]*  cmparrSZv[3];
		cmparr_dC_po_dS[15] = TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[3] -  cmparrSZv[2]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[4];
		cmparr_dC_po_dS[16] = cmparrSZv[1]*  cmparrSZv[2] -  cmparrSZv[0]*  cmparrSZv[3] ;
		cmparr_dC_po_dS[17] = cmparrSZv[0]*  cmparrSZv[2] -  cmparrSZv[1]*  cmparrSZv[1] ;
		///

		// формирование вектора градиента дискриминанта матр ГАнкеля
		TComp cmparr_dHank_po_dS[6];
		cmparr_dHank_po_dS[0] = cmparrSZv[2]*  cmparrSZv[4] -  cmparrSZv[3]*  cmparrSZv[3] ;
		cmparr_dHank_po_dS[1] = TComp(2.,0.) *cmparrSZv[2]*  cmparrSZv[3] -  TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[4] ;
		cmparr_dHank_po_dS[2] = cmparrSZv[0]*  cmparrSZv[4] + TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[3]
							  -  TComp(3.,0.) *cmparrSZv[2]*  cmparrSZv[2] ;
		cmparr_dHank_po_dS[3] = TComp(2.,0.) *cmparrSZv[1]*  cmparrSZv[2] -  TComp(2.,0.) *cmparrSZv[0]*  cmparrSZv[3] ;
		cmparr_dHank_po_dS[4] = cmparrSZv[0]*  cmparrSZv[2] -  cmparrSZv[1]*  cmparrSZv[1] ;
		cmparr_dHank_po_dS[5] = TComp(0.,0.);
		///
         	// нахождение матрицы da po dS
		TComp cmparr_da_po_dS[18], cmparrT0[18], cmparrT1[18], cmparrT2[18];
		MatrxDivideScalar(cmparr_dC_po_dS, 3, 6, cmpHanDiskr,cmparrT0);


		MtrxMultMatrx(cmparrC,3, 1, cmparr_dHank_po_dS,6, cmparrT1);
		MatrxDivideScalar(cmparrT1, 3, 6, cmpHanDiskr * cmpHanDiskr,cmparrT2);
		MtrxMinusMatrx(cmparrT0, cmparrT2,3, 6, cmparr_da_po_dS);

		// формирование матрицы Вандермонда
		TComp cmparrVander[9];
		cmparrVander[0] = TComp(1.,0.);
		cmparrVander[1] = cmparrZ[0];
		cmparrVander[2] = cmparrZ[0] * cmparrZ[0];
		cmparrVander[3] = TComp(1.,0.);
		cmparrVander[4] = cmparrZ[1];
		cmparrVander[5] = cmparrZ[1] * cmparrZ[1];
		cmparrVander[6] = TComp(1.,0.);
		cmparrVander[7] = cmparrZ[2];
		cmparrVander[8] = cmparrZ[2] * cmparrZ[2];
		///

		// формировнаие матрицы diag(1/dp(zj))
		TComp cmparrDiag[9];
		memset( cmparrDiag, 0, 9 * sizeof(TComp));
		TComp cmpMod (cmparrZ[0].modul(),0.);
		cmparrDiag [0] = TComp(1.,0.)/ (cmparrZ[0]- cmparrZ[1])/ (cmparrZ[0]- cmparrZ[2])* cmparrZ[0].Sopr()/cmpMod/cmpMod ;
		cmpMod = TComp (cmparrZ[1].modul(),0.);
		cmparrDiag [4] = TComp(1.,0.)/ (cmparrZ[1]- cmparrZ[0])/ (cmparrZ[1]- cmparrZ[2])*cmparrZ[1].Sopr()/cmpMod/cmpMod ;
		cmpMod = TComp (cmparrZ[2].modul(),0.);
		cmparrDiag [8] = TComp(1.,0.)/ (cmparrZ[2]- cmparrZ[0])/ (cmparrZ[2]- cmparrZ[1])*cmparrZ[2].Sopr()/cmpMod/cmpMod ;
		///

		// формирование матрицы C (другой)
		TComp cmparrMatrC[18], cmparrT3[9], cmparrT4[18];
		TComp cmparrMatrB[18];
		MtrxMultMatrx(cmparrDiag,3, 3, cmparrVander,3, cmparrT3);
		MtrxMultMatrx( cmparrT3,3, 3, cmparr_da_po_dS,6, cmparrT4);
		MatrxMultScalar(cmparrT4, 3, 6, TComp(-1.,0.),cmparrMatrB);
		///

		// формирование матрицы B
	   //	TComp cmparrMatrB[18];
	  //	for (int i =0 ; i < 18; i++)
	   //	{
	// TComp cmpMod ( 1./ cmparrMatrC[i].modul()/ cmparrMatrC[i].modul(), 0.);
		// cmparrMatrB[i] = cmparrMatrC[i] *  cmpMod;
	   //	}
		///

		// формирование действительной и мнимой частей матрицы B
		double arrReB[18] = {0.}, arrImB[18] = {0.};
		for (int i = 0; i < 18; i++)
		{
		 arrReB[i] = cmparrMatrB[i].m_Re;
		 arrImB[i] = cmparrMatrB[i].m_Im;
		}
		///
		//формирование диагональной матрицы дисперсий шума
		double arrMtrxDisp [36] = {0.};
		for (int i = 0; i < 6; i++)
		{
		  for (int j = 0; j < qRow; j++)
		  {
			arrMtrxDisp[ i * 6 + i] += mparrDisp[i * qRow + j];
		  }

		}
	   ///

	   // вычисление матрицы arrKmu
	   double arrKmu[9] = {0.}, arrT10[18] = {0.}, arrT11[9] = {0.}, arrT12[18] = {0.}, arrT14[9] = {0.};
	   MtrxMultMatrx(arrReB,3, 6, arrMtrxDisp,6, arrT10);
	   MtrxMultMatrxTransp(arrT10,3, 6, arrReB,3, arrT11) ;

	   MtrxMultMatrx(arrImB,3, 6, arrMtrxDisp,6, arrT12);
	   MtrxMultMatrxTransp(arrT12,3, 6, arrImB,3, arrT14) ;

	   MtrxSumMatrx(arrT11, arrT14, 3, 3, arrKmu) ;
	   ///

	   // вычисление матрицы ошибок
	   double valCoeff = mLambda / ( 2. * M_PI * valD );
	   MatrxMultScalar(arrKmu, 3, 3, valCoeff * valCoeff,arrMtrxCorr);


	   int quantTarg = 3;
	   *pvalXiSquare = calcXiSquare( quantTarg, arrEstAng, cmparrK, pcmpSZv) ;
	   TFar  FarUn(6, valD, mLambda, FacetaUn, arrDispUn);
	   *pvalXiSquare = FarUn.calcXiSquare( quantTarg, arrEstAng, cmparrUnK, cmparrSZv) ;

  return 0;
}


//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный

 // OUTPUT:
 // arrEstAng[3] - массив оценок углов цели
 // cmparrK [3] - массив оценок коеффиц отражения целей
 // arrMtrxCorr[9] - корреляционная матирица ошибок оценивания углов
 // cmparrZ- массив комплексных углов целей
 int   TFar::fncPolyn2(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ)
{
  if (m_N%4 != 0)
  {
	return -1;
  }
 // формирование замеров 4 строковых диаграмм
 TComp cmparrSZv[4];
 int qRow = m_N/4;
 double arrDisp[4] ={0.};
 for (int i = 0; i < 4; i++)
 {
   cmparrSZv[i] = TComp(0.,0.);
   for (int j =0; j < qRow; j++)
   {
	cmparrSZv[i] += pcmpSZv[ i * qRow + j];
	arrDisp[i] += mparrDisp[ i * qRow + j];
   }
 }
 ///

 double valD = qRow * m_D;

 TFaceta FacetaUn(mFaceta.m_n * qRow  , mFaceta.m_d, mFaceta.mLambda);
 TFar FarUn(4,valD , mLambda ,FacetaUn,arrDisp);
 return FarUn.fncEstimateMsd(cmparrSZv,  arrEstAng, &arrEstAng[1]
	  , cmparrK , &cmparrK[1] , arrMtrxCorr , cmparrZ , &cmparrZ[1]);

}

// ПЕРЕГРУЖЕННАЯ, В ЧИСЛО ВЫХОДНЫХ ПАРАМЕТОВ ВКЛЮЧЕН ХИ КВАДРАТ
//INPUT:
 // pcmpSZv - массив измерений дианрамм, комплексный
 // OUTPUT:
 // arrEstAng[3] - массив оценок углов цели
 // cmparrK [3] - массив оценок коеффиц отражения целей
 // arrMtrxCorr[9] - корреляционная матирица ошибок оценивания углов
 // cmparrZ- массив комплексных углов целей
 int   TFar::fncPolyn2(TComp *pcmpSZv,  double *arrEstAng
	  , TComp *cmparrK, double *arrMtrxCorr , TComp *cmparrZ,   double *pvalXiSquare)
{
  if (m_N%4 != 0)
  {
	return -1;
  }
 // формирование замеров 4 строковых диаграмм
 TComp cmparrSZv[4];
 int qRow = m_N/4;
 double arrDisp[4] ={0.};
 for (int i = 0; i < 4; i++)
 {
   cmparrSZv[i] = TComp(0.,0.);
   for (int j =0; j < qRow; j++)
   {
	cmparrSZv[i] += pcmpSZv[ i * qRow + j];
	arrDisp[i] += mparrDisp[ i * qRow + j];
   }
 }
 ///

 double valD = qRow * m_D;

 TFaceta FacetaUn(mFaceta.m_n * qRow  , mFaceta.m_d, mFaceta.mLambda);
 TFar FarUn(4,valD , mLambda ,FacetaUn,arrDisp);
 int ireturn = FarUn.fncEstimateMsd(cmparrSZv,  arrEstAng, &arrEstAng[1]
	  , cmparrK , &cmparrK[1] , arrMtrxCorr , cmparrZ , &cmparrZ[1]);
  // разность расстояний между фазовым центром  первой фасеты ФАР  FacetaUn
	   // и фазовым центром первой фасеты исходной антенны
		double delD =  (valD - m_D) /2.;
	 for (int i=0; i < 2; i++)
	 {
	  //угловая разность хода
	  TComp cmpTemp = exp_(-delD * sin(arrEstAng[i]) / mLambda * 2. * M_PI);
	  cmparrK[i] =  cmparrK[i]*cmpTemp;

	 }
 *pvalXiSquare = calcXiSquare(2, arrEstAng,cmparrK , pcmpSZv) ;
 return ireturn;

}


// коррекция массива дисперсий на кроффициенты
void TFar::fncCorrectDisp(TComp *cmparrCorectCoef)
{
	for (int i = 0; i < m_N; i++)
	{
	 mparrDisp[i] *=  cmparrCorectCoef[i].modul();
	}
}
// статистическая обработка измерений
// нахождение векиора поправочных коеффициентов
// INPUT
// quantMeas - к-во измерений  вынутых из   log файла
// pcmparrMeas - массив измерений  , исходный, правильно перевернутый
// cmparrCorectCoef - существующие коэффициенты, учитывающие к-во АМ в строках
// у нас это 0 и 3-я строки 4/3
// OUTPUT:
// cmparrCorectCoefStat  - дополнительные коэффициенты - то, что мы ищем
void TFar::calcArrayStatCoeff(int quantMeas,TComp * pcmparrMeas, TComp *cmparrCorectCoef, TComp * cmparrCorectCoefStat)
{
// текущий массив статитсических корректирующих коеффициентов  m_N =4
TComp *cmparrRo0 = ( TComp *)malloc( m_N * sizeof(TComp));

// текущий массив результирующих корректирующих коеффициентов
TComp *cmparrK = ( TComp *)malloc( m_N * sizeof(TComp));

// скоректированны массив замеров
TComp *pcmparrCorrectMeas= ( TComp *)malloc( quantMeas *m_N * sizeof(TComp));

// массив развернутых замеров
TComp *pcmparrQZv= ( TComp *)malloc( quantMeas *m_N * sizeof(TComp));

// инициализация массива  cmparrRo0
for (int i =0; i < m_N; i++)
{
 cmparrRo0[i] = TComp(1.,0.);
}
///

//массив средгих величин модкля замера
 double *parrAverModSZv = ( double *)malloc( m_N * sizeof(double));
 // транспонированный массив pcmparrQZv
 TComp *pcmparrQZvT = ( TComp *)malloc( quantMeas *m_N * sizeof(TComp));

 // вычисление среднего модуля по строке  - массив parrAverModSZv
 fncMakeArrK( cmparrCorectCoef,cmparrRo0, cmparrK);
 fncCorrectMeas(cmparrK,quantMeas, pcmparrMeas, pcmparrCorrectMeas);
  for (int i = 0; i < quantMeas ; i++)
  {
   for (int j =0; j < m_N; j++ )
   {
	 parrAverModSZv[j] += pcmparrCorrectMeas [ i * m_N + j].modul()/ ((double)quantMeas);
   }
  }
  ///
 // нахождение строки с максимальным средним модулем по этому АМ будем "равнять "  остальные
 // будем считать, что эта строка со всеми рабочими модулями
 int numArgMax = -1;
  MaxDoubleArray(parrAverModSZv
	, m_N, &numArgMax) ;
	///

// итерационный процесс нахождения корректирующих коэффициентов
// реально требует 2-3 итрации
int iC =0;
for ( iC =0; iC < 100; iC++)
{
  fncMakeArrK( cmparrCorectCoef,cmparrRo0, cmparrK); //  вычисление текущего результирующего массива коеффициентов
  fncCorrectMeas(cmparrK,quantMeas, pcmparrMeas, pcmparrCorrectMeas);  // коррекция массива замеров на эти коэффициенты
  // поворот замеров
	 for (int i = 0; i < quantMeas; i++)
	 {
	   // оценка обобщенного угла РСМ  для i-го замера
	   double valXi2 = 0.,valEstRSM = 0., valRSMDisp = 0.;
	   fncMeasureRSMProcessing(&pcmparrCorrectMeas[i * m_N], &valXi2
		 , &valEstRSM, &valRSMDisp  );
	   double valEstMuCur =  sin(valEstRSM) * m_D * 2. * M_PI/ mLambda;
	   ///
		 // повороты замеров строк i-го замера
			 for (int j = 0; j < m_N; j++)
			 {
			   pcmparrQZv[i *m_N + j] =  exp_(-valEstMuCur * ((double)j))
				  *pcmparrCorrectMeas [i *m_N + j];
			 }
		  ///
	 }
  ///
  // транспонирование массива pcmparrQZv.
  // матрица pcmparrQZvT - это транспонированная матрица pcmparrQZv
  //  quantMeas -  к-во строк матрицы pcmparrQZv
  // m_N - -  к-во столбцов матрицы pcmparrQZv
  MatrTransp(pcmparrQZv, quantMeas, m_N, pcmparrQZvT);
  ///


// вычисление корректирующих коеффициентов
cmparrCorectCoefStat[numArgMax] = TComp(1.,0.);
for (int i =0; i <  m_N; i++)
{
  if (i == numArgMax)
  {
    continue;
  }
  cmparrCorectCoefStat[i] = calcCorrectCoeff(pcmparrQZvT,quantMeas,numArgMax, i);
}
 ///

 // вычисление невязки
 double valNev = 0.;
 TComp cmp1(1.,0.);
 for (int i =0; i <  m_N; i++)
{
  TComp cmpT =  cmparrCorectCoefStat[i] -  cmp1;
  valNev += cmpT.modul() ;
}
 //memcpy(cmparrRo0  ,  cmparrCorectCoefStat, m_N * sizeof(TComp));
 for (int ii = 0; ii < m_N; ii++)
 {
   cmparrRo0[ii] *=  cmparrCorectCoefStat[ii];
 }
 if (valNev < 0.01 *((double)m_N))
 {
   break;
 }

}
memcpy( cmparrCorectCoefStat, cmparrRo0, m_N * sizeof(TComp));
free(cmparrRo0);

free(cmparrK);
free(pcmparrCorrectMeas);

free(pcmparrQZv);
free(pcmparrQZvT);

free (parrAverModSZv);


}


// pcmparrQZvT -транспонированный массив повернутых замеров (то есть, в нем замеры храняться по строкам)
// quantMeas - к-во замеров
// numArgMax - номер строки (строковой диаграммы) по которой будет осуществляться выравнивание
// numRow - номер строки, которая будет подравниваться под строку с номером numArgMax
//
//
//
TComp TFar::calcCorrectCoeff( TComp *pcmparrQZvT, int quantMeas,int numArgMax, int numRow)
{
 TComp cmpSum_qZvSigMin2  (0.,0.),cmpSumSigMin2  (0.,0.);
 for (int i=0; i < quantMeas; i++)
 {
   
   double valSigSq  = (mparrDisp [numArgMax]  + mparrDisp [numRow])
	/ pcmparrQZvT[numArgMax * quantMeas + i].modul()
	/ pcmparrQZvT[numArgMax * quantMeas + i].modul();

	TComp cmpSigMin2(1./ valSigSq, 0.);
   cmpSumSigMin2 += cmpSigMin2;

   TComp cmpTemp = pcmparrQZvT[numRow * quantMeas + i]/pcmparrQZvT[numArgMax * quantMeas + i];
   TComp cmp_qZv = cmpTemp.Ln();
   cmpSum_qZvSigMin2 +=  cmp_qZv * cmpSigMin2;
 }
 TComp cmpa =  cmpSum_qZvSigMin2/ cmpSumSigMin2;
 TComp cmpRez = exp_(cmpa);
 TComp cmp1(1.,0.);
 cmpRez = cmp1/ cmpRez;
 

 return cmpRez;
}
 //  cmparrK[i]= cmparrCorectCoef[i] * cmparrRo0[i]
 void TFar::fncMakeArrK(TComp *cmparrCorectCoef,TComp *cmparrRo0, TComp *cmparrK)
 {
   for (int i =0; i < m_N; i++)
   {
	 cmparrK[i] = cmparrCorectCoef[i] * cmparrRo0[i];
   }
 }
// коррекция исходных замеров с учетом массива корректирующих коэффициентов
// каждая строка умножается на свой коэффициент
void TFar:: fncCorrectMeas(TComp *cmparrK,int quantMeas,TComp * pcmparrMeas
   ,TComp * pcmparrCorrectMeas)
{
   for (int i =0; i < quantMeas; i++)
   {
	  for (int j =0; j < m_N; j++)
	  {
		pcmparrCorrectMeas[i * m_N + j] =  pcmparrMeas[i * m_N + j] *cmparrK[j];
	  }
   }
}

// имитация массива замера строк
void TFar::ImitateMeasureArray(double alfUMTrg,TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp, TComp *pcmpS, TComp *pcmpSZv)
{
	const double VAL_WAVE_CONST = 2. * M_PI * m_D / mLambda;
	TComp cmpFTarg( mFaceta.fncFFacetaApprox (alfUMTrg)/ ((double)m_N) ,0.);
	TComp cmpFAntp( mFaceta.fncFFacetaApprox (alfUMAntp)/ ((double)m_N) ,0.);

	for (int i =0; i < m_N ; i++)
	{
	  TComp cmpZRotTarg = exp_ (sin(alfUMTrg) * ((double)i)  * VAL_WAVE_CONST);
	  TComp cmpZRotAntp = exp_ (sin(alfUMAntp) * ((double)i) * VAL_WAVE_CONST);
	  pcmpS[i] = cmpZRotTarg *(cmpKTarg * cmpFTarg) +  cmpZRotAntp *(cmpKAntp * cmpFAntp);
		pcmpSZv[i].m_Re =( pcmpS[i].m_Re +  getGauss(0., sqrt(mparrDisp[i] /2.) ));
		pcmpSZv[i].m_Im =( pcmpS[i].m_Im +  getGauss(0., sqrt(mparrDisp[i])/2. ));

	}
}

// имитация массива замера строк ПЕРЕГРУЖЕННАЯ!!!!
//  INPUT:
// quantTarg - к-во целей
// arrUMtrg  [quantTarg] - истинные углы целей
// cmparrKTarg [quantTarg] - коэффиц отражения целей
// OUTPUT:
// pcmpS[m_N]  - массив измерений строковых диаграмм без шумов
// pcmpSZv[m_N] - массив измерений строковых диаграмм с шумами
void TFar::ImitateMeasureArray(int quantTarg, double *arrUMtrg,TComp  *cmparrKTarg
			, TComp *pcmpS, TComp *pcmpSZv)
{
	const double VAL_WAVE_CONST = 2. * M_PI * m_D / mLambda;
	TComp *cmparrFTarg = new TComp[quantTarg];
	for (int i =0 ; i < quantTarg; i++)
	{
	   cmparrFTarg [i] =  TComp ( mFaceta.fncFFacetaApprox (arrUMtrg[i])/ m_N ,0.);
	}
	for (int i =0; i < m_N ; i++)
	{
	  pcmpS[i]= TComp(0.,0.);
	  for (int j = 0; j < quantTarg; j++)
	  {
		TComp cmpZRotTarg = exp_ (sin(arrUMtrg[j]) * ((double)i)  * VAL_WAVE_CONST);
		pcmpS[i] += cmpZRotTarg * cmparrKTarg[j] * cmparrFTarg [j];
	  }
		pcmpSZv[i].m_Re =  pcmpS[i].m_Re +  getGauss(0., sqrt(mparrDisp[i] /2.) );
		pcmpSZv[i].m_Im =  pcmpS[i].m_Im +  getGauss(0., sqrt(mparrDisp[i]/2.) );

	}

	delete []cmparrFTarg;


}


// Вычисление функции хи квадрат
//  INPUT:
// quantTarg - к-во целей в гипотезе
// arrEstUMtrg  [quantTarg] - оценки углов целей
// cmparrEstKTarg [quantTarg] - оценки коэффиц отражения целей
// OUTPUT:
// pcmpS[m_N]  - массив измерений строковых диаграмм без шумов
// pcmpSZv[m_N] - массив измерений строковых диаграмм с шумами
double TFar::calcXiSquare(int quantTarg, double *arrEstUMtrg,TComp  *cmparrEstKTarg , TComp *pcmpSZv)
{
	double valXiSquare = 0.;
	const double VAL_WAVE_CONST = 2. * M_PI * m_D / mLambda;
	//TComp *cmparrDel = new TComp[m_N];
   //	memcpy(cmparrDel, pcmpSZv, sizeof(TComp)* m_N);

	TComp *cmparrFTarg = new TComp[quantTarg];
	for (int i =0 ; i < quantTarg; i++)
	{
	   cmparrFTarg [i] =  TComp ( mFaceta.fncFFacetaApprox (arrEstUMtrg[i])/ m_N ,0.);
	}
	for (int i =0; i < m_N ; i++)
	{
	  TComp cmpDel = pcmpSZv[i];
	  TComp cmpSum(0.,0.);
	  for (int j = 0; j < quantTarg; j++)
	  {
		TComp cmpZRotTarg = exp_ (sin(arrEstUMtrg[j]) * ((double)i)  * VAL_WAVE_CONST);
		cmpDel -= cmpZRotTarg * cmparrEstKTarg[j] * cmparrFTarg [j];
		cmpSum +=  cmpZRotTarg * cmparrEstKTarg[j] * cmparrFTarg [j];
	  }
	  valXiSquare += cmpDel.modul() * cmpDel.modul() /mparrDisp[i];


	}

   //	delete []cmparrDel;
	delete []cmparrFTarg;
	return valXiSquare;
}


// построенгие графиков метода статистических испытаний для фиксированных значений
// угла цели и антипода
//wchFoldName  - путь к папаеке с отчетом
// NIsp  - число испытаний
// alfUMTrg, cmpKTarg - ум цели, коеффиц отражения цели
// alfUMAntp, cmpKAntp - ум антип, коеффиц отражения антип
 void TFar::makeMonteCarloGraphs(wchar_t *wchFoldName1, const int NIsp, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp)
{

  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[9] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);

  int iNCols = 13;
  double * parrBuff = (double *) malloc(iNCols * NIsp * sizeof(double));
  TComp z1(0.,0.),z2(0.,0.);
  TURPointXY *ppntarrZ1Izm = new TURPointXY[ NIsp];
  TURPointXY *ppntarrZ2Izm = new TURPointXY[ NIsp];
  double arrEstAng[3] = {0.};
  TComp cmparrK[3];
  TComp cmparrZ[3];


  for (int i = 0; i < NIsp; i++)
	 {
		ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

		double valXiSquare= 10000000.;
		fncPolyn2(pcmpSZv,  arrEstAng, cmparrK, arrMtrxCorr , cmparrZ, &valXiSquare);
		z1 = cmparrZ[0];
		z2 = cmparrZ[1];
		valEstAngTarg =  arrEstAng[0];
		valEstAngAntp =  arrEstAng[1];
		cmpEstKTarg =  cmparrK[0];
		cmpEstKAntp  =  cmparrK[1];
	   double valXi2  =0., valEstRSM =0.,  valRSMDisp =0.;
	  fncMeasureRSMProcessing(pcmpSZv, &valXi2, &valEstRSM, &valRSMDisp  );



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
	  parrBuff[ i * iNCols + 12] = valEstRSM * 1000.;

	  ppntarrZ1Izm[i].X =  z1.m_Re;
	  ppntarrZ1Izm[i].Y =  z1.m_Im;
	  ppntarrZ2Izm[i].X =  z2.m_Re;
	  ppntarrZ2Izm[i].Y =  z2.m_Im;



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
	 wcscpy(&wcharrFileNames[30 * 12],L"EstRSM");


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
 delete []wcharrFileNames;
const double VAL_WAVE_CONST = 2. * M_PI * m_D / mFaceta.mLambda;
// построение  окружности с центром в 0
TURPointXY pointZero(0.,0.);
TURPolygon plgCircle0 = TURPolygon::fncCreateCircle(pointZero, 1., 500);
wcscpy(  wchAxesFileName,  wchFoldName);
wcscat(wchAxesFileName, L"Circle0.shp");
plgCircle0.WriteSetSHPFiles(wchAxesFileName,&plgCircle0, 1) ;
	 ///
	 ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

		fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
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



// ПЕРЕГРУЖЕННАЯ ДЛЯ СРАВНЕНИЯ ГИПОТЕЗ 2 И 3 ЦЕЛЕЙ
// построенгие графиков метода статистических испытаний для фиксированных значений
// угла цели и антипода
//INPUT:
//wchFoldName  - путь к папаеке с отчетом
// NIsp  - число испытаний
// QUantTarg - истинное число целей
// arrUMtrg - массив углов целей
// cmparrKTarg- массив коэффиц отражения целей
// В соответствии с исх данными происходит моделирование измерений ФАР
// результаты моделирования подаются на 2 алгоритма соотвнетствующих 2 гипотезам
// - гипотезе 2 целей и 3 целей
//
//
 void TFar::makeMonteCarloGraphs_for_2_and_3_Targs(wchar_t *wchFoldName1, const int NIsp
 ,const int QUantTarg, double *arrUMtrg, TComp *cmparrKTarg )
{

  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double arrMtrxCorr[4] ={0.},arrMtrxCorr3[9] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);

  int iNCols = 21;
  double * parrBuff = (double *) malloc(iNCols * NIsp * sizeof(double));
  memset(parrBuff, 0, iNCols * NIsp * sizeof(double));
  TComp z1(0.,0.),z2(0.,0.);

  double arrEstAng[2] = {0.},arrEstAng3[3] = {0.};
  TComp cmparrEstK[2], cmparrEstK3[3];
  TComp cmparrZ[2],cmparrZ3[3];
  for (int i = 0; i < NIsp; i++)
	 {

		ImitateMeasureArray(QUantTarg, arrUMtrg,cmparrKTarg, pcmpS, pcmpSZv);

		double valXiSquare = 0.;
		fncPolyn2(pcmpSZv,  arrEstAng, cmparrEstK, arrMtrxCorr , cmparrZ, &valXiSquare);
		z1 = cmparrZ[0];
		z2 = cmparrZ[1];
		fncPolyn3(pcmpSZv,  arrEstAng3, cmparrEstK3, arrMtrxCorr3 , cmparrZ3, &valXiSquare);
	   double valXi2  =0., valEstRSM =0., valRSMDisp =0.;
	  fncMeasureRSMProcessing(pcmpSZv, &valXi2, &valEstRSM, &valRSMDisp  );
	   parrBuff[ i * iNCols] = (double)(i);  // номер испытания
	  parrBuff[ i * iNCols + 1] = arrUMtrg[0] * 1000.;  // уголл цели
	  parrBuff[ i * iNCols + 2] =  arrUMtrg[1] * 1000.; // угол АНТП 1
	  parrBuff[ i * iNCols + 3] =  arrUMtrg[2] * 1000.; // угол АНТП 2

	  parrBuff[ i * iNCols + 4] = arrEstAng[0] * 1000.; // оценнка ум цели алгоритмом гипотезы 2 целей
	  parrBuff[ i * iNCols + 5] = arrEstAng[1] * 1000.; // оценнка ум антипода алгоритмом гипотезы 2 целей
	  parrBuff[ i * iNCols + 6] = sqrt(arrMtrxCorr[0]) * 1000.; // скз ошибки оценивания угла  цели по гипотезе 2-х целей
	  parrBuff[ i * iNCols + 7] = sqrt(arrMtrxCorr[3]) * 1000.; // скз ошибки оценивания угла  АНПД по гипотезе 2-х целей

	  double valSigDel= sqrt(arrMtrxCorr[0] - 2.* arrMtrxCorr[1] + arrMtrxCorr[2]);
	  parrBuff[ i * iNCols + 8] = (arrEstAng[0] - arrEstAng[1] -1.65 * valSigDel) * 1000.;// функция разрешения для гипотезы 2 целей

	  // алгоритм 3 целей

	  parrBuff[ i * iNCols + 9] =  arrEstAng3[0] * 1000.; // угол цели 0
	  parrBuff[ i * iNCols + 10] =  arrEstAng3[1] * 1000.; // угол цели 1
	  parrBuff[ i * iNCols + 11] =  arrEstAng3[2] * 1000.; // угол цели 2

	  parrBuff[ i * iNCols + 12] =  arrMtrxCorr3[0] * 1000.; // скз цели 0
	  parrBuff[ i * iNCols + 13] =  arrMtrxCorr3[4] * 1000.; // скз цели 1
	  parrBuff[ i * iNCols + 14] =  arrMtrxCorr3[8] * 1000.; // скз цели 2

	  parrBuff[ i * iNCols + 15] = (arrEstAng3[0] - arrEstAng3[1]
					-1.65 * sqrt(arrMtrxCorr3[0] - 2.* arrMtrxCorr3[1] + arrMtrxCorr3[4])) * 1000.;// функция разрешения 0 и 1 целей
	  parrBuff[ i * iNCols + 16] = (arrEstAng3[1] - arrEstAng3[2]
					-1.65 * sqrt(arrMtrxCorr3[4] - 2.* arrMtrxCorr3[5] + arrMtrxCorr3[8])) * 1000.; // функция разрешения 1 и 2 целей

	  parrBuff[ i * iNCols + 19] = (arrEstAng3[0] - arrEstAng3[1]
					-1.3 * sqrt(arrMtrxCorr3[0] - 2.* arrMtrxCorr3[1] + arrMtrxCorr3[4])) * 1000.;// функция разрешения 0 и 1 целей
	  parrBuff[ i * iNCols + 20] = (arrEstAng3[1] - arrEstAng3[2]
					-1.3 * sqrt(arrMtrxCorr3[4] - 2.* arrMtrxCorr3[5] + arrMtrxCorr3[8])) * 1000.; // функция разрешения 1 и 2 целей
	  // РСМ метод
	  parrBuff[ i * iNCols + 17] =  valEstRSM * 1000.;
	  parrBuff[ i * iNCols + 18] = valXi2 ;

	 }
	 free(pcmpSZv);
	 free(pcmpS);

	 wchar_t *wcharrFileNames = new wchar_t[iNCols * 30];
	 memset(wcharrFileNames, 0,iNCols * 30* sizeof (wchar_t));
	 wcscpy(&wcharrFileNames[0],L"n");
	 wcscpy(&wcharrFileNames[30],L"plnUMTrg0True");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnUMTrg1True");
	 wcscpy(&wcharrFileNames[30 * 3],L"plnUMTrg2True");

	 wcscpy(&wcharrFileNames[30 * 4],L"plnUMEstTrg0_Hyp2");
	 wcscpy(&wcharrFileNames[30 * 5],L"plnUMEstTrg1_Hyp2");
	 wcscpy(&wcharrFileNames[30 * 6],L"plnSKZTrg0_Hyp2");
	 wcscpy(&wcharrFileNames[30 * 7],L"plnSKZTrg1_Hyp2");

	 wcscpy(&wcharrFileNames[30 * 8],L"plnFuncRazr_Hyp2");


	 wcscpy(&wcharrFileNames[30 * 9],L"plnUMEstTrg0_Hyp3");
	 wcscpy(&wcharrFileNames[30 * 10],L"plnUMEstTrg1_Hyp3");
	 wcscpy(&wcharrFileNames[30 * 11],L"plnUMEstTrg2_Hyp3");
	 wcscpy(&wcharrFileNames[30 * 12],L"plnSKZTrg0_Hyp3");
	 wcscpy(&wcharrFileNames[30 * 13],L"plnSKZTrg1_Hyp3");
	 wcscpy(&wcharrFileNames[30 * 14],L"plnSKZTrg2_Hyp3");

	 wcscpy(&wcharrFileNames[30 * 15],L"plnFuncRazr_0_1_Hyp3_P005");
	 wcscpy(&wcharrFileNames[30 * 16],L"plnFuncRazr_1_2_Hyp3_P005");
	 wcscpy(&wcharrFileNames[30 * 19],L"plnFuncRazr_0_1_Hyp3_P01");
	 wcscpy(&wcharrFileNames[30 * 20],L"plnFuncRazr_1_2_Hyp3_P01");

	 wcscpy(&wcharrFileNames[30 * 17],L"EstRSM");
	 wcscpy(&wcharrFileNames[30 * 18],L"Xi2_RSM");



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
	 ,-15.,15., 2.) ;

 delete []wcharrFileNames;
	///
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr1.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-1.5, 1.5
	 ,-1.5, 1.5, 0.2) ;

}


// пОстрогнеие  графика иллюстрирующих изменение
// СКЗ ошибок измерения в зависимости от СКЗ шума в суммарной диаграмме
 void TFar::makeGraphs_FuncRazr_from_AngTarg2_Hyp3(wchar_t *wchFoldName, double *arrUMtrg, TComp *cmparrKTarg )

{
  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double arrEstAng3[3] = {0.}, arrMtrxCorr3[9] ={0.}, arrUMtrg_Cur[3] ={0.}, valXiSquare = 0.;
  memcpy(arrUMtrg_Cur, arrUMtrg, 3 * sizeof(double));
  TComp cmparrZ3[3],cmparrEstK3[3];
  int iNCols = 5;
  int nPoints = 1000;
   double * parrBuff = (double *) malloc(iNCols * nPoints * sizeof(double));
   ImitateMeasureArray(3, arrUMtrg,cmparrKTarg, pcmpS, pcmpSZv);


  wchar_t wcharrFileNames[5 * 30] = {0};
	 wcscpy(&wcharrFileNames[0],L"AngTarg2");
	 wcscpy(&wcharrFileNames[30],L"plnFuncRazr1");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnFuncRazr2");
	 wcscpy(&wcharrFileNames[30 * 3],L"plnFuncRazrCommon");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnFuncRazrHup2");
	 double valANgTarg2_0 = 0.001,valANgTarg2_cur = 0.001;

	 for (int i = 0; i < nPoints; i++)
	 {
	  memcpy(arrUMtrg_Cur, arrUMtrg, 3 * sizeof(double));
	  valANgTarg2_cur =  valANgTarg2_0 + ((double)i )* (0.01 - valANgTarg2_0)/ ((double)nPoints);
	  arrUMtrg_Cur[1] = valANgTarg2_cur;
	 // arrUMtrg_Cur[2] = -arrUMtrg_Cur[1];
	  ImitateMeasureArray(3, arrUMtrg_Cur,cmparrKTarg, pcmpS, pcmpSZv);
	  fncPolyn3(pcmpS,  arrEstAng3, cmparrEstK3, arrMtrxCorr3 , cmparrZ3, &valXiSquare);
	  parrBuff[ i* iNCols] = valANgTarg2_cur * 1000.;
	  parrBuff[ i * iNCols + 1] = (arrEstAng3[0] - arrEstAng3[1]
					-1.65 * sqrt(arrMtrxCorr3[0] - 2.* arrMtrxCorr3[1] + arrMtrxCorr3[4])) * 1000.;// функция разрешения 0 и 1 целей
	  parrBuff[ i * iNCols + 2] = (arrEstAng3[1] - arrEstAng3[2]
					-1.65 * sqrt(arrMtrxCorr3[4] - 2.* arrMtrxCorr3[5] + arrMtrxCorr3[8])) * 1000.; // функция разрешения 1 и 2 целей


	   arrUMtrg_Cur[2] = -valANgTarg2_cur;
	   ImitateMeasureArray(3, arrUMtrg_Cur,cmparrKTarg, pcmpS, pcmpSZv);
	   fncPolyn3(pcmpS,  arrEstAng3, cmparrEstK3, arrMtrxCorr3 , cmparrZ3, &valXiSquare);
	   parrBuff[ i * iNCols + 3] = (arrEstAng3[1] - arrEstAng3[2]
					-1.65 * sqrt(arrMtrxCorr3[4] - 2.* arrMtrxCorr3[5] + arrMtrxCorr3[8])) * 1000.; // функция разрешения 1 и 2 целей

   // функция разрешеничя 2 целей
   ImitateMeasureArray(2, arrUMtrg_Cur,cmparrKTarg, pcmpS, pcmpSZv);
   fncPolyn2(pcmpS,  arrEstAng3, cmparrEstK3, arrMtrxCorr3 , cmparrZ3, &valXiSquare);
   parrBuff[ i * iNCols + 4] = (arrEstAng3[0] - arrEstAng3[1]
					-1.65 * sqrt(arrMtrxCorr3[0] - 2.* arrMtrxCorr3[1] + arrMtrxCorr3[3])) * 1000.; // функция разрешения 1 и 2 целей

	 }

	 wchar_t wchFoldName1[300] ={0.};
		wcscpy(wchFoldName1, wchFoldName);
			wcscat(wchFoldName1, L"\\");
	for (int i = 1; i < iNCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName1  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,iNCols// - к-во переменных о корорых накоплена информация в буфере
								  ,nPoints //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , 1.  //  масштаб по оси X
								  , 1.  // масштаб по оси Y
								   ) ;
	 }

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName1);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 12.
	 ,-11.,11., 0.5) ;

   free(pcmpSZv);
   free(pcmpS);
   free(parrBuff);
}

// пОстрогнеие группы графиков иллюстрирующих изменение
// СКЗ ошибок измерения в зависимости от СКЗ шума в суммарной диаграмме
 void TFar::makeGraph_SKZ_from_Sig(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp,  double  valCoefRange)
{
  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[4] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);
   double * parrBuff = (double *) malloc(6 * 2 * sizeof(double));
   ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

   fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	  , &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr );
  wchar_t wcharrFileNames[6 * 30] = {0};
	 wcscpy(&wcharrFileNames[0],L"Sig");
	 wcscpy(&wcharrFileNames[30],L"plnSkzUMTrg");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnSkzUMAntp");
	 wcscpy(&wcharrFileNames[30 * 3],L"pln_1_65_SkzDifferm");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnDiffer");
	 wcscpy(&wcharrFileNames[30 * 5],L"pln_1_28_SkzDifferm");

	 double scalex = 6.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
   parrBuff[0] = 0.;
   parrBuff[1] = 0.;
   parrBuff[2] = 0.;
   parrBuff[3] = 0.;
   parrBuff[4] = (alfUMTrg - alfUMAntp)* 1000.;
   parrBuff[5] = 0.;

   parrBuff[6] = valCoefRange;
   parrBuff[7] = sqrt(arrMtrxCorr[0]) * valCoefRange * 1000.;
   parrBuff[8] = sqrt(arrMtrxCorr[3]) * valCoefRange * 1000.;
   parrBuff[9] = 1.65* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valCoefRange* 1000.;
   parrBuff[10] = (alfUMTrg - alfUMAntp)* 1000.;
   parrBuff[11] = 1.28 * sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valCoefRange * 1000.;
	for (int i = 1; i < 6; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,6// - к-во переменных о корорых накоплена информация в буфере
								  ,2 //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , scalex  //  масштаб по оси X
								  , scaley  // масштаб по оси Y
								   ) ;
	 }

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., valCoefRange *scalex  * 1.2
	 ,0.,parrBuff[9] * 1.1, valCoefRange * scalex *0.05) ;

   free(pcmpSZv);
   free(pcmpS);
   free(parrBuff);
}


// пОстрогнеие группы графиков иллюстрирующих изменение
// СКЗ ошибок измерения в зависимости от отношение СКЗ шума в суммарной диаграмме
// в децибелах
 void TFar::makeGraph_SKZ_from_dBell(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,double  alfUMAntp,TComp cmpKAntp,  double  valCoefRange)
{
  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[4] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);
  const int iNUmPoints = 1000;
   double * parrBuff = (double *) malloc(6 * iNUmPoints * sizeof(double));
   ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

   fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	  , &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr );
  wchar_t wcharrFileNames[6 * 30] = {0};
	 wcscpy(&wcharrFileNames[0],L"Sig");
	 wcscpy(&wcharrFileNames[30],L"plnSkzUMTrg");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnSkzUMAntp");
	 wcscpy(&wcharrFileNames[30 * 3],L"pln2SkzDifferm");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnDiffer");
	 wcscpy(&wcharrFileNames[30 * 5],L"pln1_2SkzDifferm");

	 double scalex = 3.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
   parrBuff[0] = 0.;
   parrBuff[1] = 0.;
   parrBuff[2] = 0.;
   parrBuff[3] = 0.;
   parrBuff[4] = (alfUMTrg - alfUMAntp)* 1000.;
   parrBuff[5] = 0.;

   parrBuff[6] = valCoefRange;
   parrBuff[7] = sqrt(arrMtrxCorr[0]) * valCoefRange * 1000.;
   parrBuff[8] = sqrt(arrMtrxCorr[3]) * valCoefRange * 1000.;
   parrBuff[9] = 2.* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valCoefRange* 1000.;
   parrBuff[10] = (alfUMTrg - alfUMAntp)* 3000./ M_PI;
   parrBuff[11] = 1.2 * sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valCoefRange * 1000.;
   double step = 0.25;
	for (int i = 0; i < iNUmPoints; i++)
	{
	  double valNoiseSig = ((double)(i +1))* step;
	  parrBuff[i * 6] = 10. * log(cmpKTarg.modul()/valNoiseSig)/log(10.);
	  parrBuff[i * 6 + 1] = sqrt(arrMtrxCorr[0]) * valNoiseSig * 1000.;
	  parrBuff[i * 6 + 2] = sqrt(arrMtrxCorr[3]) * valNoiseSig * 1000.;
	  parrBuff[i * 6 + 3] =  2.* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valNoiseSig* 1000.;
	  parrBuff[i * 6 + 4] = (alfUMTrg - alfUMAntp)* 1000.;
	  parrBuff[i * 6 + 5] = 1.2 * sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) * valNoiseSig * 1000.;

	}
	for (int i = 1; i < 6; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  ,6// - к-во переменных о корорых накоплена информация в буфере
								  ,iNUmPoints //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , scalex  //  масштаб по оси X
								  , scaley  // масштаб по оси Y
								   ) ;
	 }

	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., 120.
	 ,0.,70., 120 *0.05) ;
   free(pcmpSZv);
   free(pcmpS);
   free(parrBuff);
}

// Построение группы графиков иллюстрирующих изменение СКЗ ошибок измерения
// в зависимости  от положения  АНП
// угол АНП меняеется от -80 тд до  alfUMTrg
 void TFar::makeGraph_SKZ_from_AngDiffer(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
			,TComp cmpKAntp)
{
  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[4] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);

  const int NCols = 11;
  const int NIsp = 1000;
  double * parrBuff = (double *) malloc(NCols * NIsp * sizeof(double));
  double alfUMAntp = 0.;
  double valStep =  (alfUMTrg + 83. * M_PI/ 3000.)/  ((double)(NIsp+2));
  for (int i = 0; i < NIsp; i++)
	 {
		alfUMAntp =  -83. * M_PI/ 3000. + ((double)(i + 1)) * valStep;
		ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);


		fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	  , &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr );

			double valXi2  =0., valEstRSM =0.,  valRSMDisp =0.;

	  fncMeasureRSMProcessing(pcmpS, &valXi2, &valEstRSM, &valRSMDisp  );


	   const double ValAngParam = 0.01;
	   double valEstWhite =0., valWhiteDisp  =0.;
	  fncMeasureWhiteProcessing(ValAngParam,pcmpS
 , &valXi2, &valEstWhite,&valWhiteDisp  );

	  parrBuff[ i * NCols] = (alfUMTrg - alfUMAntp)* 1000.;
	  parrBuff[ i * NCols + 1] =  1.65* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3])  * 1000.;
	  parrBuff[ i * NCols + 2] =  alfUMTrg * 1000.;
	  parrBuff[ i * NCols + 3] =  alfUMAntp * 1000.;

	  parrBuff[ i * NCols + 4] = sqrt(arrMtrxCorr[0]) * 1000.;
	  parrBuff[ i * NCols + 5] = sqrt(arrMtrxCorr[3]) * 1000.;

	  parrBuff[ i * NCols + 6] = valEstRSM * 1000.;
	  parrBuff[ i * NCols + 7] = parrBuff[ i * NCols];
	  parrBuff[ i * NCols + 8] = 1.28* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3])  * 1000.;
	  parrBuff[ i * NCols + 9] =  valXi2;
	  parrBuff[ i * NCols + 10] =  valEstWhite* 1000.;


	 }
	 free(pcmpSZv);
	 free(pcmpS);

	 wchar_t wcharrFileNames[11 * 30] = {0};
	 wcscpy(&wcharrFileNames[0],L"AngDiffer"); //  угловое расстояние НЛЦ - АНТП
	 wcscpy(&wcharrFileNames[30],L"1_65_SigDiffer");
	 wcscpy(&wcharrFileNames[30 * 2],L"plnUMTargTrue");
	 wcscpy(&wcharrFileNames[30 * 3],L"plnUMAntpTrue");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnSigTarg");
	 wcscpy(&wcharrFileNames[30 * 5],L"plnSigAntp");
	 wcscpy(&wcharrFileNames[30 * 6],L"plnUM_RSMIzm");
	 wcscpy(&wcharrFileNames[30 * 7],L"AngDiffer");

	 wcscpy(&wcharrFileNames[30 * 8],L"1_28_SigDiffer");
	 wcscpy(&wcharrFileNames[30 * 9],L"Xi2");
	 wcscpy(&wcharrFileNames[30 * 10],L"plnUM_WhiteIzm");

	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i <  NCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , NCols // - к-во переменных о корорых накоплена информация в буфере
								  ,NIsp //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , scalex  //  масштаб по оси X
								  , scaley  // масштаб по оси Y
								   ) ;
	 }



	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,0., parrBuff[0] * 1.1
	 ,-parrBuff[0] * 1.1,parrBuff[0] * 1.1, parrBuff[0] *0.05) ;
	 free(parrBuff);
}

// пОстрогнеие группы графиков иллюстрирующих изменение СКЗ ошибок измерения
// в зависимости  от фазового угла сигнала антипода
//фаз  угол АНП меняеется от -пи тд до  пи
//строятся графики:
// ско цели
// ско антп
// разность углов (горизонтальная линия)
// коефф корреляции
// левая часть неравенства разрешения
 void TFar::makeGraph_SKZ_from_AntpPhAng(wchar_t *wchFoldName1, double alfUMTrg, TComp  cmpKTarg
	   ,double alfUMAntp,TComp cmpKAntp)
{
  TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
  TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
  double valEstAngTarg = 0., valEstAngAntp =0., arrMtrxCorr[4] ={0.};
  TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.),  cmpKAntpCur, cmpKAntpMod( cmpKAntp.modul(),0.);;

  const int NCols = 6;
  const int NIsp = 360;
  double * parrBuff = (double *) malloc(NCols * NIsp * sizeof(double));


  for (int i = 0; i < NIsp; i++)
	 {

		cmpKAntpCur = exp_(-M_PI +  ((double)i)/ 180. * M_PI) *cmpKAntpMod;
		ImitateMeasureArray( alfUMTrg,  cmpKTarg
			,  alfUMAntp, cmpKAntpCur, pcmpS, pcmpSZv);

		fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	  , &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr );



	  parrBuff[ i * NCols] = (-M_PI +  ((double)i)/ 180. * M_PI)/2.;
	  parrBuff[ i * NCols + 1] =  2.* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3])  * 1000.;
	  parrBuff[ i * NCols + 2] = fabs( alfUMTrg -alfUMAntp)* 1000.;
	  parrBuff[ i * NCols + 3] = sqrt(arrMtrxCorr[0]) * 1000.;
	  parrBuff[ i * NCols + 4] = sqrt(arrMtrxCorr[3]) * 1000.;

	  parrBuff[ i * NCols + 5] =arrMtrxCorr[1]/sqrt(arrMtrxCorr[0])/sqrt(arrMtrxCorr[3]);

	 }
	 free(pcmpSZv);
	 free(pcmpS);

	 wchar_t wcharrFileNames[10 * 30] = {0};
	 wcscpy(&wcharrFileNames[0],L"AntpPhAng"); //  угловое расстояние НЛЦ - АНТП
	 wcscpy(&wcharrFileNames[30],L"2SigDiffer");
	 wcscpy(&wcharrFileNames[30 * 2],L"AngDiffer");
	 wcscpy(&wcharrFileNames[30 * 3],L"plnSigTarg");
	 wcscpy(&wcharrFileNames[30 * 4],L"plnSigAntp");
	 wcscpy(&wcharrFileNames[30 * 5],L"plnCoefCor");


	 double scalex = 1.;
	 double scaley = 1.;
		wchar_t wchFoldName[300] ={0};
	wcscpy(  wchFoldName,  wchFoldName1);
	wcscat(wchFoldName, L"\\");
	 for (int i = 1; i <  NCols; i++)
	 {
	  TYrWriteShapeFile::WriteOneReport(wchFoldName  // путь к папке
								  , parrBuff // массив с информацией - матрица nBuffRows x nBuffCols
								  , NCols // - к-во переменных о корорых накоплена информация в буфере
								  ,NIsp //  - к-во точек
								  ,wcharrFileNames //матрица с именаими переменных - матрица nBuffCols x lenName
								  ,30 // максимальная длина имени переменной
								  ,0  // номер переменной по оси X
								  ,i   // номер переменной по оси Y
								  , scalex  //  масштаб по оси X
								  , scaley  // масштаб по оси Y
								   ) ;
	 }



	wchar_t wchAxesFileName[300] ={0};
	wcscpy(  wchAxesFileName,  wchFoldName);
	wcscat(wchAxesFileName, L"AxesArr.shp");
	TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName,-M_PI/2. * 1.1, M_PI/2. * 1.1
	 ,-1.2, 1.2, 0.1) ;

	 free(parrBuff);
}


//METHOD WHITE
 void TFar::fncMeasureWhiteProcessing(const double ValAngParam, TComp *pcmpSZv
 , double *valXi2, double *valEstWhite, double *valWhiteDisp  )
{
 // формирование суммарной диаграммы
 TComp cmpSSumZv (0.,0.);
 for (int i = 0; i < m_N; i++)
 {
   cmpSSumZv += pcmpSZv[i] ;
 }
 ///
 // измерения квадратичной диграммы
 TComp cmpSSquareZv = createSquareMeasure(ValAngParam, pcmpSZv);
 ///

 // отношение
  TComp cmpYZv = cmpSSquareZv/ cmpSSumZv;
  double valYZv = cmpYZv.m_Re;
  ///

  double valx0 = ValAngParam /2., valToler = 0.00001, val_dx_po_dY = 0.;
 findRootWhiteEquation(&valx0, ValAngParam,valYZv,valToler, &val_dx_po_dY );
 *valEstWhite = valx0;
 double valSig2YZv = calcSig2Y_MethodWhite(ValAngParam, pcmpSZv);
 *valXi2 = 0.;
 *valWhiteDisp = valSig2YZv * val_dx_po_dY * val_dx_po_dY;

 
}

// вычисление дисперсии отношения диаграмм для метода Уайта
double TFar::calcSig2Y_MethodWhite(const double ValAngPar, TComp *pcmpSZv)
{
   double valrez = 0.;

   TDiagr diagrPlus(*this, ValAngPar);
   TDiagr diagrMinus(*this, -ValAngPar);
   TDiagr diagrSquare(*this, -ValAngPar);
	for (int i = 0; i < m_N; i++)
	{
	diagrSquare.mparrCmpElCurrent[i] =  diagrPlus.mparrCmpElCurrent[i]  + diagrMinus.mparrCmpElCurrent[i] ;

	}
	// формирование суммарной диаграммы
 TComp cmpSSumZv (0.,0.);
 for (int i = 0; i < m_N; i++)
 {
   cmpSSumZv += pcmpSZv[i] ;
 }
 // формирование измерения квадратичной диаграммы
	TComp cmpYSqZv = createSquareMeasure( ValAngPar, pcmpSZv);
	///

   for (int i = 0; i < m_N; i++)
   {
	 TComp cmp_dYZv_po_dS = (diagrSquare.mparrCmpElCurrent[i] *  cmpSSumZv - cmpYSqZv )/ ( cmpSSumZv* cmpSSumZv);
	 valrez +=  cmp_dYZv_po_dS.modul() * cmp_dYZv_po_dS.modul() * mparrDisp[i] ;
   }
  return valrez;
}

 TComp TFar::createSquareMeasure(const double ValAngPar, TComp *pcmpSZv)
{


	TComp *pcmpSZvCur =  (TComp *)malloc(m_N * sizeof(TComp ));
	for (int i = 0; i < 4; i++)
	{
	  pcmpSZvCur[i] = pcmpSZv[i] * TComp(4.,0.);
	}
	TDiagr diagrPlus(*this, ValAngPar);
	TDiagr diagrMinus(*this, -ValAngPar);
	TDiagr diagrSquare(*this, -ValAngPar);
	//TComp cmp2(2.,0.);
	for (int i =0; i < m_N; i++)
	{
	  diagrSquare.mparrCmpElCurrent[i] =  diagrPlus.mparrCmpElCurrent[i] + diagrMinus.mparrCmpElCurrent[i];
	 // diagrSquare.mparrCmpElCurrent[i] = diagrSquare.mparrCmpElCurrent[i]/ cmp2;
	}

	TComp cmpRez (0.,0.);
	for (int i =0; i < m_N; i++)
	{
	  cmpRez +=  diagrSquare.mparrCmpElCurrent[i] * pcmpSZvCur[i] ;
    }
	free(pcmpSZvCur);

	return cmpRez;
}

// решение уравнения уайта методоим ньютонп
bool TFar::findRootWhiteEquation(double *valx0,const double ValAngPar
	 ,double valYZv, const double valToler, double *val_dx_po_dY )
{ bool bsucceed = false;
	TDiagr diagrPlus(*this, ValAngPar);
	TDiagr diagrMinus(*this, -ValAngPar);
	TDiagr diagrSum(*this, 0.);
	double valFGr =0., val_dFgr = 0.;

	TComp cmpYZv(valYZv, 0.);
	for (int i =0; i < 100; i++)
	{
	 
	  valFGr =  diagrPlus.fncF(*valx0 ).m_Re+ diagrMinus.fncF(*valx0 ).m_Re
	   - valYZv * diagrSum.fncF(*valx0).m_Re;
	   val_dFgr = diagrPlus.dF_po_dTet(*valx0 ).m_Re + diagrMinus.dF_po_dTet(*valx0 ).m_Re
	   - valYZv * diagrSum.dF_po_dTet(*valx0).m_Re;
	   double valDel = - valFGr/val_dFgr;
	   *valx0  += valDel ;
	   if (fabs(valDel ) <valToler)
	   {
	   bsucceed = true;
		break;
	   }

	}
	if ( bsucceed)
	{
	   double	val_dF0_po_dx = diagrPlus.d2F_po_dTet2(*valx0 ).m_Re+ diagrMinus.d2F_po_dTet2(*valx0 ).m_Re
	   - valYZv * diagrSum.d2F_po_dTet2(*valx0).m_Re;
	  double val_dF0_po_dy = - diagrSum.fncF(*valx0).m_Re;
	   *val_dx_po_dY =  val_dF0_po_dy / val_dF0_po_dx;
    }
  return bsucceed;
}



// Вычисление теоретического значения СКЗ ошибки измерения угла цели и разрешения цели и антипода
// квадратичным методом
// INPUT:
// alfUMTrg  - ум цели
// alfUMAntp - ум АНТП
// cmpKTarg  - компл сишгнал цели
// cmpKAntp - компл сишгнал АНТП
// OUTPUT:
// pbRazreshenie - признак разрешения
// pvalSigTarg  - скз ошибки измерения ум цели
 void TFar::calcTheoretical_LATarg_SKZ(bool *pbRazreshenie,double *pvalSigTarg
	, double alfUMTrg,  double alfUMAntp, TComp  cmpKTarg, TComp cmpKAntp)
{
	TComp *pcmpSZv = (TComp *)malloc( m_N * sizeof( TComp));
	TComp *pcmpS = (TComp *)malloc( m_N * sizeof( TComp));
	double  arrMtrxCorr[4] ={0.};
	TComp cmpEstKTarg (0.,0.),cmpEstKAntp (0.,0.);

	ImitateMeasureArray( alfUMTrg,  cmpKTarg
	,  alfUMAntp, cmpKAntp, pcmpS, pcmpSZv);

   //	fncEstimateMsd(pcmpS,  &valEstAngTarg, &valEstAngAntp
	//, &cmpEstKTarg , &cmpEstKAntp , arrMtrxCorr );
	double arrEstAng[2] = {0.};
	TComp cmparrK[2];
	TComp cmparrZ[2];
	fncPolyn2(pcmpS,  arrEstAng
	  , cmparrK, arrMtrxCorr , cmparrZ);


	if (1.65* sqrt(arrMtrxCorr[0]- 2. * arrMtrxCorr[1] + arrMtrxCorr[3]) <=  fabs(alfUMTrg - alfUMAntp))
	{
	*pbRazreshenie = true;
	}
	else
	{
	*pbRazreshenie = false;
	}

	*pvalSigTarg  = (arrMtrxCorr[0] < arrMtrxCorr[3])?arrMtrxCorr[0]:arrMtrxCorr[3];
	*pvalSigTarg = sqrt(*pvalSigTarg) ;
	free(pcmpSZv);
	free(pcmpS);
}

// Вычисление СКЗ ошибки применения РСМ когда принимаемый сигнал является суммой сигнала цели
// и сигнала едиственного антипода
// Сначала рассчитывается дисперсия флуктуационной ошибки
// затем систематическая. Складываются под корнем
// INPUT:
// VAlUMTarg ,  VAlUMAntp   - УМ цели и антипода, рад
// CMpKTarg,  CMpKAntp - коэффиц отражения цели и антипода
// OUTPUT:
// возвращает СКЗ ошибки оценивания угла цели
//  calc_Guaranted_SKZ_RSM_For_TwoTargs
double  TFar::calc_Guaranted_SKZ_RSM_For_TwoTargs ( double valUMTarg, double valUMAntp
	,  TComp  cmpKTarg,  TComp cmpKAntp)
{
	TComp pcmpS[4], pcmpSZv[4], cmpAntp1(0.,0.);
	ImitateMeasureArray(valUMTarg, cmpKTarg	,valUMAntp, cmpAntp1, pcmpS, pcmpSZv);
	double valXi2  =0., valEstRSM  =0., valRSMDisp =0.;
	fncMeasureRSMProcessing(pcmpS, &valXi2, &valEstRSM, &valRSMDisp  );
	if (fabs(valUMTarg - valUMAntp) > 2. * mLambda / ( m_D * ((double)m_N)) )
	{
	return sqrt(valRSMDisp);
	}

	TComp cmpCentre ;
	double valRadius =0., valArgFi =0.;
	double valSysErr = calcGuarantSystError_RSM_For_2Targs ( valUMTarg,  valUMAntp
	, cmpKTarg.modul()  , cmpKAntp.modul() , &cmpCentre, &valRadius, &valArgFi) ;
	double valSigE = sqrt(valRSMDisp + valSysErr * valSysErr) ;
	return valSigE;
}

// Вычисление СКЗ ошибки применения РСМ когда принимаемый сигнал является суммой сигнала цели
// и сигнала едиственного антипода
// Сначала рассчитывается дисперсия флуктуационной ошибки
// затем средняя систематическая. Складываются под корнем
// INPUT:
// VAlUMTarg ,  VAlUMAntp   - УМ цели и антипода, рад
// CMpKTarg,  CMpKAntp - коэффиц отражения цели и антипода
// OUTPUT:
// возвращает СКЗ ошибки оценивания угла цели
//  calc_Guaranted_SKZ_RSM_For_TwoTargs
double  TFar::calc_Mean_SKZ_RSM_For_TwoTargs ( double valUMTarg, double valUMAntp
	,  TComp  cmpKTarg,  TComp cmpKAntp)
{
	TComp pcmpS[4], pcmpSZv[4], cmpAntp1(0.,0.);
	ImitateMeasureArray(valUMTarg, cmpKTarg	,valUMAntp, cmpAntp1, pcmpS, pcmpSZv);
	double valXi2  =0., valEstRSM  =0., valRSMDisp =0.;
	fncMeasureRSMProcessing(pcmpS, &valXi2, &valEstRSM, &valRSMDisp  );

	TComp cmpCentre ;
	double valRadius =0., valArgFi =0.;
	double valSysErr = calcMeanSystError_RSM_For_2Targs ( valUMTarg,  valUMAntp
	, cmpKTarg.modul()  , cmpKAntp.modul() , &cmpCentre, &valRadius, &valArgFi) ;
	double valSigE = sqrt(valRSMDisp + valSysErr * valSysErr) ;
	return valSigE;
}

// вычисление дисперсии ошибки измерения угла РСМ
// INPUT:
// VAlNWaveCur - отношение шум/сигнал
// VAlNAppert - апертура
// VAlLamb - длина волны
// VAlTetta - угол                                ImitateMeasureArray
double TFar::calcTheorDisp_RSM(const double VAlNWaveCur, const double VAlNAppert
, const double VAlLamb, const double VAlTetta)
{
	//double valMu = VAlNAppert * sin( VAlTetta) * M_PI/2./ VAlLamb;
	//double valTemp = cos (valMu/2.);
	//double valTemp1 = VAlLamb / M_PI / VAlNAppert/ (valTemp * valTemp) * VAlNWaveCur;
	//return 2. * valTemp1 * valTemp1 ;
	// по новому
	// TSincDgr SincDgr (VAlNAppert);
	// double valFDiagr = SincDgr.fncDiagrFromRad(VAlTetta, VAlLamb ) ;
	// double valTemp2 = VAlLamb / M_PI / VAlNAppert/ valFDiagr * VAlNWaveCur;
	// double b = 4. * valTemp2 * valTemp2;
	// return b;


	 double vald = VAlNAppert /2.;
	 double temp =  VAlLamb/ M_PI/ vald;
	 double fi = 2. * M_PI * vald / VAlLamb * VAlTetta ;
	 double valFDiagr = fncDiagrSinx_div_x(fi);
	 double sigSq = VAlNWaveCur *VAlNWaveCur * cos(fi/2.)* cos(fi/2.)/valFDiagr/valFDiagr * temp* temp/2.;

	 return sigSq;

}

// нахождение максимальной сиситематической ошибки измерения угла цели
// в присуьтсвии помехи - сигнала от второй цели
// INPUT:
// VAlUMTarg, VAlUMAntp - углы цели и антипода
// VAlTargAmp, VAlAntpAmp_Ro - амплитуды сигналов цели и антипода
// OUTPUT:
// pcmpCircleCentre - центр образа единичного круга
// *pRadius - его радиус
// *pvalDeltaFi - фазовый угол сигнала антипода, реализующий максимум сиситематическойц ошибки
// Возвращает максимальную ошибку в рад
double  TFar::calcGuarantSystError_RSM_For_2Targs (const double VAlUMTarg, const double VAlUMAntp
	, const double VAlTargAmp, const double VAlAntpAmp_Ro, TComp *pcmpCentre, double *pvalRadius, double *pvalDeltaFi)
{
	const double VAlAppert = m_D * ((double)m_N);
	TSincDgr SincDgr (VAlAppert /2. );
	double valFTarg =  SincDgr.fncDiagrFromRad(VAlUMTarg,  mLambda ) ;
	double valFAntp =  SincDgr.fncDiagrFromRad(VAlUMAntp,  mLambda ) ;
	double valMuTarg = VAlUMTarg * VAlAppert* M_PI/ mLambda;
	//SincDgr.transformAngToGeneralizedAng (VAlUMTarg,  mLambda );
	//double valMuAntp = SincDgr.transformAngToGeneralizedAng (VAlUMAntp,  mLambda );
	double valMuAntp = VAlUMAntp * VAlAppert* M_PI/ mLambda;

	TComp cmpZTarg = exp_(valMuTarg);
	TComp cmpZAntp = exp_(valMuAntp);
	TComp cmp_1(1.,0.);

	TComp cmp0 =  cmp_1 -  cmpZTarg;
	TComp cmp1 =  cmp_1 -  cmpZAntp;

	TComp cmp4 =  cmp_1 + cmpZTarg;
	TComp cmp5 =  cmp_1 + cmpZAntp;

	TComp cmp2(VAlTargAmp * valFTarg, 0.);
	TComp cmp3(VAlAntpAmp_Ro * valFAntp, 0.);

	TComp cmpb =  cmp2 * cmp0;
	TComp cmpa =  cmp3 * cmp1;


	TComp cmpd =  cmp2 * cmp4;
	TComp cmpc =  cmp3 * cmp5;

	TComp cmpArg0(1.,0.), cmpArg1(0.,1.),cmpArg2(-1.,0.);
	TComp cmpImage0(0.,0.),cmpImage1(0.,0.),cmpImage2(0.,0.);

	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg0, &cmpImage0);
	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg1, &cmpImage1);
	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg2, &cmpImage2);

	findCircleParams (cmpImage0,cmpImage1,cmpImage2, pcmpCentre, pvalRadius);

	TComp cmpMax0 = *pcmpCentre;
	cmpMax0.m_Im += (*pvalRadius );

	TComp cmpMax1 = *pcmpCentre;
	cmpMax1.m_Im -= (*pvalRadius );
	double val0 = -2. * atan(cmpMax0.m_Im);
	double val1 = -2. * atan(cmpMax1.m_Im) ;
	double valDelMu0 = fabs(val0 - valMuTarg);
	double valDelMu1 = fabs(val1 - valMuTarg);
	if (valDelMu0 < valDelMu1)
	{
		valDelMu0 =  valDelMu1;

		cmpMax0 = cmpMax1;
	}
	TComp cmpb1 =  cmpb * TComp(-1.,0.) ;
	TComp cmpc1 =  cmpc * TComp(-1.,0.) ;
	 TComp cmpArgMaxZ;
	fncLinFrac( cmpd,  cmpb1,  cmpc1,  cmpa, cmpMax0, &cmpArgMaxZ);

	*pvalDeltaFi = cmpArgMaxZ.phase() ;

	

	return valDelMu0 * mLambda /M_PI/ VAlAppert;
}

// нахождение средней сиситематической ошибки измерения угла цели
// в присуьтсвии помехи - сигнала от второй цели
// INPUT:
// VAlUMTarg, VAlUMAntp - углы цели и антипода
// VAlTargAmp, VAlAntpAmp_Ro - амплитуды сигналов цели и антипода
// OUTPUT:
// pcmpCircleCentre - центр образа единичного круга
// *pRadius - его радиус
// *pvalDeltaFi - фазовый угол сигнала антипода, реализующий максимум сиситематическойц ошибки
// Возвращает максимальную ошибку в рад
double  TFar::calcMeanSystError_RSM_For_2Targs (const double VAlUMTarg, const double VAlUMAntp
	, const double VAlTargAmp, const double VAlAntpAmp_Ro, TComp *pcmpCentre, double *pvalRadius, double *pvalDeltaFi)
{
	const double VAlAppert = m_D * ((double)m_N);
	TSincDgr SincDgr (VAlAppert /2. );
	double valFTarg =  SincDgr.fncDiagrFromRad(VAlUMTarg,  mLambda ) ;
	double valFAntp =  SincDgr.fncDiagrFromRad(VAlUMAntp,  mLambda ) ;
	double valMuTarg = VAlUMTarg * VAlAppert* M_PI/ mLambda;
	//SincDgr.transformAngToGeneralizedAng (VAlUMTarg,  mLambda );
	//double valMuAntp = SincDgr.transformAngToGeneralizedAng (VAlUMAntp,  mLambda );
	double valMuAntp = VAlUMAntp * VAlAppert* M_PI/ mLambda;

	TComp cmpZTarg = exp_(valMuTarg);
	TComp cmpZAntp = exp_(valMuAntp);
	TComp cmp_1(1.,0.);

	TComp cmp0 =  cmp_1 -  cmpZTarg;
	TComp cmp1 =  cmp_1 -  cmpZAntp;

	TComp cmp4 =  cmp_1 + cmpZTarg;
	TComp cmp5 =  cmp_1 + cmpZAntp;

	TComp cmp2(VAlTargAmp * valFTarg, 0.);
	TComp cmp3(VAlAntpAmp_Ro * valFAntp, 0.);

	TComp cmpb =  cmp2 * cmp0;
	TComp cmpa =  cmp3 * cmp1;


	TComp cmpd =  cmp2 * cmp4;
	TComp cmpc =  cmp3 * cmp5;

	TComp cmpArg0(1.,0.), cmpArg1(0.,1.),cmpArg2(-1.,0.);
	TComp cmpImage0(0.,0.),cmpImage1(0.,0.),cmpImage2(0.,0.);

	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg0, &cmpImage0);
	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg1, &cmpImage1);
	fncLinFrac( cmpa,  cmpb,  cmpc,  cmpd, cmpArg2, &cmpImage2);

	findCircleParams (cmpImage0,cmpImage1,cmpImage2, pcmpCentre, pvalRadius);
  /*
	TComp cmpMax0 = *pcmpCentre;
	cmpMax0.m_Re += (*pvalRadius );


	double val0 = -2. * atan(cmpMax0.m_Im);

	double valDelMu0 = fabs(val0 - valMuTarg);


	TComp cmpb1 =  cmpb * TComp(-1.,0.) ;
	TComp cmpc1 =  cmpc * TComp(-1.,0.) ;
	 TComp cmpArgMaxZ;
	fncLinFrac( cmpd,  cmpb1,  cmpc1,  cmpa, cmpMax0, &cmpArgMaxZ);

	*pvalDeltaFi = cmpArgMaxZ.phase() ;
	*/

	 double temp =	-2. * atan((*pcmpCentre).m_Im )* mLambda/M_PI/ VAlAppert - VAlUMTarg;
	return temp;
}

double max_d(double a, double b)
{
	return (a > b)?a:b;
}

double calcDisp_RSM(const double VAlNWaveCur, const double VAlNAppert
, const double VAlLamb, const double VAlTetta)
{
	//double valMu = VAlNAppert * sin( VAlTetta) * M_PI/2./ VAlLamb;
	//double valTemp = cos (valMu/2.);
	//double valTemp1 = VAlLamb / M_PI / VAlNAppert/ (valTemp * valTemp) * VAlNWaveCur;
	//return 2. * valTemp1 * valTemp1 ;
	// по новому
	 TSincDgr SincDgr (VAlNAppert);
	 double valFDiagr = SincDgr.fncDiagrFromRad(VAlTetta, VAlLamb ) ;
	 double valTemp2 = VAlLamb / M_PI / VAlNAppert/ valFDiagr * VAlNWaveCur;
	 double b = 4. * valTemp2 * valTemp2;
	 return b;
}

//----------------------------------------
// приблтженный расчет ширины диаграммы через sinx/x
 double TFar::findDiagrWidthApprox()
{
 double vala = 2. * M_PI / mLambda * m_D * ((double)m_N)/ 2.;
 double tet1 = TET0707 / vala;
 return tet1;
}

 // приблтженный расчет  диаграммы через sinx/x
double  TFar::fncFFarApprox (const double valTetta)
{
	double valWidth = findDiagrWidthApprox();
	return fncDiagrSimple(valWidth, valTetta );

}
#pragma package(smart_init)

