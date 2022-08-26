//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FormVessel_.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"


#include <math.h>
#include <dir.h>


#include "YrString.h"
#include <stdio.h>
#include <stdlib.h>
#include "AM_2D.h"
#include "Far_2D.h"
#include "Comp.h"
#include "Faceta.h"
#include "Far.h"
#include "Diagrams.h"
#include "Gauss.h"

#include "MatrixProccess.h"
#include "SincDgr.h"


#include "Target.h"
#include "YrWrite.h"

#include "Sins.h"
#include "DriverMech.h"

#include "Measurer.h"
#include "YrWriteShapeFile.h"

#include "InitTargData.h"
#include "CalcCorMatrx.h"

#include "Environment.h"
#include "Equations.h"
#include "MyShellTraj.h"

#include "SincDgr.h"
#include "URPointZ.h"
#include "ShellBody.h"
#include "Bius.h"
#include "Line.h"



extern const double VAL_C;
TForm5 *Form5;
//---------------------------------------------------------------------------
__fastcall TForm5::TForm5(TComponent* Owner)
	: TForm(Owner)
{
  for (int i = 0; i < MAX_QUANT_AM; i++)
	{
		mbarrAM[i] = true;
	}
	StringGrid1->Cells[0][1] = L"ФАР";
	StringGrid1->Cells[0][2] = L"АУ";
	StringGrid1->Cells[1][0]  = L"X";
	StringGrid1->Cells[2][0]  = L"Y";
	StringGrid1->Cells[3][0]  = L"Z";
	create5P10();
}
//---------------------------------------------------------------------------

void __fastcall TForm5::fncInputData()
{
	String s_22 = LabeledEdit2->Text;
	if (wcschr(s_22.w_str(), L','))
	{

		// длина волны
		mLambda = VAL_C / StrTo_Dbl_(LabeledEdit15->Text) * 100. / 1000000000.;
		LabeledEdit13->Text = mLambda;

		// расстояние между излучателями
		mdEmitCol = mLambda * 0.55;
		mdEmitRow = mLambda * 0.55;
		LabeledEdit12->Text = mdEmitCol;

		// ЭТАЛОННЫЙ СИГНАЛ ДЛЯ РАСЧЕТА АМПЛИТУДЫ
		// амплитуда
		mEtalonAmp = StrTo_Dbl_(LabeledEdit2->Text);
		// дальность
		mEtalonDist = StrTo_Dbl_(LabeledEdit3->Text) * 1000.;
		// ЭПР
		mEtalonAPR = StrTo_Dbl_(LabeledEdit4->Text);
		// СКО внутр шума суммарной диаграммы 5П10
		mNoiseSKZ_5P10 = StrTo_Dbl_(LabeledEdit16->Text);
		// СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
		mEtalonSigAmplFact_5P10 = StrTo_Dbl_(LabeledEdit18->Text);

		// мощность на передачу
		mEtalonPowerPrd = StrTo_Dbl_(LabeledEdit36->Text);
		// КУ на передачу
		mEtalonKYPrd = StrTo_Dbl_(LabeledEdit38->Text);
		// мощность на прием
		mEtalonKYPriem = StrTo_Dbl_(LabeledEdit37->Text);


		mEtalonSign = TEtalonSign(mEtalonAmp, mEtalonDist, mEtalonAPR,
			mNoiseSKZ_5P10, mEtalonSigAmplFact_5P10, mEtalonPowerPrd, mEtalonKYPrd
          , mEtalonKYPriem);
		///

		//  антенна на передачу
        	// мощность на передачу
		mPowerPrd = StrTo_Dbl_(LabeledEdit40->Text);
		// КУ на передачу
		 mKYPrd = StrTo_Dbl_(LabeledEdit39->Text);

		 mTransmitAnt.mKYPrd  =mKYPrd;
		 mTransmitAnt.mPowerPrd = mPowerPrd;
		///


		///




		mSigSins  =  StrTo_Dbl_(LabeledEdit35->Text) / 1000.;

		mSig_d_po_dt_Sins =   StrTo_Dbl_(LabeledEdit54->Text) / 1000.;

			// интервал между измерениями
		mMeasT = 1. / StrTo_Dbl_(LabeledEdit21->Text);
		mSinsDelayT = StrTo_Dbl_(LabeledEdit48->Text);

		mBius = TBius(mMeasT, mSinsDelayT );




	}
	else
	{
		// ширинадиаграммы
		// mWidthDgr =  StrTo_Dbl_(LabeledEdit3->Text) * M_PI/3000./2. ;
	}



	// к-во столбцов  излучателей

	mNumEmitCols = StrToInt(LabeledEdit9->Text);
	// к-во строк излучателей
	mNumEmitRows = StrToInt(LabeledEdit11->Text);
	// к-во столбцов АМ  (по горизонтали )
	mNumAMCols = StrToInt(LabeledEdit5->Text);
	// к-во строк АМ (по вертикали )
	mNumAMRows = StrToInt(LabeledEdit6->Text);
	// расстояние между АМ по вериткали
	mdAMRow = ((double)mNumEmitRows) * mdEmitRow;
	// расстояние между АМ по горизонтали
	mdAMCol = ((double)mNumEmitCols) * mdEmitCol;
	LabeledEdit8->Text = mdAMRow;
	LabeledEdit7->Text = mdAMCol;
	///








	wchar_t wcsStr[2] = {
		0
	};



	const double SigEmitNoise = mNoiseSKZ_5P10 / sqrt(((double)(8 * 16 * 28)));
	const double SigEmitAmplFact = mEtalonSigAmplFact_5P10 / sqrt
		(((double)(8 * 16 * 28)));
	mAM_2D = TAM_2D(mNumEmitCols, mNumEmitRows, mdEmitCol, mdEmitRow,
		SigEmitNoise, SigEmitAmplFact);

	mFar_2D = TFar_2D(mNumAMCols, mNumAMRows, mLambda, mdAMCol, mdAMRow,
		mAM_2D, mbarrAM);

	for (int i = 0; i < mNumAMCols * mNumAMRows; i++) {

		if (!mbarrAM[i]) {
			mFar_2D.mpAm2D[i].mOtklCoefUs = 1.;
			mFar_2D.mpAm2D[i].mSigEmitNoise = 0.;
		}
	}

   //
 //*************************************************************************************
 //******* СОЗДАНИЕ КОРАБЛЯ ******************************************************************************
 //*************************************************************************************
	  // 3.3 создание привода
	mDriverSigBet  =      0.0021 ;
	mDriverSigEps  =      0.00021 ;
	mDriverDynamicSigBet =      0.003141;
	mDriverDynamicSigEps =      0.003141;

	 // 3.4  параметры корабля
	   // константы
	mVesselWidth = 0; // ширина(м)
	mVesselLength = 0; // длина(м)



	// парамеитры движения
	mQ0 = StrTo_Dbl_(LabeledEdit25->Text) * M_PI / 180. ; // генеральный курс
	mVVess = StrTo_Dbl_(LabeledEdit26->Text)* 0.514  ;// скорость корабля своего 20 узлов


	mMaxAmp_AftFlexure  = StrTo_Dbl_(LabeledEdit22->Text)* M_PI/180.;
	// период колебаний кормового изгиба
	mT_AftFlexure =  StrTo_Dbl_(LabeledEdit24->Text);
	//максимальная амплитуда бортового изгиба корабля в рад на 100 м
	mMaxAmp_BoardFlexure =   StrTo_Dbl_(LabeledEdit20->Text)* M_PI/180.;
	// период колебаний бортового изгиба
	mT_BoardFlexure = StrTo_Dbl_(LabeledEdit23->Text);


	/////////

   	 // 3.1 создание СИНС
	mMaxSig_Q      =     mSigSins ; //0.000582;
	mMaxSig_Psi    =      mSigSins ; //0.00145;
	mMaxSig_Tet    =    mSigSins ; // 0.00145;
	mMaxSig_dQdt   =      mSig_d_po_dt_Sins ;
	mMaxSig_dPsidt =      mSig_d_po_dt_Sins ;
	mMaxSig_dTetdt =      mSig_d_po_dt_Sins ;


		for (int i = 0; i < 3; i++)
		{
			marrFarParallacs[i] = StrTo_Dbl_(StringGrid1->Cells[i + 1][1] );
			marrArtParral[i] = StrTo_Dbl_(StringGrid1->Cells[i + 1][2]);
		}


}

//---------------------------------------------------------------------------
void __fastcall TForm5::Button1Click(TObject *Sender)
{
Close();
}
//---------------------------------------------------------------------------
void __fastcall TForm5::create5P10()
{
			// длина волны
		 mLambda = VAL_C /8./ 10000000.;
		// расстояние между излучателями
		mdEmitCol = mLambda * 0.55;
		mdEmitRow = mLambda * 0.55;


		// ЭТАЛОННЫЙ СИГНАЛ ДЛЯ РАСЧЕТА АМПЛИТУДЫ
		// амплитуда
		mEtalonAmp =   500.;
		// дальность
		mEtalonDist = 12000.;
		// ЭПР
		mEtalonAPR = 1.;
			// СКО внутр шума суммарной диаграммы 5П10
		mNoiseSKZ_5P10 = 18.735;
			// СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
		mEtalonSigAmplFact_5P10 = 0.01;
		// мощность на передачу
		mEtalonPowerPrd = 4000.;
    	// КУ на передачу
		mEtalonKYPrd = 840.;
			// КУ на прием
		mEtalonKYPriem =  5200.;


		mEtalonSign = TEtalonSign(mEtalonAmp, mEtalonDist, mEtalonAPR,
			mNoiseSKZ_5P10, mEtalonSigAmplFact_5P10, mEtalonPowerPrd, mEtalonKYPrd
          , mEtalonKYPriem);
		///

		//  антенна на передачу
					// мощность на передачу
					mPowerPrd = 4000.;
					// КУ на передачу
					mKYPrd =  3000.;
		 mTransmitAnt.mKYPrd  =mKYPrd;
		 mTransmitAnt.mPowerPrd = mPowerPrd;
		///


		///

		mSigSins  = 0.41 * 0.001;
		mSig_d_po_dt_Sins = 1.16 * 0.001;



			// интервал между измерениями
			mMeasT = 0.02;
			// задержка СИНС
			mSinsDelayT = 0.02;


		mBius = TBius(mMeasT, mSinsDelayT );

		marrArtParral[0] = 0.;
		marrArtParral[1] = 10.;
		marrArtParral[2] = 20.;

		marrFarParallacs[0] = 5.;
		marrFarParallacs[1] = 5.;
		marrFarParallacs[2] = 20.;
		for (int i = 0; i < 3; i++)
		{
			StringGrid1->Cells[i + 1][1] = marrFarParallacs[i];
			StringGrid1->Cells[i + 1][2] = marrArtParral[i];
		}






	// к-во столбцов  излучателей

	mNumEmitCols = 8;
	// к-во строк излучателей
	mNumEmitRows = 8;
	// к-во столбцов АМ  (по горизонтали )
	mNumAMCols = 8;
	// к-во строк АМ (по вертикали )
	mNumAMRows = 8;
	// расстояние между АМ по вериткали
	mdAMRow = ((double)mNumEmitRows) * mdEmitRow;
	// расстояние между АМ по горизонтали
	mdAMCol = ((double)mNumEmitCols) * mdEmitCol;

	///








	wchar_t wcsStr[2] = {
		0
	};



	const double SigEmitNoise = mNoiseSKZ_5P10 / sqrt(((double)(8 * 16 * 28)));
	const double SigEmitAmplFact = mEtalonSigAmplFact_5P10 / sqrt
		(((double)(8 * 16 * 28)));
	mAM_2D = TAM_2D(mNumEmitCols, mNumEmitRows, mdEmitCol, mdEmitRow,
		SigEmitNoise, SigEmitAmplFact);

	mFar_2D = TFar_2D(mNumAMCols, mNumAMRows, mLambda, mdAMCol, mdAMRow,
		mAM_2D, mbarrAM);

	for (int i = 0; i < mNumAMCols * mNumAMRows; i++) {

		if (!mbarrAM[i]) {
			mFar_2D.mpAm2D[i].mOtklCoefUs = 1.;
			mFar_2D.mpAm2D[i].mSigEmitNoise = 0.;
		}
	}









   //
 //*************************************************************************************
 //******* СОЗДАНИЕ КОРАБЛЯ ******************************************************************************
 //*************************************************************************************



		// 3.3 создание привода
	mDriverSigBet  =      0.0021 ;
	mDriverSigEps  =      0.00021 ;
	mDriverDynamicSigBet =      0.003141;
	mDriverDynamicSigEps =      0.003141;

	 // 3.4  параметры корабля
	   // константы
	mVesselWidth = 200.; // ширина(м)
	mVesselLength = 30.; // длина(м)

	memset(marrFarParallacs ,0, 3 * sizeof(double)) ;//  вектор параллакса
	mMaxQ =     3./180.*M_PI; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	mT_Q = 18.; // период рыскания
	mMaxPsi =      3./180.*M_PI;// максимальный угол килевой качки(амплитуда)
	mT_Psi = 12; // период килевой качки
	mMaxTet =      12./180.*M_PI; //максимальный угол боротовой качки(амплитуда)
	mT_Tet = 6; // период бортовой качки
	mMaxVert =     1. ;

	// парамеитры движения
	mQ0 = 0. ; // генеральный курс
	mVVess =  0.514 * 20 ;// скорость корабля своего 20 узлов


	mMaxAmp_AftFlexure  = 1. * M_PI/180.;
	// период колебаний кормового изгиба
	mT_AftFlexure = 4.;
	//максимальная амплитуда бортового изгиба корабля в рад на 100 м
	mMaxAmp_BoardFlexure =  1. * M_PI/180.;
	// период колебаний бортового изгиба
	mT_BoardFlexure = 2.;


	/////////

	 // 3.1 создание СИНС
	 mMaxSig_Q      =     mSigSins ; //0.000582;
	mMaxSig_Psi    =      mSigSins ; //0.00145;
	mMaxSig_Tet    =    mSigSins ; // 0.00145;
	mMaxSig_dQdt   =      mSig_d_po_dt_Sins ;
	mMaxSig_dPsidt =      mSig_d_po_dt_Sins ;
	mMaxSig_dTetdt =      mSig_d_po_dt_Sins ;
	mK1         = 0.01 ;
	mSigV       =      0.2 * sqrt(2.) ;
	mSigH       =      0.1 ;
	mMaxSig_H =     0.1 ;
	mMaxSig_VH =     0.05 ;


	//	mLambda = VAL_C / StrTo_Dbl_(LabeledEdit15->Text) * 100. / 1000000000.;
		 LabeledEdit15->Text = VAL_C /  mLambda / 1000000000.;

		LabeledEdit12->Text = mdEmitCol;

		LabeledEdit2->Text = mEtalonAmp ;
		// дальность
		LabeledEdit3->Text = mEtalonDist / 1000.;
		// ЭПР
		LabeledEdit4->Text = mEtalonAPR ;

		LabeledEdit16->Text = mNoiseSKZ_5P10 ;
		// СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
		LabeledEdit18->Text = mEtalonSigAmplFact_5P10 ;

		// мощность на передачу
		LabeledEdit36->Text = mEtalonPowerPrd ;
		// КУ на передачу
		LabeledEdit38->Text = mEtalonKYPrd ;
		// мощность на прием
		LabeledEdit37->Text = mEtalonKYPriem ;

		LabeledEdit40->Text = mPowerPrd ;
		// КУ на передачу
		 LabeledEdit39->Text = mKYPrd;

		 LabeledEdit35->Text =	mSigSins  * 1000.;

		 LabeledEdit54->Text = mSig_d_po_dt_Sins * 1000.;

		LabeledEdit21->Text = 1. /mMeasT ;
		LabeledEdit48->Text = mSinsDelayT;


		// к-во столбцов  излучателей

	LabeledEdit9->Text = mNumEmitCols ;
	// к-во строк излучателей
	LabeledEdit11->Text = mNumEmitRows ;
	// к-во столбцов АМ  (по горизонтали )
	LabeledEdit5->Text = mNumAMCols ;
	// к-во строк АМ (по вертикали )
	LabeledEdit6->Text = mNumAMRows ;

	LabeledEdit8->Text = mdAMRow;
	LabeledEdit7->Text = mdAMCol;
	///


	LabeledEdit7->Text =	mDriverSigBet * 1000. ;
	LabeledEdit10->Text = mDriverSigEps  * 1000. ;
	LabeledEdit17->Text =	mDriverDynamicSigBet * 1000. ;
	LabeledEdit19->Text =	mDriverDynamicSigEps * 1000. ;

	 // 3.4  параметры корабля
	   // константы

	// парамеитры движения
	LabeledEdit19->Text = mQ0 * 180. / M_PI ; // генеральный курс
	LabeledEdit26->Text = mVVess/ 0.514 ;// скорость корабля своего 20 узлов


	LabeledEdit22->Text = mMaxAmp_AftFlexure  * 180. / M_PI ;
	// период колебаний кормового изгиба
	LabeledEdit20->Text = mMaxAmp_BoardFlexure * 180. / M_PI ;// mT_AftFlexure;
	//максимальная амплитуда бортового изгиба корабля в рад на 100 м
	LabeledEdit24->Text =  mT_AftFlexure;
	// период колебаний бортового изгиба
	LabeledEdit24->Text = mT_BoardFlexure ;



}
//---------------------------------------------------------------------------

