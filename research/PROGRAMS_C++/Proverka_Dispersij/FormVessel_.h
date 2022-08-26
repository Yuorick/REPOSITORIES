//---------------------------------------------------------------------------

#ifndef FormVessel_H
#define FormVessel_H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------


#include "SingleSign.h"
#include "InitTargData.h"
#include "EtalonSign.h"
#include "Target.h"
#include "TransmitAnt.h"
#include "AM_2D.h"
#include "Far_2D.h"
#include "URPointZ.h"

#include "Environment.h"
#include "MyShellTraj.h"
#include "Bius.h"
#include "ArtCannon.h"
#include "ArtComplex.h"
#include "Vessel.h"
#include "Fight.h"
#include "ArtCannon.h"
#include <Grids.hpp>



#define MAX_QUANT_AM 1600

//---------------------------------------------------------------------------
class TAM_2D;
class TFar_2D;
class TURPointZ;
class  TSingleSign;

class TEtalonSign;
class TTransmitAnt;
class TWind;
class TEnvironment;
class TBius;
class TVessel;

class TForm5 : public TForm
{
__published:	// IDE-managed Components
	TPanel *Panel6;
	TStaticText *StaticText1;
	TPanel *Panel3;
	TLabel *Label2;
	TPanel *Panel4;
	TLabeledEdit *LabeledEdit5;
	TLabeledEdit *LabeledEdit6;
	TLabeledEdit *LabeledEdit7;
	TLabeledEdit *LabeledEdit8;
	TLabeledEdit *LabeledEdit14;
	TPanel *Panel5;
	TLabeledEdit *LabeledEdit9;
	TLabeledEdit *LabeledEdit11;
	TLabeledEdit *LabeledEdit12;
	TLabeledEdit *LabeledEdit13;
	TLabeledEdit *LabeledEdit15;
	TLabeledEdit *LabeledEdit42;
	TLabeledEdit *LabeledEdit34;
	TPanel *Panel7;
	TLabel *Label4;
	TLabeledEdit *LabeledEdit39;
	TLabeledEdit *LabeledEdit40;
	TPanel *Panel2;
	TLabel *Label1;
	TLabeledEdit *LabeledEdit2;
	TLabeledEdit *LabeledEdit3;
	TLabeledEdit *LabeledEdit4;
	TLabeledEdit *LabeledEdit16;
	TLabeledEdit *LabeledEdit18;
	TLabeledEdit *LabeledEdit38;
	TLabeledEdit *LabeledEdit37;
	TLabeledEdit *LabeledEdit36;
	TPanel *Panel10;
	TLabel *Label5;
	TLabeledEdit *LabeledEdit54;
	TLabeledEdit *LabeledEdit35;
	TPanel *Panel9;
	TLabel *Label9;
	TLabeledEdit *LabeledEdit48;
	TLabeledEdit *LabeledEdit21;
	TPanel *Panel12;
	TLabel *Label6;
	TLabeledEdit *LabeledEdit1;
	TLabeledEdit *LabeledEdit10;
	TLabeledEdit *LabeledEdit17;
	TLabeledEdit *LabeledEdit19;
	TPanel *Panel13;
	TLabel *Label7;
	TLabeledEdit *LabeledEdit20;
	TLabeledEdit *LabeledEdit22;
	TLabeledEdit *LabeledEdit23;
	TLabeledEdit *LabeledEdit24;
	TLabeledEdit *LabeledEdit25;
	TLabeledEdit *LabeledEdit26;
	TButton *Button1;
	TPanel *Panel1;
	TLabel *Label3;
	TStringGrid *StringGrid1;
	void __fastcall Button1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TForm5(TComponent* Owner);

		// ЭТАЛОННЫЙ СИГНАЛ ДЛЯ РАСЧЕТА АМПЛИТУДЫ
	//амплитуда
	double mEtalonAmp;
	//дальность
	double mEtalonDist;
	//ЭПР
	double mEtalonAPR;
	//СКО внутр шума суммарной диаграммы 5П10
	double mNoiseSKZ_5P10;
	// СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
	double mEtalonSigAmplFact_5P10;
	//
	// мощность на передачу
	double mEtalonPowerPrd;
	// КУ на передачу
	double mEtalonKYPrd;
	// мощность на прием
	double mEtalonKYPriem;
	TEtalonSign mEtalonSign;
	///

	// АНТЕННА
	// данные по АМ
	// к-во излучателей по горизонтали
	int mNumEmitCols;
	// к-во излучателей по вертикали
	int mNumEmitRows;
	// расстояние между излучателями по горизонтали
	double mdEmitCol;
	// расстояние между излучателями по вериткали
	double mdEmitRow;
	// длина волны
	double mLambda;
	// антенный модуль
	TAM_2D mAM_2D;

	// данные по ФАР
	// к-во АМ  по горизонтали
	int mNumAMCols;
	// к-во АМ по вертикали
	int mNumAMRows;
	// расстояние между АМ по горизонтали
	double mdAMCol;
	// расстояние между АМ по вериткали
	double mdAMRow;
	// СКЗ шума в суммарной диаграмме
	double mSigNoise;

	// ФАР заложенная в расчеты
	TFar_2D mFar_2D;

	// мощность на передачу
	double mPowerPrd;
	// КУ на передачу
	double mKYPrd;
	// антенна на передачу
	TTransmitAnt mTransmitAnt ;
	// индикаторы рабочих модуленй в модельной ФАР (которая заложена в расчетах)
	bool mbarrAM[MAX_QUANT_AM];
  //-------------------------------------------------------------------
	// СКЗ угловой ошибки СИНС (углов качек)
	double	mSigSins;

	// СКЗ угловой ошибки СИНСпо скорости углов (углов качек)
	double	mSig_d_po_dt_Sins;

	// корабль наш
	double mVesselWidth ; // ширина(м)
	double mVesselLength ;
	double marrFarParallacs[3] ;

	double mMaxQ ; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	double mT_Q; // период рыскания
	double mMaxPsi ;// максимальный угол килевой качки(амплитуда)
	double mT_Psi ; // период килевой качки
	double mMaxTet ; //максимальный угол боротовой качки(амплитуда)
	double mT_Tet; // период бортовой качки
	double mMaxVert  ;

	// парамеитры движения  корабля нашего
	double mQ0  ; // генеральный курс
	double mVVess  ;// скорость корабля своего 20 узлов
	double marrDelt[4] ;//  начальные фазы

	//максимальная амплитуда кормового изгиба корабля в рад на 100 м
	double mMaxAmp_AftFlexure;
	// период колебаний кормового изгиба
	double mT_AftFlexure;
	//максимальная амплитуда бортового изгиба корабля в рад на 100 м
	double mMaxAmp_BoardFlexure;
	// период колебаний бортового изгиба
	double mT_BoardFlexure;

	// 3.1 создание СИНС
	TSins mSins ;
	double mMaxSig_Q ;
	double mMaxSig_Psi  ;
	double mMaxSig_Tet  ;
	double mMaxSig_dQdt ;
	double mMaxSig_dPsidt ;
	double mMaxSig_dTetdt ;
	double mK1         ;
	double mSigV      ;
	double mSigH     ;
	double mMaxSig_H ;
	double mMaxSig_VH ;
	// тип алгоритма фильтрации
	//	bool mbSkaliga;

	// привод
	double mDriverSigBet ;// точность измерения угла Bet привода
	double mDriverSigEps ;// точность измерения угла Eps  привода (угла места)
	double mDriverDynamicSigBet ;// точность отработки угла курса  привода
	double mDriverDynamicSigEps ;// точность  привода отработки угла места



	//

	//Темп фильтрации
	// интервал между измерениями
	double mMeasT;

	// Задержка СИНС
	double mSinsDelayT;


	TBius mBius;



	double marrArtParral[3]; // вектор параллакса АУ



	// 	точногсть отработки приводом скорости  углов
	double  mSigDrivAY_dU_po_dt ;

	void __fastcall fncInputData();

	void __fastcall create5P10();


};
//---------------------------------------------------------------------------
extern PACKAGE TForm5 *Form5;
//---------------------------------------------------------------------------
#endif
