//---------------------------------------------------------------------------

#ifndef FormVessH
#define FormVessH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Grids.hpp>
//---------------------------------------------------------------------------




//---------------------------------------------------------------------------

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
	TPanel *Panel1;
	TLabel *Label3;
	TStringGrid *StringGrid1;
	TButton *Button1;
	void __fastcall LabeledEdit26Change(TObject *Sender);
	void __fastcall Button1Click(TObject *Sender);
private:	// User declarations
public:		// User declarations

	__fastcall TForm5(TComponent* Owner);

			// ЭТАЛОННЫЙ СИГНАЛ ДЛЯ РАСЧЕТА АМПЛИТУДЫ
	//амплитуда
	double *mpEtalonAmp;
	//дальность
	double *mpEtalonDist;
	//ЭПР
	double *mpEtalonAPR;
	//СКО внутр шума суммарной диаграммы 5П10
	double *mpNoiseSKZ_5P10;
	// СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
	double *mpEtalonSigAmplFact_5P10;
	//
	// мощность на передачу
	double *mpEtalonPowerPrd;
	// КУ на передачу
	double *mpEtalonKYPrd;
	// мощность на прием
	double *mpEtalonKYPriem;

	///

	// АНТЕННА
	// данные по АМ
	// к-во излучателей по горизонтали
	int *mpNumEmitCols;
	// к-во излучателей по вертикали
	int *mpNumEmitRows;
	// расстояние между излучателями по горизонтали
	double *mpdEmitCol;
	// расстояние между излучателями по вериткали
	double *mpdEmitRow;
	// длина волны
	double *mpLambda;


	// данные по ФАР
	// к-во АМ  по горизонтали
	int *mpNumAMCols;
	// к-во АМ по вертикали
	int *mpNumAMRows;
	// расстояние между АМ по горизонтали
	double *mpdAMCol;
	// расстояние между АМ по вериткали
	double *mpdAMRow;
	// СКЗ шума в суммарной диаграмме
	double *mpSigNoise;



	// мощность на передачу
	double *mpPowerPrd;
	// КУ на передачу
	double *mpKYPrd;

	// индикаторы рабочих модуленй в модельной ФАР (которая заложена в расчетах)
	bool *mbarrAM;
  //-------------------------------------------------------------------
	// СКЗ угловой ошибки СИНС (углов качек)
	double *mpSigSins;

	// СКЗ угловой ошибки СИНСпо скорости углов (углов качек)
	double *mpSig_d_po_dt_Sins;

	// корабль наш
	double *mpVesselWidth ; // ширина(м)
	double *mpVesselLength ;
	double *mparrFarParallacs ;

	double *mpMaxQ ; /// максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	double *mpT_Q; // период рыскания
	double *mpMaxPsi ;// максимальный угол килевой качки(амплитуда)
	double *mpT_Psi ; // период килевой качки
	double *mpMaxTet ; //максимальный угол боротовой качки(амплитуда)
	double *mpT_Tet; // период бортовой качки
	double *mpMaxVert  ;

	// парамеитры движения  корабля нашего
	double *mpQ0  ; // генеральный курс
	double *mpVVess  ;// скорость корабля своего 20 узлов
	double *mparrDelt ;//  начальные фазы

	//максимальная амплитуда кормового изгиба корабля в рад на 100 м
	double *mpMaxAmp_AftFlexure;
	// период колебаний кормового изгиба
	double *mpT_AftFlexure;
	//максимальная амплитуда бортового изгиба корабля в рад на 100 м
	double *mpMaxAmp_BoardFlexure;
	// период колебаний бортового изгиба
	double *mpT_BoardFlexure;

	// 3.1 создание СИНС

	double *mpMaxSig_Q ;
	double *mpMaxSig_Psi  ;
	double *mpMaxSig_Tet  ;
	double *mpMaxSig_dQdt ;
	double *mpMaxSig_dPsidt ;
	double *mpMaxSig_dTetdt ;
	double *mpK1         ;
	double *mpSigV      ;
	double *mpSigH     ;
	double *mpMaxSig_H ;
	double *mpMaxSig_VH ;
	// тип алгоритма фильтрации
	//	bool mbSkaliga;

	// привод
	double *mpDriverSigBet ;// точность измерения угла Bet привода
	double *mpDriverSigEps ;// точность измерения угла Eps  привода (угла места)
	double *mpDriverDynamicSigBet ;// точность отработки угла курса  привода
	double *mpDriverDynamicSigEps ;// точность  привода отработки угла места



	//

	//Темп фильтрации
	// интервал между измерениями
	double *mpMeasT;

	// Задержка СИНС
	double *mpSinsDelayT;
	double *mparrArtParral; // вектор параллакса АУ



	// 	точногсть отработки приводом скорости  углов
	double  *mpSigDrivAY_dU_po_dt ;



};
//---------------------------------------------------------------------------
extern PACKAGE TForm5 *Form5;
//---------------------------------------------------------------------------
#endif
