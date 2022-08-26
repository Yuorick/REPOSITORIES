//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <CheckLst.hpp>
#include "AM_2D.h"
#include "Far_2D.h"
#include "URPointZ.h"
#include "SingleSign.h"
#include <Buttons.hpp>
#include <Dialogs.hpp>
//#include "TargBearing0.h"
#include "InitTargData.h"
#include "EtalonSign.h"
#include "Target.h"
#include "TransmitAnt.h"
#include <ComCtrls.hpp>
//#include "ManulEffectImproved_v4.h"

#define MAX_QUANT_AM 1600

//---------------------------------------------------------------------------
class TAM_2D;
class TFar_2D;





class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TLabeledEdit *LabeledEdit1;
	TPanel *Panel3;
	TButton *Button2;
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
	TLabel *Label2;
	TPanel *Panel2;
	TLabel *Label1;
	TLabeledEdit *LabeledEdit2;
	TLabeledEdit *LabeledEdit3;
	TLabeledEdit *LabeledEdit16;
	TLabeledEdit *LabeledEdit21;
	TPanel *Panel7;
	TLabeledEdit *LabeledEdit22;
	TLabeledEdit *LabeledEdit23;
	TLabeledEdit *LabeledEdit24;
	TLabeledEdit *LabeledEdit18;
	TLabeledEdit *LabeledEdit30;
	TLabeledEdit *LabeledEdit31;
	TLabeledEdit *LabeledEdit42;
	TLabeledEdit *LabeledEdit39;
	TLabeledEdit *LabeledEdit40;
	TLabeledEdit *LabeledEdit33;
	TLabeledEdit *LabeledEdit34;
	TPanel *Panel1;
	TOpenDialog *OpenDialog1;
	TButton *Button3;
	TButton *Button4;
	TEdit *Edit2;
	TLabel *Label3;
	TLabeledEdit *LabeledEdit4;
	TLabeledEdit *LabeledEdit32;
	TLabeledEdit *LabeledEdit36;
	TLabeledEdit *LabeledEdit37;
	TLabel *Label4;
	TLabeledEdit *LabeledEdit25;
	void __fastcall Button2Click(TObject *Sender);
	void __fastcall FormShow(TObject *Sender);
	void __fastcall LabeledEdit2Exit(TObject *Sender);
	void __fastcall LabeledEdit3Exit(TObject *Sender);
	void __fastcall LabeledEdit4Exit(TObject *Sender);
	void __fastcall LabeledEdit16Exit(TObject *Sender);
	void __fastcall LabeledEdit18Exit(TObject *Sender);
	void __fastcall LabeledEdit36Exit(TObject *Sender);
	void __fastcall LabeledEdit38Exit(TObject *Sender);
	void __fastcall LabeledEdit37Exit(TObject *Sender);
	void __fastcall ComboBox1Exit(TObject *Sender);
	void __fastcall ComboBox2Exit(TObject *Sender);
	void __fastcall LabeledEdit10Exit(TObject *Sender);
	void __fastcall LabeledEdit19Exit(TObject *Sender);
	void __fastcall LabeledEdit20Exit(TObject *Sender);
	void __fastcall LabeledEdit17Exit(TObject *Sender);
	void __fastcall LabeledEdit28Exit(TObject *Sender);
	void __fastcall LabeledEdit35Exit(TObject *Sender);
	void __fastcall ComboBox3Exit(TObject *Sender);
	void __fastcall LabeledEdit24Exit(TObject *Sender);
	void __fastcall LabeledEdit23Exit(TObject *Sender);
	void __fastcall LabeledEdit22Exit(TObject *Sender);
	void __fastcall LabeledEdit25Exit(TObject *Sender);
	void __fastcall LabeledEdit26Exit(TObject *Sender);
	void __fastcall LabeledEdit29Exit(TObject *Sender);
	void __fastcall LabeledEdit9Exit(TObject *Sender);
	void __fastcall LabeledEdit11Exit(TObject *Sender);
	void __fastcall LabeledEdit15Exit(TObject *Sender);
	void __fastcall LabeledEdit12Exit(TObject *Sender);
	void __fastcall LabeledEdit13Exit(TObject *Sender);
	void __fastcall LabeledEdit5Exit(TObject *Sender);
	void __fastcall LabeledEdit6Exit(TObject *Sender);
	void __fastcall LabeledEdit14Exit(TObject *Sender);
	void __fastcall LabeledEdit8Exit(TObject *Sender);
	void __fastcall LabeledEdit7Exit(TObject *Sender);
	void __fastcall LabeledEdit40Exit(TObject *Sender);
	void __fastcall LabeledEdit39Exit(TObject *Sender);
	void __fastcall LabeledEdit21Exit(TObject *Sender);
	void __fastcall LabeledEdit30Exit(TObject *Sender);

	void __fastcall Button3Click(TObject *Sender);
	void __fastcall Button4Click(TObject *Sender);


	

private:	// User declarations
public:		// User declarations
	__fastcall TForm1(TComponent* Owner);

		// путь к папке с нрафиками
	 wchar_t *mpwchOutFile0;
	 // путь к папке с графиками  результатов расчета диаграмм
	wchar_t mwchOutFold[400]; //
 // ЭТАЛОННЫЙ СИГНАЛ ДЛЯ РАСЧЕТА АМПЛИТУДЫ
 //амплитуда  сигнала цели
 double mAmpTarg;
 //амплитуда  сигнала антипода
 double mAmpAntp;

 //СКО внутр шума суммарной диаграммы 5П10
 double mNoiseSKZ_5P10;
 // СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
 double mEtalonSigAmplFact_5P10;
 //
 	// истинное угловое положение цели по УМ
	double mEpsTargTrue;
// истинное угловое положение цели по КУ
	double mBetTargTrue;
	// истинное угловое положение ант по УМ
	double mEpsAntpTrue;
// истинное угловое положение ант по КУ
	double mBetAntpTrue;




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



	// индикаторы рабочих модуленй в модельной ФАР (которая заложена в расчетах)
	bool mbarrAM[MAX_QUANT_AM];








	void __fastcall fncInputData();
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------

#endif
