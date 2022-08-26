//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Dialogs.hpp>
#include <Grids.hpp>
#include "Comp.h"
class TComp;
class TDiagrSet;
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TPanel *Panel2;
	TEdit *Edit1;
	TButton *Button1;
	TButton *Button2;
	TLabel *Label2;
	TSaveDialog *SaveDialog1;
	TPanel *Panel1;
	TLabel *Label1;
	TLabel *Label4;
	TLabel *����������;
	TLabeledEdit *LabeledEdit1;
	TLabeledEdit *LabeledEdit2;
	TLabeledEdit *LabeledEdit3;
	TLabeledEdit *LabeledEdit4;
	TLabeledEdit *LabeledEdit9;
	TLabeledEdit *LabeledEdit10;
	TStringGrid *StringGrid2;
	TEdit *Edit2;
	TEdit *Edit3;
	TStringGrid *StringGrid3;
	TStringGrid *StringGrid4;
	TStringGrid *StringGrid5;
	TComboBox *ComboBox1;
	TLabel *Label3;
	void __fastcall Button2Click(TObject *Sender);
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall LabeledEdit4Change(TObject *Sender);
	void __fastcall LabeledEdit5Change(TObject *Sender);
	void __fastcall LabeledEdit6Change(TObject *Sender);
	void __fastcall LabeledEdit11Change(TObject *Sender);


	void __fastcall LabeledEdit15Change(TObject *Sender);
	void __fastcall LabeledEdit16Change(TObject *Sender);
	void __fastcall LabeledEdit17Change(TObject *Sender);
	void __fastcall LabeledEdit18Change(TObject *Sender);

	void __fastcall Edit3Exit(TObject *Sender);
	void __fastcall StringGrid2Exit(TObject *Sender);
	void __fastcall createInputGridForTarg_and_Antp();
	void __fastcall createOutputGridForTarg_and_Antp();
	void __fastcall createOutputGridForTrueAnglesUM_Targ_and_Antp();
	void __fastcall createOutputGridForRezAnglesUM_Targ_and_Antp();


	void __fastcall LabeledEdit1Exit(TObject *Sender);
	void __fastcall LabeledEdit2Exit(TObject *Sender);
	void __fastcall LabeledEdit3Exit(TObject *Sender);
	void __fastcall LabeledEdit4Exit(TObject *Sender);

private:	// User declarations
public:		// User declarations
// ���� � ����� � ���������
wchar_t *mpwchOutFile0; //

// ���� � ����� � ���������
wchar_t mwchOutFold[400]; //

//1. �������� ������� ��������� �� ������ 0,707, �.�.
double mWidthDgr;
// 2.��� ����
double mNoiseSkz;
//3. ������ �������
double mAltCoord ;
//4. �-�� ���������
int mNumDiagr;
//5. ������ ����
double  mAltTarg;
//6.������ ��������
double  mAltAntp;
// 7.��������� ����
double mDistTarg;
// 8.������ ����, �������
double mSignLevTarg;
// 9.������ ��������, �������
double  mSignLevAntp;
// 10. ������ ����� ����� ��������
double marrAlfaDiagr[1000];

// 11. ������� ���������� ����
TComp mcmpKTarg;
//12. ����� ��������� ��������
TComp mcmpKAntp;
//13. ���� ���� ��������� ����
double mPhaseTarg;
//14. ���� ���� ��������� ��������
double mPhaseAntp;
// 15.����� �����
double mLambda;
// 16. ������ ������� ��������� ������ ��������
double mUpAngVeer;
// 16. ������� ������� ��������� ������ ��������
double mDownAngVeer;

// ��� ������ ��������
int miTypeDiagrSet;

	__fastcall TForm1(TComponent* Owner);
	void __fastcall fncInputData();

double  __fastcall calcCritUMTarg();
// ���� ������ ���� ������ ����, �� ������� �������� �� �����
double  __fastcall calcCritElevTarg() ;

 // ���������� ��������� ���� ����� �������� ������������ ������������
double  __fastcall calcMinUMAntp()  ;
// ���������� ���� ����� ���� ��������� ������������ ������������
double  __fastcall calcUMTarg();

// ���������� ���� ����� ��������� ��������� ������������ ������������
double  __fastcall calcUMAntp() ;

void __fastcall fillInputGrids();
//void __fastcall createInputGridForVeerDiagrams() ;
//void __fastcall fillVeerDiagrams();
//void __fastcall updateVeerDiagrams();
void __fastcall updatRezGrids();
void __fastcall updatRezTrueAngUMGrid() ;
void __fastcall InputGridDataTarg_and_Antp();
void __fastcall CreateDiagrSet( TDiagrSet *DiagrSet0);


};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif