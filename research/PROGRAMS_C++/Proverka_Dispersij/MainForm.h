//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Dialogs.hpp>
#include <ExtCtrls.hpp>

#include "FormCoastTargs.h"
#include "FormAeroTargs.h"
#include "FormSeaTargs.h"
//#include "FormVessel.h"
//#include "FormVess.h"

#include "AM_2D.h"
#include "Far_2D.h"
#include "URPointZ.h"
#include "SingleSign.h"
#include <Buttons.hpp>
#include <Dialogs.hpp>

#include "InitTargData.h"
#include "EtalonSign.h"
//#include "Target.h"
#include "TransmitAnt.h"


//#include "Environment.h"
//#include "MyShellTraj.h"
//#include "ControlSyst.h"
//#include "ArtCannon.h"
//#include "ArtComplex.h"
//#include "Vessel.h"
//#include "Fight.h"
//#include "ArtCannon.h"
#include <ButtonGroup.hpp>
#include <Grids.hpp>


#define MAX_QUANT_AM 1600
 enum enumTargetType;
 enum enumShellType;
 enum enumCannonType;
//---------------------------------------------------------------------------
class TAM_2D;
class TFar_2D;
class TURPointZ;
class  TSingleSign;
class TInitTargData;
class TEtalonSign;
class TTransmitAnt;
//class TWind;
//class TEnvironment;
//class TControlSyst;
//class TArtCannon;
//class TArtComplex;
//class TVessel;
//class TFight;
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TOpenDialog *OpenDialog1;
	TPanel *Panel1;
	TLabel *Label3;
	TButton *Button3;
	TButton *Button4;
	TEdit *Edit2;
	TPanel *Panel6;
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
	TLabeledEdit *LabeledEdit39;
	TLabeledEdit *LabeledEdit40;
	TPanel *Panel8;
	TLabel *Label1;
	TLabeledEdit *LabeledEdit4;
	TLabeledEdit *LabeledEdit10;
	TLabeledEdit *LabeledEdit16;
	TLabeledEdit *LabeledEdit17;
	TLabeledEdit *LabeledEdit18;
	TLabeledEdit *LabeledEdit38;
	TLabeledEdit *LabeledEdit37;
	TLabeledEdit *LabeledEdit36;
	TLabeledEdit *LabeledEdit28;
	TLabeledEdit *LabeledEdit29;
	TLabel *Label4;
	TLabeledEdit *LabeledEdit1;
	TLabeledEdit *LabeledEdit2;
	TLabeledEdit *LabeledEdit3;
	TLabeledEdit *LabeledEdit19;
	TButton *Button1;
	TLabeledEdit *LabeledEdit20;
	TLabeledEdit *LabeledEdit21;
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Button3Click(TObject *Sender);
	void __fastcall Button4Click(TObject *Sender);
	void __fastcall Panel7Click(TObject *Sender);
	void __fastcall Panel6Click(TObject *Sender);
	void __fastcall Panel2Click(TObject *Sender);
	void __fastcall Panel3Click(TObject *Sender);
	void __fastcall Panel8Click(TObject *Sender);
	void __fastcall Panel9Click(TObject *Sender);
	void __fastcall Panel10Click(TObject *Sender);


	void __fastcall Button2Click(TObject *Sender);

	void __fastcall Button7Click(TObject *Sender);
	void __fastcall LabeledEdit15Change(TObject *Sender);
	void __fastcall LabeledEdit15Exit(TObject *Sender);
	void __fastcall LabeledEdit31Change(TObject *Sender);
	void __fastcall LabeledEdit21Change(TObject *Sender);
	void __fastcall LabeledEdit48Change(TObject *Sender);
	void __fastcall LabeledEdit54Change(TObject *Sender);
	void __fastcall LabeledEdit35Change(TObject *Sender);
	void __fastcall LabeledEdit40Change(TObject *Sender);
	void __fastcall LabeledEdit39Change(TObject *Sender);
	void __fastcall LabeledEdit29Change(TObject *Sender);
	void __fastcall LabeledEdit28Change(TObject *Sender);
	void __fastcall LabeledEdit1Change(TObject *Sender);

private:	// User declarations
public:		// User declarations
	__fastcall TForm1(TComponent* Owner);

	// ���� � ����� � ���������
	wchar_t *mpwchOutFile0;
	// ���� � ����� � ���������  ����������� ������� ��������
	wchar_t mwchOutFold[400]; //
	// ��������� ������ ��� ������� ���������
				// ��������� ������ ��� ������� ���������
	//���������
	double mEtalonAmp;
	//���������
	double mEtalonDist;
	//���
	double mEtalonAPR;
	//��� ����� ���� ��������� ��������� 5�10
	double mNoiseSKZ_5P10;
	// ��� �������� ������� �������� ��������� ��������� 5�10
	double mEtalonSigAmplFact_5P10;
	//
	// �������� �� ��������
	double mEtalonPowerPrd;
	// �� �� ��������
	double mEtalonKYPrd;
	// �������� �� �����
	double mEtalonKYPriem;

	///

	// �������
	// ������ �� ��
	// �-�� ����������� �� �����������
	int mNumEmitCols;
	// �-�� ����������� �� ���������
	int mNumEmitRows;
	// ���������� ����� ������������ �� �����������
	double mdEmitCol;
	// ���������� ����� ������������ �� ���������
	double mdEmitRow;
	// ����� �����
	double mLambda;


	// ������ �� ���
	// �-�� ��  �� �����������
	int mNumAMCols;
	// �-�� �� �� ���������
	int mNumAMRows;
	// ���������� ����� �� �� �����������
	double mdAMCol;
	// ���������� ����� �� �� ���������
	double mdAMRow;
	// ��� ���� � ��������� ���������
	double mSigNoise;
	 bool mbarrAM[MAX_QUANT_AM];
 // ����� R
 double mSigmaR;

	 double marrArtParral[3]; // ������ ���������� ��


	// �������� �� ��������
	double mPowerPrd;
	// �� �� ��������
	double mKYPrd;


	TEtalonSign mEtalonSign;
	///


	// �������� ������
	TAM_2D mAM_2D;


	// ��� ���������� � �������
	TFar_2D mFar_2D;


	// ������� �� ��������
	TTransmitAnt mTransmitAnt ;

	// ���� ��������
	double mAntpCoef;

	// ������ ����������
	double mTrajectH;

	// ��� ����
	double mEPR;

	// ������� �������
	double mAntH;

	//-------------------------------------------------------------------


	void __fastcall fncInputData();

 //	void __fastcall ExchangeData();

	void __fastcall create5P10();

	void __fastcall fillVessel();

	double  __fastcall calcDispNew(double eps, double dispSum, double A, double lamb,double d
, double valFDiagr) ;

};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif