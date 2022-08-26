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

			// ��������� ������ ��� ������� ���������
	//���������
	double *mpEtalonAmp;
	//���������
	double *mpEtalonDist;
	//���
	double *mpEtalonAPR;
	//��� ����� ���� ��������� ��������� 5�10
	double *mpNoiseSKZ_5P10;
	// ��� �������� ������� �������� ��������� ��������� 5�10
	double *mpEtalonSigAmplFact_5P10;
	//
	// �������� �� ��������
	double *mpEtalonPowerPrd;
	// �� �� ��������
	double *mpEtalonKYPrd;
	// �������� �� �����
	double *mpEtalonKYPriem;

	///

	// �������
	// ������ �� ��
	// �-�� ����������� �� �����������
	int *mpNumEmitCols;
	// �-�� ����������� �� ���������
	int *mpNumEmitRows;
	// ���������� ����� ������������ �� �����������
	double *mpdEmitCol;
	// ���������� ����� ������������ �� ���������
	double *mpdEmitRow;
	// ����� �����
	double *mpLambda;


	// ������ �� ���
	// �-�� ��  �� �����������
	int *mpNumAMCols;
	// �-�� �� �� ���������
	int *mpNumAMRows;
	// ���������� ����� �� �� �����������
	double *mpdAMCol;
	// ���������� ����� �� �� ���������
	double *mpdAMRow;
	// ��� ���� � ��������� ���������
	double *mpSigNoise;



	// �������� �� ��������
	double *mpPowerPrd;
	// �� �� ��������
	double *mpKYPrd;

	// ���������� ������� �������� � ��������� ��� (������� �������� � ��������)
	bool *mbarrAM;
  //-------------------------------------------------------------------
	// ��� ������� ������ ���� (����� �����)
	double *mpSigSins;

	// ��� ������� ������ ������ �������� ����� (����� �����)
	double *mpSig_d_po_dt_Sins;

	// ������� ���
	double *mpVesselWidth ; // ������(�)
	double *mpVesselLength ;
	double *mparrFarParallacs ;

	double *mpMaxQ ; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	double *mpT_Q; // ������ ��������
	double *mpMaxPsi ;// ������������ ���� ������� �����(���������)
	double *mpT_Psi ; // ������ ������� �����
	double *mpMaxTet ; //������������ ���� ��������� �����(���������)
	double *mpT_Tet; // ������ �������� �����
	double *mpMaxVert  ;

	// ���������� ��������  ������� ������
	double *mpQ0  ; // ����������� ����
	double *mpVVess  ;// �������� ������� ������ 20 �����
	double *mparrDelt ;//  ��������� ����

	//������������ ��������� ��������� ������ ������� � ��� �� 100 �
	double *mpMaxAmp_AftFlexure;
	// ������ ��������� ��������� ������
	double *mpT_AftFlexure;
	//������������ ��������� ��������� ������ ������� � ��� �� 100 �
	double *mpMaxAmp_BoardFlexure;
	// ������ ��������� ��������� ������
	double *mpT_BoardFlexure;

	// 3.1 �������� ����

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
	// ��� ��������� ����������
	//	bool mbSkaliga;

	// ������
	double *mpDriverSigBet ;// �������� ��������� ���� Bet �������
	double *mpDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
	double *mpDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
	double *mpDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����



	//

	//���� ����������
	// �������� ����� �����������
	double *mpMeasT;

	// �������� ����
	double *mpSinsDelayT;
	double *mparrArtParral; // ������ ���������� ��



	// 	��������� ��������� �������� ��������  �����
	double  *mpSigDrivAY_dU_po_dt ;



};
//---------------------------------------------------------------------------
extern PACKAGE TForm5 *Form5;
//---------------------------------------------------------------------------
#endif
