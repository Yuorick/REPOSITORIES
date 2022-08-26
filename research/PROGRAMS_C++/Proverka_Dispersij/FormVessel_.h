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
	TEtalonSign mEtalonSign;
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
	// �������� ������
	TAM_2D mAM_2D;

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

	// ��� ���������� � �������
	TFar_2D mFar_2D;

	// �������� �� ��������
	double mPowerPrd;
	// �� �� ��������
	double mKYPrd;
	// ������� �� ��������
	TTransmitAnt mTransmitAnt ;
	// ���������� ������� �������� � ��������� ��� (������� �������� � ��������)
	bool mbarrAM[MAX_QUANT_AM];
  //-------------------------------------------------------------------
	// ��� ������� ������ ���� (����� �����)
	double	mSigSins;

	// ��� ������� ������ ������ �������� ����� (����� �����)
	double	mSig_d_po_dt_Sins;

	// ������� ���
	double mVesselWidth ; // ������(�)
	double mVesselLength ;
	double marrFarParallacs[3] ;

	double mMaxQ ; /// ������������ ���� ���������� �� ������������ �����(��������� ���� ��������)
	double mT_Q; // ������ ��������
	double mMaxPsi ;// ������������ ���� ������� �����(���������)
	double mT_Psi ; // ������ ������� �����
	double mMaxTet ; //������������ ���� ��������� �����(���������)
	double mT_Tet; // ������ �������� �����
	double mMaxVert  ;

	// ���������� ��������  ������� ������
	double mQ0  ; // ����������� ����
	double mVVess  ;// �������� ������� ������ 20 �����
	double marrDelt[4] ;//  ��������� ����

	//������������ ��������� ��������� ������ ������� � ��� �� 100 �
	double mMaxAmp_AftFlexure;
	// ������ ��������� ��������� ������
	double mT_AftFlexure;
	//������������ ��������� ��������� ������ ������� � ��� �� 100 �
	double mMaxAmp_BoardFlexure;
	// ������ ��������� ��������� ������
	double mT_BoardFlexure;

	// 3.1 �������� ����
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
	// ��� ��������� ����������
	//	bool mbSkaliga;

	// ������
	double mDriverSigBet ;// �������� ��������� ���� Bet �������
	double mDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
	double mDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
	double mDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����



	//

	//���� ����������
	// �������� ����� �����������
	double mMeasT;

	// �������� ����
	double mSinsDelayT;


	TBius mBius;



	double marrArtParral[3]; // ������ ���������� ��



	// 	��������� ��������� �������� ��������  �����
	double  mSigDrivAY_dU_po_dt ;

	void __fastcall fncInputData();

	void __fastcall create5P10();


};
//---------------------------------------------------------------------------
extern PACKAGE TForm5 *Form5;
//---------------------------------------------------------------------------
#endif
