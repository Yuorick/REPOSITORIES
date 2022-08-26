//---------------------------------------------------------------------------

#ifndef FormSeaTargsH
#define FormSeaTargsH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Dialogs.hpp>
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
#include "ControlSyst.h"
#include "ArtCannon.h"
#include "ArtComplex.h"
#include "Vessel.h"
#include "Fight.h"
#include "ArtCannon.h"
#include <ButtonGroup.hpp>



//#define MAX_QUANT_AM 1600
// enum enumTargetType;
// enum enumShellType;
// enum enumCannonType;
//---------------------------------------------------------------------------
class TAM_2D;
class TFar_2D;
class TURPointZ;
class  TSingleSign;
class TInitTargData;
class TEtalonSign;
class TTransmitAnt;
class TWind;
class TEnvironment;
class TControlSyst;
class TArtCannon;
class TArtComplex;
class TVessel;
class TFight;
//----------------------------------------------------------------------------------------
class TForm4 : public TForm
{
__published:	// IDE-managed Components
	TPanel *Panel3;
	TLabel *Label9;
	TLabel *Label7;
	TLabel *Label4;
	TLabel *Label5;
	TComboBox *ComboBox3;
	TComboBox *ComboBox1;
	TComboBox *ComboBox2;
	TOpenDialog *OpenDialog1;
	TPanel *Panel2;
	TLabel *Label2;
	TLabeledEdit *LabeledEdit46;
	TLabeledEdit *LabeledEdit34;
	TPanel *Panel10;
	TLabel *Label1;
	TLabeledEdit *LabeledEdit49;
	TLabeledEdit *LabeledEdit51;
	TLabeledEdit *LabeledEdit52;
	TLabeledEdit *LabeledEdit53;
	TLabeledEdit *LabeledEdit50;
	TLabeledEdit *LabeledEdit2;
	TPanel *Panel7;
	TLabel *Label8;
	TLabeledEdit *LabeledEdit25;
	TLabeledEdit *LabeledEdit27;
	TLabeledEdit *LabeledEdit29;
	TLabeledEdit *LabeledEdit23;
	TLabeledEdit *LabeledEdit28;
	TButton *Button1;
	TLabeledEdit *LabeledEdit10;
	TComboBox *ComboBox5;
	TLabel *Label6;
	TPanel *Panel5;
	TSaveDialog *SaveDialog1;
	TButton *Button2;
	TButton *Button5;
	TEdit *Edit1;
	TLabel *Label11;
	TPanel *Panel1;
	TLabel *Label3;
	TButton *Button4;
	TEdit *Edit2;
	TPanel *Panel4;
	TButton *Button3;
	TButton *Button6;
	void __fastcall Button4Click(TObject *Sender);
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Button3Click(TObject *Sender);
	void __fastcall Button2Click(TObject *Sender);

	void __fastcall Button5Click(TObject *Sender);
	void __fastcall Button6Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TForm4(TComponent* Owner);
		// ���� � ����� � ���������
	wchar_t *mpwchOutFile0;
	// ���� � ����� � ���������  ����������� ������� ��������
	wchar_t mwchOutFold[400]; //

		// ���� � � ����� � �� �� ������
   wchar_t *mpwchOutDataFileTaran;
	// ������

	//��� ����������
	enumDetonatorType mDetonatorType;
	enumShellType mEnumShellType;
	enumCannonType menumCannonType;
	//
	//�-�� ��������
	int mQuantShells;
	
	///

	// ����

	//��� ����
	enumTargetType  mEnumTargType;
	//���� �������, ����
	double mBearing0;
	//������, �
	double mElev0;
	//��������, �/�
	double mVelocity0;
	//���������, �
	double mDist0;

	//��� ���� ��������, �/�/�
 	double mWSkz;
	//
	double mTargZenitAng0;
	//���� �����
	double mTargCourse0;
	// ��� ����
	double mTargEPR ;

	// �������� � ������ ��
	double mSigAUDelayT;

	// �����������
	// ���� ��������
	double mRateOfFire;

	// 	��������� ��������� ������� �����
	double  mSigDrivAY_U ;

	// 	��������� ��������� �������� ��������  �����
	double  mSigDrivAY_dU_po_dt ;

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


	// ������
	double mDriverSigBet ;// �������� ��������� ���� Bet �������
	double mDriverSigEps ;// �������� ��������� ���� Eps  ������� (���� �����)
	double mDriverDynamicSigBet ;// �������� ��������� ���� �����  �������
	double mDriverDynamicSigEps ;// ��������  ������� ��������� ���� �����

	//
	 double marrArtParral[3]; // ������ ���������� ��
	//���� ����������
	// �������� ����� �����������
	double mMeasT;

	// �������� ����
	double mSinsDelayT;
  // ���� ������� ���
  double mRzvT;

	TControlSyst mControlSyst;
	TFar_2D mFar_2D;
	TArtCannon mArtCannon ;
	TArtComplex mArtComplex;  // ��
	TTransmitAnt mTransmitAnt;
	TEtalonSign mEtalonSign;
	TVessel mVessel;
	TFight mFight;
	TEnvironment mEnvironment;
	TInitTargData mInitTargData;
	TShellBody mShellBody;

	void __fastcall fncInputData();
	void __fastcall create5P10();
};
//---------------------------------------------------------------------------
extern PACKAGE TForm4 *Form4;
//---------------------------------------------------------------------------
#endif