//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
//---------------------------------------------------------------------------
#include "URPointZ.h"
#include "Environment.h"
#include "MyShellTraj.h"
#include "LearnShellBody.h"
//#include "URPolyLine.h"
#include <Dialogs.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
 class URPolyLine;
 class TFar_2D;
class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TOpenDialog *OpenDialog1;
	TPanel *Panel9;
	TOpenDialog *OpenDialog2;
	TOpenDialog *OpenDialog3;
	TComboBox *ComboBox4;
	TLabel *Label14;
	TLabeledEdit *LabeledEdit16;
	TPanel *Panel3;
	TButton *Button1;
	TButton *Button5;
	TEdit *Edit2;
	TLabel *Label5;
	TPanel *Panel5;
	TLabel *Label6;
	TLabeledEdit *LabeledEdit7;
	TLabeledEdit *LabeledEdit6;
	TPanel *Panel6;
	TLabel *Label7;
	TLabeledEdit *LabeledEdit9;
	TLabeledEdit *LabeledEdit11;
	TPanel *Panel1;
	TPanel *Panel7;
	TLabel *Label8;
	TLabeledEdit *LabeledEdit17;
	TLabeledEdit *LabeledEdit19;
	TLabeledEdit *LabeledEdit18;
	TButton *Button2;
	TLabeledEdit *LabeledEdit1;
	TLabeledEdit *LabeledEdit2;







	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Button5Click(TObject *Sender);
	void __fastcall Button2Click(TObject *Sender);

	
   //	void __fastcall Button7Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TForm1(TComponent* Owner);


	// ���� � ����� � ��������� ����� �������
	wchar_t *mpwchOutFileAppointmentPoints;


	// ���� � ����� � ���������
	wchar_t *mpwchOutFile0;
	// ���� � ����� � ���������  ����������� ������� ��������
	wchar_t mwchOutFold[400]; //



		// ���� � ����� � ���������  ���������� � ���������
	wchar_t *mpwchOutFileTraj0;
	wchar_t mpwchOutFileTraj[400];

   // csv ���� � ��������
   wchar_t mpwchTableFile[400];
   //	wchar_t *mpwchTableFile0;

	// ���������� � ���������
	wchar_t mpwchGraphDir[400];
	wchar_t *mpwchGraphDir0;

	// ���������� � ���������  ���������
	wchar_t mpwchDispGraphDir[400];
	wchar_t *mpwchDispGraphDir0;
	//���� ������������ � ���
	double mEpsGSK;
	double mBetGSK;

	// �������� �������
	double mVVess;

	// �������� �������
	double mVShell;

	//


	void __fastcall fncInputData();

};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
