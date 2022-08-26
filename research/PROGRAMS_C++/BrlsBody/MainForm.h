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
//---------------------------------------------------------------------------
class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TOpenDialog *OpenDialog1;
	TOpenDialog *OpenDialog2;
	TOpenDialog *OpenDialog3;
	TOpenDialog *OpenDialog4;
	TOpenDialog *OpenDialog5;
	TOpenDialog *OpenDialog6;
	TOpenDialog *OpenDialog7;
	TOpenDialog *OpenDialog8;
	TOpenDialog *OpenDialog9;
	TOpenDialog *OpenDialog10;
	TPanel *Panel3;
	TPanel *Panel9;
	TLabel *Label11;
	TLabel *Label1;
	TLabel *Label2;
	TButton *Button5;
	TEdit *Edit6;
	TButton *Button2;
	TButton *Button3;
	TEdit *Edit1;
	TButton *Button4;
	TEdit *Edit2;
	TLabel *Label3;
	TPanel *Panel1;
	TLabel *Label4;
	TLabel *Label5;
	TLabel *Label6;
	TLabel *Label7;
	TLabel *Label8;
	TLabel *Label9;
	TButton *Button6;
	TEdit *Edit3;
	TButton *Button7;
	TButton *Button8;
	TEdit *Edit4;
	TButton *Button9;
	TEdit *Edit5;
	TButton *Button10;
	TEdit *Edit7;
	TEdit *Edit8;
	TButton *Button11;
	TPanel *Panel2;
	TLabel *Label10;
	TLabel *Label12;
	TLabel *Label14;
	TButton *Button12;
	TEdit *Edit9;
	TButton *Button13;
	TButton *Button14;
	TEdit *Edit10;
	TLabel *Label13;
	TPanel *Panel4;
	TButton *Button1;
	TButton *Button15;
	TEdit *Edit11;
	TLabel *Label15;
	TOpenDialog *OpenDialog11;
	TSaveDialog *SaveDialog1;
	TLabel *Label16;
	TEdit *Edit12;
	TButton *Button16;
	TButton *Button17;
	TEdit *Edit13;
	TLabel *Label17;
	TButton *Button18;
	TEdit *Edit14;
	TLabel *Label18;
	TOpenDialog *OpenDialog12;
	TOpenDialog *OpenDialog13;
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Button5Click(TObject *Sender);
	void __fastcall Button2Click(TObject *Sender);
	void __fastcall Button3Click(TObject *Sender);
	void __fastcall Button4Click(TObject *Sender);
	void __fastcall Button6Click(TObject *Sender);
	void __fastcall Button8Click(TObject *Sender);
	void __fastcall Button9Click(TObject *Sender);
	void __fastcall Button10Click(TObject *Sender);
	void __fastcall Button7Click(TObject *Sender);
	void __fastcall Button11Click(TObject *Sender);
	void __fastcall Button12Click(TObject *Sender);
	void __fastcall Button14Click(TObject *Sender);
	void __fastcall Button13Click(TObject *Sender);
	void __fastcall Button15Click(TObject *Sender);
	void __fastcall Button16Click(TObject *Sender);
	void __fastcall Button17Click(TObject *Sender);
	void __fastcall Button18Click(TObject *Sender);
private:	// User declarations
public:		// User declarations
	__fastcall TForm1(TComponent* Owner);
	wchar_t *mpwchSHP_HorProj;
	wchar_t *mpwchSHP_Wing;
	wchar_t *mpwchSHP_Stab;

	wchar_t *mpwchSHP_VertProj;
	wchar_t *mpwchSHP_Rule;
	wchar_t *mpwchSHP_Shaft;
	wchar_t *mpwchSHP_Sharnirs;
	wchar_t *mpwchSHP_SGF_Line;
	wchar_t *mpwchSHP_Wing_Line;
	wchar_t *mpwchSHP_Stab_Line;

	wchar_t *mpwchSHP_FaceProj;
	wchar_t *mpwchSHP_Horiz_Line;


	wchar_t *mpwchFileInp;
	wchar_t  mwchFoldInp[300];
	wchar_t *mpwchFileOut;
	wchar_t  mwchFoldOut[300];
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
