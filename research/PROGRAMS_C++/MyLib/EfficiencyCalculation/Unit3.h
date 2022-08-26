//---------------------------------------------------------------------------

#ifndef Unit3H
#define Unit3H
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <Dialogs.hpp>
//---------------------------------------------------------------------------
class TTkachenkoForm : public TForm
{
__published:	// IDE-managed Components
	TLabel *Label1;
	TEdit *Edit1;
	TButton *Button1;
	TLabel *Label2;
	TComboBox *CB3;
	TButton *Button2;
	TComboBox *CB1;
	TLabel *Label3;
	TOpenDialog *OpenDlg;
	TComboBox *CB4;
	TLabel *Label4;
	TComboBox *CB2;
	TLabel *Label5;
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Button2Click(TObject *Sender);
private:	// User declarations
public:		// User declarations

  UnicodeString ResultDir;

	__fastcall TTkachenkoForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TTkachenkoForm *TkachenkoForm;
//---------------------------------------------------------------------------
#endif
