//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "Unit3.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TTkachenkoForm *TkachenkoForm;
//---------------------------------------------------------------------------
__fastcall TTkachenkoForm::TTkachenkoForm(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TTkachenkoForm::Button1Click(TObject *Sender)
{ AnsiString str;
  int i;

	OpenDlg->FileName= "[выбирете директорию]";

	if(mrOk==OpenDlg->Execute())
	{
		str= OpenDlg->FileName;
		for(i=str.Length()-1;i>0;i--)
		{
			if( str.c_str()[i]=='\\') break;

		}
		str.SetLength(i+1);
		ResultDir= str;
		Edit1->Text= ResultDir;
  }
}
//---------------------------------------------------------------------------
void __fastcall TTkachenkoForm::Button2Click(TObject *Sender)
{


  ModalResult= mrOk;
}
//---------------------------------------------------------------------------



