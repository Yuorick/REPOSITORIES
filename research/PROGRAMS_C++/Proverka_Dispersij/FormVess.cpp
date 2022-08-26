//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
 #include "MainForm.h"
#include "FormVess.h"
#include "YrString.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm5 *Form5;
//---------------------------------------------------------------------------
__fastcall TForm5::TForm5(TComponent* Owner)
	: TForm(Owner)
{
   int iii = 0;
}
//---------------------------------------------------------------------------
void __fastcall TForm5::LabeledEdit26Change(TObject *Sender)
{
double t = StrTo_Dbl_(LabeledEdit26->Text) ;
*mpVesselWidth =  t;
}
//---------------------------------------------------------------------------

void __fastcall TForm5::Button1Click(TObject *Sender)
{
double t = StrTo_Dbl_(LabeledEdit26->Text) ;
*mpVesselWidth =  t;
Close();
}
//---------------------------------------------------------------------------


