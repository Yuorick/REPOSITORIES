//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "MainForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
#include "YrRastr.h"


#include "YrWrite.h"
#include "YrRead.h"
TForm3 *Form3;
class TYrRead;
//---------------------------------------------------------------------------
__fastcall TForm3::TForm3(TComponent* Owner)
	: TForm(Owner)
{
}

//---------------------------------------------------------------------------

void __fastcall TForm3::Button2Click(TObject *Sender)
{
  // Edit1->Text = L"-999";
	 // Edit4->Text = L"-999";
		OpenDialog1->Filter = L"???? c ??? (*.csv)|*.csv" ;
	if(OpenDialog1->Execute())
		Edit1->Text = OpenDialog1->FileName ;
		// sFileN = OpenDialog1->FileName ;
		// ShowMessage(sFileN);
		 pchInpFileNameDEM = (OpenDialog1->FileName).w_str() ;
}
//---------------------------------------------------------------------------

void __fastcall TForm3::Button3Click(TObject *Sender)
{

   double arrMass[NN] ;
   int  quantPoints = -1;
   mNum?ols = StrToInt(LabeledEdit19->Text);

   int nrows = -1;

   TYrRead::YrReadCSV_(pchInpFileNameDEM ,  mNum?ols,&nrows,arrMass);

	switch(ComboBox4->ItemIndex)   // ????
	{
	case 0:
	TYrWrite::WriteReportForIntMassiveTXT_Comma(pchInpFileNameDigPoints_Rez,arrMass,mNum?ols, nrows);
   //	mWSkz = 0.002;
	break;
	case 1:
	TYrWrite::WriteReportForIntMassiveTXT_Semicolon(pchInpFileNameDigPoints_Rez,arrMass,mNum?ols, nrows);
	break ;



	default:
	break;

	}
 


  int a =1;

}
//---------------------------------------------------------------------------










void __fastcall TForm3::Button8Click(TObject *Sender)
{
	 Application->Terminate();
}
//---------------------------------------------------------------------------



void __fastcall TForm3::Button7Click(TObject *Sender)
{


   SaveDialog1->Filter = L"TXT ????? (*.txt)|*.txt" ;


  if (SaveDialog1->Execute())
  {
   //	ShowMessage( (SaveDialog1->FileName).w_str()) ;
   pchInpFileNameDigPoints_Rez =  (SaveDialog1->FileName).w_str();

  }


int iii = 0;

}
//---------------------------------------------------------------------------






