//---------------------------------------------------------------------------

#ifndef MainFormH
#define MainFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>

#include "Vessel.h"
//---------------------------------------------------------------------------
class TVessel;
class TForm1 : public TForm
{
__published:	// IDE-managed Components
	TButton *Button1;
private:	// User declarations
public:		// User declarations
	__fastcall TForm1(TComponent* Owner);

	TVessel mVessel;
};
//---------------------------------------------------------------------------
extern PACKAGE TForm1 *Form1;
//---------------------------------------------------------------------------
#endif
