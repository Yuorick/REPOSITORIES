//---------------------------------------------------------------------------

#ifndef FarFormH
#define FarFormH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TForm2 : public TForm
{
__published:	// IDE-managed Components
	TButton *Button1;
	TPanel *Panel1;
	TImage *Image1;
	TLabel *Label1;
	void __fastcall FormShow(TObject *Sender);
	void __fastcall Image1MouseDown(TObject *Sender, TMouseButton Button, TShiftState Shift,
          int X, int Y);
	void __fastcall Button1Click(TObject *Sender);
	void __fastcall Image2MouseDown(TObject *Sender, TMouseButton Button, TShiftState Shift,
          int X, int Y);
private:	// User declarations
public:		// User declarations


// размеры решетки
  int mNumRows;
  int mNumCols;
  ///

//указатель на массив индикатор рабочих АМ  у модельной ФАР
  bool *mpbarrWorking;

  //указатель на массив индикатор рабочих АМ  у реальной ФАР
  bool *mpbarrRealWorking;

// признак 5П10
  bool *mpb5P10;

  // размеры окна для АМ
  int miWidth;
  int miHeight ;

   // размер АМ на картинке
  double mdx ;
  double mdy;
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
 // int nX,nY;	// размеры решетки
//  bool* fl[100];     // массив исправности решетки

//  int W,H;
 // double dx,dy;


	__fastcall TForm2(TComponent* Owner);
	void __fastcall DrawCell(int x,int y) ;
	void __fastcall DrawCell_(int i,int j);
	void __fastcall drawFar();
	void __fastcall drawRealFar();

	void __fastcall _DrawCell_(TImage *Image1, bool *pbarrWorking, int i,int j);

};
//---------------------------------------------------------------------------
extern PACKAGE TForm2 *Form2;
//---------------------------------------------------------------------------
#endif
