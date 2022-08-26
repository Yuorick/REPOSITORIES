//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "FarForm.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm2 *Form2;
//---------------------------------------------------------------------------
__fastcall TForm2::TForm2(TComponent* Owner)
	: TForm(Owner)
{
  ;
}
//---------------------------------------------------------------------------
void __fastcall TForm2::FormShow(TObject *Sender)
{
  drawFar();
 // drawRealFar();
 
}

void __fastcall TForm2::drawFar()
{

      Image1->Canvas->Brush->Color=clSilver;//clRed;
	  Image1->Canvas->Pen->Width=1;
	  Image1->Canvas->Pen->Color= clBlack;//clWhite;


	  miWidth = Image1->Width;
	  miHeight =Image1->Height;

	  mdx=(double)miWidth/mNumCols;
	  mdy=(double)miHeight/mNumRows;

	  TRect rec(0,miHeight,miWidth,0);

	  Image1->Canvas->FillRect(rec);

	//  TRect rec1(0,mNumRows *miHeight, mNumCols *miWidth, 0);

 /* if(!(*mpb5P10))
  {




  }
  else
  {

  }*/
  for( int i =0; i < mNumRows; i++)
	  for(int j=0; j < mNumCols; j++)
	  {
		  DrawCell_(i, j);
	  }
	  Image1->Canvas->Pen->Width= 5;
	  Image1->Canvas->MoveTo(rec.Left,rec.Top);
	  Image1->Canvas->LineTo(rec.Right,rec.Top);
	  Image1->Canvas->LineTo(rec.Right,rec.Bottom);
	  Image1->Canvas->LineTo(rec.Left,rec.Bottom);
	  Image1->Canvas->LineTo(rec.Left,rec.Top);

	  ///
  	  Image1->Canvas->Pen->Width= 7;
	  Image1->Canvas->MoveTo(rec.Left,rec.Top/2);
	  Image1->Canvas->LineTo(rec.Right,rec.Top/2);
	  Image1->Canvas->MoveTo(rec.Right/2,rec.Top);
	  Image1->Canvas->LineTo(rec.Right/2,rec.Bottom);
	  Image1->Canvas->Pen->Width= 1;
	 // memcpy( mpbarrRealWorking, mpbarrWorking, sizeof(bool) * mNumCols * mNumRows);

}


 /*
void __fastcall TForm2::drawRealFar()
{

	  Image2->Canvas->Brush->Color=clSilver;//clRed;
	  Image2->Canvas->Pen->Width=1;
	  Image2->Canvas->Pen->Color= clBlack;//clWhite;


	  miWidth = Image2->Width;
	  miHeight =Image2->Height;

	  mdx=(double)miWidth/mNumCols;
	  mdy=(double)miHeight/mNumRows;

	  TRect rec(0,miHeight,miWidth,0);

	  Image2->Canvas->FillRect(rec);



  if(!(*mpb5P10))
  {


	  for( int i =0; i < mNumRows; i++)
	  for(int j=0; j < mNumCols; j++)
	  {
		  _DrawCell_(Image2, mpbarrRealWorking, i, j);
	  }
	  Image2->Canvas->Pen->Width= 5;
	  Image2->Canvas->MoveTo(rec.Left,rec.Top);
	  Image2->Canvas->LineTo(rec.Right,rec.Top);
	  Image2->Canvas->LineTo(rec.Right,rec.Bottom);
	  Image2->Canvas->LineTo(rec.Left,rec.Bottom);
	  Image2->Canvas->LineTo(rec.Left,rec.Top);


  }
  else
  {
	  for(int j=1; j < (mNumCols -1); j++)
	  {
		  _DrawCell_(Image2,mpbarrRealWorking,0 , j);
	  }
	  for( int i =1; i < 3; i++)
	  for(int j=0; j < mNumCols; j++)
	  {
		  _DrawCell_(Image2, mpbarrRealWorking, i, j);
	  }
	  for(int j= 1; j < (mNumCols -1); j++)
	  {
		  _DrawCell_(Image2,mpbarrRealWorking, 3 , j);
	  }

	  Image2->Canvas->Pen->Width= 5;
	  Image2->Canvas->MoveTo(0,mdy);
	  Image2->Canvas->LineTo(0,3 * mdy);
	  Image2->Canvas->LineTo(mdx,3 * mdy);
	  Image2->Canvas->LineTo(mdx,4* mdy);
	  Image2->Canvas->LineTo(7 *mdx,4* mdy);
	  Image2->Canvas->LineTo(7 *mdx,3* mdy);
	  Image2->Canvas->LineTo(8 *mdx,3* mdy);
	  Image2->Canvas->LineTo(8 *mdx, mdy);
	  Image2->Canvas->LineTo(7 *mdx, mdy);
	  Image2->Canvas->LineTo(7 *mdx, 0);
	  Image2->Canvas->LineTo(mdx, 0);
	  Image2->Canvas->LineTo(mdx, mdy);
	  Image2->Canvas->LineTo(0, mdy);
	  Image2->Canvas->LineTo(0,3 * mdy);


  }
  	  Image2->Canvas->Pen->Width= 7;
	  Image2->Canvas->MoveTo(rec.Left,rec.Top/2);
	  Image2->Canvas->LineTo(rec.Right,rec.Top/2);
	  Image2->Canvas->MoveTo(rec.Right/2,rec.Top);
	  Image2->Canvas->LineTo(rec.Right/2,rec.Bottom);
	  Image2->Canvas->Pen->Width= 1;


}
*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
/*
void __fastcall TForm2::DrawCell(int x,int y)
{ TRect r;

  r.Left  = (int)(dx*x);
  r.Top   = (int)(dy*y);
  r.Right = (int)(dx*(x+1));
  r.Bottom= (int)(dy*(y+1));
  Image1->Canvas->FillRect(r);

  if(fl[x][y])
  { Image1->Canvas->MoveTo(r.Left,r.Top);
	Image1->Canvas->LineTo(r.Right,r.Bottom);
		Image1->Canvas->MoveTo(r.Left,r.Bottom);
	Image1->Canvas->LineTo(r.Right,r.Top);  }

	Image1->Canvas->MoveTo((int)(dx*x),(int)(dy*y));
  Image1->Canvas->LineTo((int)(dx*x),(int)(dy*(y+1)));
	Image1->Canvas->MoveTo((int)(dx*(x+1)),(int)(dy*y));
  Image1->Canvas->LineTo((int)(dx*(x+1)),(int)(dy*(y+1)));

	Image1->Canvas->MoveTo((int)(dx*x),(int)(dy*y));
  Image1->Canvas->LineTo((int)(dx*(x+1)),(int)(dy*y));
	Image1->Canvas->MoveTo((int)(dx*x),(int)(dy*(y+1)));
  Image1->Canvas->LineTo((int)(dx*(x+1)),(int)(dy*(y+1)));

}
 */
//---------------------------------------------------------------------------
void __fastcall TForm2::DrawCell_(int i,int j)
{ TRect r;

  r.Left  = (int)(mdx * j);
 // r.Top   = (int)(mdy * j);
 r.Top   = (int)(mdy*(i+1));
  r.Right = (int)(mdx *(j+1));
 // r.Bottom= (int)(mdy*(j+1));
 r.Bottom= (int)(mdy * i);
  Image1->Canvas->FillRect(r);

  if(!mpbarrWorking[ i * mNumCols + j])
  { Image1->Canvas->MoveTo(r.Left,r.Top);
	Image1->Canvas->LineTo(r.Right,r.Bottom);
	 	Image1->Canvas->MoveTo(r.Left,r.Bottom);
	Image1->Canvas->LineTo(r.Right,r.Top);
  }

	Image1->Canvas->MoveTo((int)(mdx * j),(int)(mdy * i));
  Image1->Canvas->LineTo((int)(mdx * j),(int)(mdy *(i+1)));
	Image1->Canvas->MoveTo((int)(mdx *(j+1)),(int)(mdy*i));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*(i+1)));

	Image1->Canvas->MoveTo((int)(mdx*j),(int)(mdy*i));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*i));
	Image1->Canvas->MoveTo((int)(mdx*j),(int)(mdy*(i+1)));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*(i+1)));

}

void __fastcall TForm2::_DrawCell_(TImage *Image1, bool *pbarrWorking, int i,int j)
{ TRect r;

  r.Left  = (int)(mdx * j);
 // r.Top   = (int)(mdy * j);
 r.Top   = (int)(mdy*(i+1));
  r.Right = (int)(mdx *(j+1));
 // r.Bottom= (int)(mdy*(j+1));
 r.Bottom= (int)(mdy * i);
  Image1->Canvas->FillRect(r);

  if(!pbarrWorking[ i * mNumCols + j])
  { Image1->Canvas->MoveTo(r.Left,r.Top);
	Image1->Canvas->LineTo(r.Right,r.Bottom);
	 	Image1->Canvas->MoveTo(r.Left,r.Bottom);
	Image1->Canvas->LineTo(r.Right,r.Top);
  }

	Image1->Canvas->MoveTo((int)(mdx * j),(int)(mdy * i));
  Image1->Canvas->LineTo((int)(mdx * j),(int)(mdy *(i+1)));
	Image1->Canvas->MoveTo((int)(mdx *(j+1)),(int)(mdy*i));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*(i+1)));

	Image1->Canvas->MoveTo((int)(mdx*j),(int)(mdy*i));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*i));
	Image1->Canvas->MoveTo((int)(mdx*j),(int)(mdy*(i+1)));
  Image1->Canvas->LineTo((int)(mdx*(j+1)),(int)(mdy*(i+1)));

}
/*
void __fastcall TForm2::Image1MouseDown(TObject *Sender, TMouseButton Button, TShiftState Shift,
          int X, int Y)
{
 int x,y;

  for(x=0;x<nX;x++)
  { if((X>(int)(dx*x))&&(X<=(int)(dx*(x+1)))) break;
  }
  for(y=0;y<nY;y++)
  { if((Y>(int)(dy*y))&&(Y<=(int)(dy*(y+1)))) break;
  }

  fl[x][y]=!fl[x][y];
  DrawCell(x,y);

}
*/
void __fastcall TForm2::Image1MouseDown(TObject *Sender, TMouseButton Button, TShiftState Shift,
          int X, int Y)
{
 int numRow,numCol;

  for(numRow = 0;numRow < mNumRows; numRow++)
  {
   if((Y>(int)(mdy * numRow))&&(Y<=(int)(mdy*(numRow +1 )))) break;
  }
  for(numCol =0; numCol < mNumCols ;numCol++)
  {
   if(( X > (int)(mdx * numCol))&&( X <= (int)(mdx*(numCol + 1)))) break;
  }

/*  if(*mpb5P10)
  {
  if( ((numRow == 0)&& (numCol == 0))
	||((numRow == 0)&& (numCol == 7))
	||((numRow == 3)&& (numCol == 0))
	||((numRow == 3)&& (numCol == 7))
	)
	{
	  return;
	}
  }*/
  mpbarrWorking[numRow * mNumCols +  numCol] = !mpbarrWorking[numRow * mNumCols +  numCol];
  drawFar();
//  memcpy( mpbarrRealWorking, mpbarrWorking, sizeof(bool) * mNumCols * mNumRows);
 // drawRealFar();
 // DrawCell_(numRow,numCol);

}
//---------------------------------------------------------------------------
void __fastcall TForm2::Button1Click(TObject *Sender)
{
Close();
}
//---------------------------------------------------------------------------


void __fastcall TForm2::Image2MouseDown(TObject *Sender, TMouseButton Button, TShiftState Shift,
          int X, int Y)
{
int numRow,numCol;

  for(numRow = 0;numRow < mNumRows; numRow++)
  {
   if((Y>(int)(mdy * numRow))&&(Y<=(int)(mdy*(numRow +1 )))) break;
  }
  for(numCol =0; numCol < mNumCols ;numCol++)
  {
   if(( X > (int)(mdx * numCol))&&( X <= (int)(mdx*(numCol + 1)))) break;
  }

  if(*mpb5P10)
  {
  if( ((numRow == 0)&& (numCol == 0))
	||((numRow == 0)&& (numCol == 7))
	||((numRow == 3)&& (numCol == 0))
	||((numRow == 3)&& (numCol == 7))
	)
	{
	  return;
	}
  }
 // mpbarrRealWorking[numRow * mNumCols +  numCol] = !mpbarrRealWorking[numRow * mNumCols +  numCol];
 // drawRealFar();
}
//---------------------------------------------------------------------------


