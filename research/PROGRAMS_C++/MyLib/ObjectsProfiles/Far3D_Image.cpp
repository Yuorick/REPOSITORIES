//---------------------------------------------------------------------------


#pragma hdrstop
#include <math.h>
 #include <string.h>
 #include <stdlib.h>
#include "URPolygon.h"
#include "UrPointXY.h"
#include "URPolyLine.h"
#include "Far3D_Image.h"

// INPUT:
// wchFoldName - путь к папке с файлами
// Emitter:
// ValLength0 - длина горизонтальной части эмиттера
// ValLength1 - длина наклонной части эмиттера
// ValFi - угол между горизонтальной и наклонной частями эмиттера
//  AM:
// QuantEmit - к-во эмиттеров
// DistEmit - расстояние мжде эмиттерами
// ValLengthBoardAM - расстояние между границами АМ
// SideWidthAM - толщина боковых частей полигона кожуха эмиттеров
// BackWidthAM - толщина задней части  полигона кожуха эмиттеров
// LenghPlnOutAM  - длина линии выхода с АМ
// FAR:
// QuantAM - к-во АМ
//
//

void createFar3DPict(wchar_t *wchFoldName, const int QuantRows, const int QuantCols,const double VALx,const double VALy
 ,const double VAL_al1  , const double VALh,const double VAL_al2)
{
  const int QuantParall =  QuantRows *   QuantCols + 2;
  TURPolygon *parrPLg = new TURPolygon[ QuantParall];
  TURPointXY *ppntarrLeftLow  = new TURPointXY  [ QuantParall];
  TURPointXY *ppntarrPhCentrs  = new TURPointXY  [ QuantRows *   QuantCols];

  // начальная точка паралеллограммов - левая нижняя
  TURPointXY pntStart(-((double)(QuantCols/2))* VALx * cos(VAL_al1)
  , -((double)(QuantCols/2))* VALx * sin(VAL_al1) - ((double)(QuantRows/2)) *VALy);
  for (int i =0; i < QuantRows; i++)
  {
  TURPointXY pntStartRow(pntStart.X,pntStart.Y + ((double)i) * VALy);
  for (int j =0; j < QuantCols; j++)
  {
	ppntarrLeftLow[i * QuantCols + j] = TURPointXY( pntStartRow.X +  ((double)j)* VALx * cos(VAL_al1)
	,pntStartRow.Y +  ((double)j)* VALx * sin(VAL_al1));
	ppntarrPhCentrs[i * QuantCols + j].X = ppntarrLeftLow[i * QuantCols + j].X + VALx * cos (VAL_al1)/ 2.;
	ppntarrPhCentrs[i * QuantCols + j].Y = ppntarrLeftLow[i * QuantCols + j].Y + VALx * sin (VAL_al1)/ 2. + VALy /2.;


  }
  }

  ppntarrLeftLow[QuantParall-2] = pntStart;
  ppntarrLeftLow[QuantParall-1] = TURPointXY(pntStart.X, pntStart.Y + ((double)QuantRows) * VALy); ;

  for (int i =0; i < QuantRows; i++)
  for (int j =0; j < QuantCols; j++)
  {
	parrPLg [i * QuantCols + j] = TURPolygon::fncCreateParallelogramm(ppntarrLeftLow[i * QuantCols + j]
	  , VALx, VALy ,VAL_al1 );
  }

   parrPLg[QuantParall-2] = TURPolygon(5);
   parrPLg[QuantParall-1] = TURPolygon(5);
   parrPLg[QuantParall-2].Points [0] = pntStart;
   parrPLg[QuantParall-2].Points [4] = pntStart;
   parrPLg[QuantParall-2].Points [1] = TURPointXY(pntStart.X - VALh * cos(VAL_al2),pntStart.Y + VALh * sin(VAL_al2)) ;
   parrPLg[QuantParall-2].Points [2] = parrPLg[QuantParall-2].Points [1];
   parrPLg[QuantParall-2].Points [2].Y +=  ((double)QuantRows) *  VALy;

   parrPLg[QuantParall-2].Points [3] = TURPointXY(pntStart.X, pntStart.Y + ((double)QuantRows) *  VALy);

   parrPLg[QuantParall-1].Points [0] =parrPLg[QuantParall-2].Points [3];
   parrPLg[QuantParall-1].Points [4] =parrPLg[QuantParall-2].Points [3];
   parrPLg[QuantParall-1].Points [1] =parrPLg[QuantParall-2].Points [2];
   parrPLg[QuantParall-1].Points [3].X  =parrPLg[QuantParall-1].Points [0].X + ((double)QuantCols) * VALx * cos(VAL_al1);
   parrPLg[QuantParall-1].Points [3].Y  =parrPLg[QuantParall-1].Points [0].Y + ((double)QuantCols) * VALx * sin(VAL_al1);

   parrPLg[QuantParall-1].Points [2].X = parrPLg[QuantParall-1].Points [3].X
		   + parrPLg[QuantParall-1].Points [1].X-parrPLg[QuantParall-1].Points [0].X;
   parrPLg[QuantParall-1].Points [2].Y = parrPLg[QuantParall-1].Points [3].Y
		   + parrPLg[QuantParall-1].Points [1].Y - parrPLg[QuantParall-1].Points [0].Y;

	 wchar_t wchFileName[300]= {0};
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\FAR.shp");
   TURPolygon::WriteSetSHPFiles(wchFileName,parrPLg, QuantParall ) ;


	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\PhCentres.shp");
	TURPointXY::WriteSetSHPFiles(wchFileName,ppntarrPhCentrs, QuantRows *   QuantCols ) ;
	delete []parrPLg;
	delete []ppntarrLeftLow;
	delete [] ppntarrPhCentrs ;

}
#pragma package(smart_init)
