//---------------------------------------------------------------------------

#ifndef PicturesForReport1H
#define PicturesForReport1H
class TURPointXY;
void  _fastcall createPicture1(wchar_t *wchFoldName ) ;
bool findMirrowPnt(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, const double valx0, TURPointXY *pntMirrow);

void   createSetGraphs(wchar_t *wchFoldName,const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY pntMirrow, TURPointXY &pntAntipod);

void	 findAllMirrPnts(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY *arrPntMirr, int &quantMirrPoints);

void   calcAntipod(const double valAmpl ,const double  valOmega
	, const double  valPh0,const TURPointXY pntRadar, const TURPointXY pntTagr
	, TURPointXY pntMirrow, TURPointXY *pntAntipod);

void  _fastcall createPicture1(wchar_t *wchFoldName, double valVessRast );
#endif
