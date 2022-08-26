//---------------------------------------------------------------------------

#ifndef LineH
#define LineH
#include "URPolyLine.h"
class TURPointXY;
class TLine:public TURPolyLine
{
public:


	TLine();

	TLine::TLine(TURPointXY * pPoints);

	TLine::TLine(const TURPointXY  pnt1, const TURPointXY  pnt2);

	bool  IsVertical();

	bool TLine::calcTang(double *pvalTang);

	bool  calcY(const double x, double *py);

	bool  calcPointXY(const double x, TURPointXY *pPointXY) ;


};
#endif
