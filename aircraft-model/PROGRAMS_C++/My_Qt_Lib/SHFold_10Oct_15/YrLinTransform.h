//---------------------------------------------------------------------------

#ifndef YrLinTransformH
#define YrLinTransformH
//---------------------------------------------------------------------------


void LinTransf(const double x0, const double y0,
				const double alf0,
				const double x, const double y
				, double *xTr,  double *yTr
				);
double dist2(const double x1,const double y1,const double x2,const double y2);
double AngleBetweenVectors(const double x1,const double y1,const double x2,const double y2);
#endif
