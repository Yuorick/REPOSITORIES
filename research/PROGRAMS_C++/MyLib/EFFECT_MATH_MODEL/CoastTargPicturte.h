//---------------------------------------------------------------------------

#ifndef CoastTargPicturteH
#define CoastTargPicturteH
class TURPolygon;
class TURFigure;
void createNibourAppointmPointPictureForCoastTarg(wchar_t *wchOutFold, double *arrShellScatteringsCorMtrxPos_SSK
		 , TURPolygon Polygon, const double VAlTargCourse0, const double VAlPsi);

void createNibourAppointmPointPictureForCoastTarg(wchar_t *wchOutFold, double *arrElK
		 , TURFigure *pFigure, const double VAlTargCourse0, const double VAlPsi);

void createNibourAppointmPointPictureForCoastTarg(wchar_t *wchOutFold, double *arrMtrxCorrFluct
		 , double *arrMtrxCorrSyst, TURFigure *pFigure, const double VAlTargCourse0, const double VAlPsi);

void createNibourAppointmPointPictureForCoastTarg(wchar_t *wchOutFold, double *arrMtrxCorrFluct
		 , double *arrMtrxCorrSyst,TURPolygon Polygon, const double VAlTargCourse0, const double VAlPsi);
#endif
