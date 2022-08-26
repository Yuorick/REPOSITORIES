//---------------------------------------------------------------------------

#ifndef SeaTargPicturteH
#define SeaTargPicturteH
class TTarget;
class TURPolygon;
class TURPointXY;
void createNibourAppointmPointPictureForSeaTarg(wchar_t *wchOutFold
  ,double *arrCoMtrx_1_and_3_Groups, double *arrCoMtrx_2_Group,double *arrVectMiss_GSK, const TTarget Target, const int QuantShells);

TURPolygon createEllipsTochkiPadenia(double *arrCorrMatrxVectMiss_GSK, const TURPointXY PNtSdvig, const double VAlLevel);
#endif
