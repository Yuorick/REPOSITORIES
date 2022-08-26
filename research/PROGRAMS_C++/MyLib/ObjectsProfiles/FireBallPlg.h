//---------------------------------------------------------------------------

#ifndef FireBallPlgH
#define FireBallPlgH
class TURPointXY;
int NumPartsFB =  1 ;
int NumPointsFB =  9 ;
double PointsFB[] = {
-0.008670 , 0.516276
, 316.433494 , 16.543041
, 290.253522 , -8.434130
, 572.147330 , 0.405071
, 572.136906 , -0.365480
, 251.756775 , -18.729782
, 284.762077 , 7.631175
, 0.005852 , -0.449408
, -0.008670 , 0.516276
};


int PartsFB [] = {0};

void drawFireBall( wchar_t *pwcharrFileOut,const double  valAngRotation, const TURPointXY pntSdvig,const double valRastigenie);
#endif
