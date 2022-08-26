//---------------------------------------------------------------------------

#ifndef SimpleBody_3DH
#define SimpleBody_3DH
#include "Plane.h"

class TPlane;
class TSimpleBody_3D
{
public:
	 TPlane mPlane; // плоскость
	 double mM; // масса

	__fastcall TSimpleBody_3D() ;

	__fastcall TSimpleBody_3D(const TPlane Plane, const double M);

	TSimpleBody_3D  (const TSimpleBody_3D &R);

	 TSimpleBody_3D operator=(TSimpleBody_3D  R);

   //	virtual __fastcall ~TSimpleBody_3D() ;

	virtual double calcCapacity() ;

	virtual void calcCentreOfGravity(double *arrCentreGrav)  ;

	virtual  void calcInertiaMtrx(double *arrInertTens) ;

	void  calcInertiaMtrxMass(double *arrInertTens);

	static void calcInertiaMtrxShteinerSdvig(const double VAlM, double *arrS
 , double *arrInertMtrxInp, double *arrInertMtrxOut);

   void calcCentreOfGravityComplicatedAxes(double *arrCentreGrav);



};
#endif
