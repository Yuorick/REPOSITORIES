//---------------------------------------------------------------------------

#ifndef Part_Body_3DH
#define Part_Body_3DH
#include "Plane.h"
#include "SimpleBody_3D.h"
class TPart_Body_3D
{
public:

 TPlane mPlane; // плоскость
 TSimpleBody_3D mSimpleBody_3D; //полигон в плоскости
  TPart_Body_3D();

};

#endif
