//---------------------------------------------------------------------------

#ifndef SeaPlgH
#define SeaPlgH
void drawSeaSin( wchar_t *pwchShapeFileOut, const double valAmpl
   ,const double valOmega, const double valPh0, const double valXMin
   , const double valXMax, const double valDeep, const int numPoints );
void drawSeaPlane( wchar_t *pwchShapeFileOut, const double valXMin
   , const double valXMax, const double valDeep );
#endif
