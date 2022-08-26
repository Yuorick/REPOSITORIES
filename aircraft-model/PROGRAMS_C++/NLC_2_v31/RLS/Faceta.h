//---------------------------------------------------------------------------

#ifndef FacetaH
#define FacetaH
//#include "Far.h"

class TFaceta
{
public:

// количество излучателей в фасете
 int m_n;
 // расстояние между излучателями в модуле
 double  m_d;
 // длина волны
 double	 mLambda ;


 __fastcall  TFaceta ::TFaceta() ;
// Конструктор копирования
__fastcall  TFaceta (const TFaceta &R2) ;
 // парам констр
 __fastcall TFaceta(const int n,const double d,const double Lambda) ;
 // парам констр
 __fastcall TFaceta(const int n);

// оператор присваивания
TFaceta   operator=(TFaceta  R2) ;
double  fncFSource (const double valTetta);
double  fnc_dFSource_po_dTetta (const double valTetta);
double  fncFFaceta (const double valTetta);
double  fnc_dFFaceta_po_dTet (const double tet);



static   double fnc_dF3_po_tet(const double val_gam, const double tet) ;

static   double fnc_F3(const double val_gam, const double tet) ;

static   double fnc_F4(const double v,const  double u, const double v1,const  double u1) ;

double findDiagrWidth();

void createDiagrGraphs(wchar_t *wchFoldName1 );

 double findDiagrWidthApprox();

 double  fncFFacetaApprox (const double valTetta) ;

 double  fnc_dFFacetaApprox_po_dTet (const double tet);

 double fnc_d2FFacetaApprox_po_dTet2(const double tet);
};
#endif
