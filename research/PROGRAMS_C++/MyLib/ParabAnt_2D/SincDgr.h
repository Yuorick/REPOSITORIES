//---------------------------------------------------------------------------

#ifndef SincDgrH
#define SincDgrH
class TComp;
class TURPointXY;
class TSincDgr
{
public:

 // ????????
 double  mAppert;

 //????????? ???? ?  ????????e
 double mNoiseDisp;

 // ??????? ??????? ???????? k ? ??????????? ?????????
 // ?? ????, ??????? ???????? ????? P = (1 + k)P0
double mAmplFactSig;




 TSincDgr() ;


// ??????????? ???????????
 TSincDgr(const TSincDgr &R) ;
 TSincDgr &operator=(const TSincDgr  &R2) ;
 // ????? ??????
 TSincDgr (const double RadTet05, const double NoiseSkz, const double AmplFactSig
	 , const double Lambda);

	 // ????? ??????   2

	TSincDgr(const double  Appert);

 double transformAngToGeneralizedAng (const double  valAng, const double Lambda  );

 double transformGeneralizedAngToAng (const double  GeneralizedAng , const double Lambda );

 void ImitateMeasure( double valTargGen, TComp cmpKTarg, double  valAntpGen,TComp cmpKAntp
		, TComp *pcmpS, TComp *pcmpSZv);

double fncDiagrFromRad(const double VAlTetRad, const double Lambda ) ;

void createSHP_Graph(wchar_t *Fold, const double Lambda, int numRoot
	   , TURPointXY pntSdvig, double scalex, double scaley);

double fncDerivDgr_po_dTet(const double VAlTetRad, const double Lambda );

double fncDelta05(const double Lambda );

double fnc_qu( const double Lambda ) ;

double fncDerivDgr_po_dDelta05(const double VAlTetRad, const double Lambda ) ;

double fncDerivDgr_po_dLambda(const double VAlTetRad, const double Lambda )   ;

void createSHP_Graph_SKZ_From_Tet(wchar_t *wchFileName, const double Lambda,const double VAlSig0);


};
#endif
