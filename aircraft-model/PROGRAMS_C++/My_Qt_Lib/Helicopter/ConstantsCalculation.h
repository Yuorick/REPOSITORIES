//---------------------------------------------------------------------------

#ifndef ConstantsCalculationH
#define ConstantsCalculationH
double calcCt(const double TKel0,const double HeliH, const double BladeOmega
	, const double BladeR,const	double HelicMass) ;

double calc_rQ(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAPofile_d1 ) ;

double calc_X_BladeCentreMass(const double VAl_L, const double VAlPofile_d0,const double VAlPofile_d1 ) ;

double calc_X_StatMoment_Sg(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAlPofile_d1 , const double VAlM );

double calcInertiaMoment0(const double VAlL, const double VAlPofile_d0,const double VAlPofile_d1,
  const double VAlM );

double calcInertiaMoment(const double VAlR, const double VAlRadHorizHsarnir
	, const double VAlPofile_d0,const double VAlPofile_d1 , const double VAlM  );

double calcInertiaMomentHorSharnir(const double VAlL, const double VAlPofile_d0,const double VAlPofile_d1,
  const double VAlM );

double calc_C_y_alfa(const double TKel0,const double HeliH, const double BladeOmega
	, const double BladeR , const double VAlRadHorizHsarnir, const double VAlBlade_b
	, const int QUantBlades, const	double HelicMass);
#endif
