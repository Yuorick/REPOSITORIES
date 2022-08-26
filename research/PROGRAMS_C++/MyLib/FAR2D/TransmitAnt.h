//---------------------------------------------------------------------------

#ifndef TransmitAntH
#define TransmitAntH
// передающая антенна
class TTransmitAnt
{
public:

 // мощность на передачу
	double mPowerPrd;
	// КУ на передачу
	double mKYPrd;




 __fastcall  TTransmitAnt() ;
// Конструктор копирования
__fastcall  TTransmitAnt (const TTransmitAnt &R2) ;

 // оператор присваивания
 TTransmitAnt   &operator=(const TTransmitAnt  &R2) ;

  // парам констр
 __fastcall TTransmitAnt(const double PowerPrd,const double KYPrd);





};
#endif
