//---------------------------------------------------------------------------

#ifndef EtalonSignH
#define EtalonSignH
class TEtalonSign
{
public:
 // амлитуда
 double mEtalonAmp;
 //дальность
 double mEtalonDist;
 //ЭПР
 double mEtalonEPR;
 //СКО внутр шума суммарной диаграммы 5П10
 double mNoiseSKZ_5P10;
 // СКЗ разброса коэффиц усиления суммарной диаграммы 5П10
 double mEtalonSigAmplFact_5P10;
 // мощность на передачу
 double mEtalonPowerPrd;
 // КУ на передачу
 double mEtalonKYPrd;
 // мощность на прием
 double mEtalonKYPriem;


 __fastcall  TEtalonSign() ;
// Конструктор копирования
__fastcall  TEtalonSign (const TEtalonSign &R2) ;
 // парам констр
  __fastcall TEtalonSign(const double EtalonAmp,const double EtalonDist, const double EtalonAPR
   , const double NoiseSKZ_5P10, const double EtalonSigAmplFact_5P10) ;

  // парам констр 2
 __fastcall TEtalonSign(const double EtalonAmp,const double EtalonDist, const double EtalonAPR
   , const double NoiseSKZ_5P10, const double EtalonSigAmplFact_5P10, const double EtalonPowerPrd,  const double  mEtalonKYPrd
   , const double);


 // оператор присваивания
 TEtalonSign   &operator=(const TEtalonSign  &R2) ;

 double calcNWaveEtalon();

};
#endif
