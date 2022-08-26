//---------------------------------------------------------------------------

#ifndef ArtComplexH
#define ArtComplexH
#include "ArtCannon.h"
class TArtCannon;
class TArtComplex
{
public:
// АУ
	TArtCannon mArtCannon;
//  точность отработки  углов
	double mSigU;
// точность отработки  скорости углов
	double mSig_dU_po_dt;

	// конструктор по умолчанию
	TArtComplex () ;
	// конструктор копирования
	TArtComplex  (const TArtComplex  &R) ;
	// оператор присваивания
	TArtComplex  &operator=(const TArtComplex  &R2) ;

	// парам конструктор1
	TArtComplex (const TArtCannon ArtCannon, const double SigU, const double Sig_dU_po_dt)  ;

}  ;
#endif
