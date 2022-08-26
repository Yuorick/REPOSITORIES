//---------------------------------------------------------------------------

#ifndef ShipTargH
#define ShipTargH
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#include "PartData.h"
#include "ShipTraj.h"
#include "InitTargData.h"
#include "TargData.h"

class TPartData ;
class TTraject;
class TInitTargData ;
class TTargData;

//extern const int QUANT_TYPES_OF_TRAJECT = 3;
class TShipTarg
{
public:

// данные по количеству участков траектории  mquantParts<= 10
   int mquantParts ;
// массив параметров, описывающих участки траектории(не более 10 шт)
   TPartData marrPartData[10];
   // Текущее время
   long double mT;
   // траектория
   TShipTraj mTraject ;

   // данные по уязвимости корабля
   TTargData mTargData;

	 // для отчета
	 // к-аво точек траектории
	int mQuantPntReport ;

	 //параметр характеризующий заразарвироанную память под буфер памяти
	int mLenMemoryAlloc ;

	// буфер памяти
	 double *mparrBuff    ;

	// путь к папке с отчетом
	wchar_t *mpwcharrFoldReport ;


     // деструктор
	 ~TShipTarg() ;
	// конструктор по умолчанию
	TShipTarg () ;
	// конструктор копирования
	TShipTarg  (const TShipTarg  &R) ;
	// оператор присваивания
	TShipTarg  operator=(TShipTarg   R2) ;

	// парам конструктор
	TShipTarg ( const long double Bearing, const long double TargCourse
	, const long double TargZenitAng,  const long double V, const long double H ,
	const long double R
	,const long double valT
	,const int quantParts, TPartData *arrPartData, wchar_t* pwcharrFoldReport);;
		// парам конструктор
	TShipTarg (const TInitTargData InitData
	,const int quantParts, TPartData *arrPartData, wchar_t* pwcharrFoldReport);

	 // парам конструктор  3
 TShipTarg (const TInitTargData InitData,const int quantParts, TPartData *arrPartData,TTargData TargData, wchar_t* pwcharrFoldReport);

	bool recalcTrajPoint(const long double tNext);
	int  getNumCurrentPart(const long double valt) ;
	long double  getTimePartStarted(const long double valt) ;
	long double  getTimePartWillFinish(const long double valt) ;
	int getNumMаneuvreType(const long double valt);
	 void WriteReport() ;

	 void updateReportData() ;



}  ;
#endif
