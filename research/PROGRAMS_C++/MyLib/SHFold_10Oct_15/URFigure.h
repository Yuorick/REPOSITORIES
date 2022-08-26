//---------------------------------------------------------------------------

#ifndef URFigureH
#define URFigureH


//include <Graphics.hpp>

//#include "SHPStruct.h"
//#include "Layer.h"

//class TLayer ;
//---------------------------------------------------------------------------
class TURPointXY;
class TURPolygon;
class TURFigure
{
public:

	__fastcall TURFigure() ;

  // 	virtual __fastcall ~TURFigure() ;

   //	virtual  int  calcQuantsNetPoints(const double  VAlCellSize);

	void find_Objects_Type_And_Quant(wchar_t *wchFileName, int *ipShapeType, int*ipQuant);

	void ChangeByteOrder(int * pi0);

	virtual void createTargPointsArray( const int valTargCellSize, TURPointXY **ppTargPntArray, int *lenTargPntArray)  ;

	virtual void calcBoundBox();

	virtual void WriteSetSHPFiles(wchar_t *wchFileName);

	void createAimPointsArray(const int valTargCellSize, TURPointXY **ppAimPntArray, int *lenAimPntArray) ;

	virtual void createUnatedPointsArray(TURPointXY **ppunatedPoints, int *quantUnitedPoints) ;

	virtual void   LinearTransformation(const double  valAng , const TURPointXY pntSdvig,const double valRastigenie );

	int  calcDimension(const double VAlTolerance);

	virtual void   ConvexShell(TURPolygon *pPolgConv);

	virtual TURPolygon   Buffer(const double VAlBufferX,const double VAlBufferY);

	virtual TURPolygon  createBuffPolygonInAccordanceWithCorrMtrx(double *arrMtrxCorrSyst, const double VAlCoeff);


	/*
	int RecNumber ;		// Big, Номер записи
	int RecLength ;		// Big, длина записи
	ShapeType Type ; 	// тип объекта
	TLayer *Layer ;     // Слой, которому принадлежит фигура
	TColor FigureColor ;// Цвет фигуры



	static void __fastcall SwapInt(int &Data) ;
	virtual bool __fastcall LoadFromStream(TStream *Stream) ;
	virtual bool __fastcall SaveToStream(TStream *Stream) ;
	virtual void __fastcall Draw(TCanvas *Canvas, bool bSelect = false) ;
	virtual TURFigure *__fastcall PtInFigure(const TPointXY &pt) ;
	virtual TURFigure *__fastcall PtNearVisibleFigure(const TPointXY &pt) ;
	virtual void __fastcall SetOrgExt(const TPointXY &Offset, double kExt) ;

	virtual double __fastcall GetArea() ;
	virtual double __fastcall GetLength() ;
	virtual int __fastcall GetSize() ;

	__property double Area = {read = GetArea};
	__property double Length = {read = GetLength};
	__property int Size = {read = GetSize};*/
};

#endif
