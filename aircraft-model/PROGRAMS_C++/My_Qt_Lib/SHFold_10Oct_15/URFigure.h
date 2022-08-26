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

      TURFigure() ;

   virtual   ~TURFigure() ;

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



    static void   SwapInt(int &Data) ;
    virtual bool   LoadFromStream(TStream *Stream) ;
    virtual bool   SaveToStream(TStream *Stream) ;
    virtual void   Draw(TCanvas *Canvas, bool bSelect = false) ;
    virtual TURFigure *  PtInFigure(const TPointXY &pt) ;
    virtual TURFigure *  PtNearVisibleFigure(const TPointXY &pt) ;
    virtual void   SetOrgExt(const TPointXY &Offset, double kExt) ;

    virtual double   GetArea() ;
    virtual double   GetLength() ;
    virtual int   GetSize() ;

	__property double Area = {read = GetArea};
	__property double Length = {read = GetLength};
	__property int Size = {read = GetSize};*/
};

#endif
