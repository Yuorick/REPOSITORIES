//---------------------------------------------------------------------------

#ifndef QPlaneH
#define QPlaneH

class TURPointXY;
class QPlane
{
public:

 double  marrS0[3]; // точка на плоскости - начало координат прямоугольной системы координат на плоскости  (СКП)
 double  marrF[9]; //ортогональная  матрица перехода из сисиемы координат плоскости в исходную систему координат (скоростную)
            // состоит из столбцов координат ортов системы координат плоскости. Первый столбец -
            // это ось X СКП, последнгий столбец - вектор единичной нормали к плоскости
            // Преобразование координат точки из СКП в исходную СК выглядит так:
            // S = arrF * Sплоск +  marrS0
     QPlane();
 // конструктор копирования
 QPlane (const QPlane &R) ;
 // оператор присваивания
 QPlane &operator=(const QPlane  &R);
 // параметр констр
 QPlane( double* arrS0, double* arrF);

bool findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking
                                     , double *pXintersect, double *pYintersect) ;

void fillNormalVect(double *arrN);

void transform_xyzSSK_to_xyzSKP(double *arrPointINtersect_PrSK, double *arrPointINtersect_SKP);

void transform_xyzSKP_to_xyzSSK(double *arrPoint_SKP, double *arrPoint_PrSK);

bool findIntersectingPoint_with_Line(double *arrPosWorking, double *arrVeloWorking
                             , double *pXintersect, double *pYintersect, double *palf
                                             , double *pt);

int createInputDataReport(wchar_t*FileName, const bool bHeader);

bool findIntersectingPoint_with_Line_PrSK(double *arrPosWorking, double *arrVeloWorking
                             , double *arrPointINtersect_PrSK, double *palf
                                             , double *pt);


};
#endif
