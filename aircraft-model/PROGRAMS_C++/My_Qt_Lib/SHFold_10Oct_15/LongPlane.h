#ifndef LONGPLANE_H
#define LONGPLANE_H


class TLongPlane
{
public:

 long double  marrS0[3]; // точка на плоскости - начало координат прямоугольной системы координат на плоскости  (СКП)
 long double  marrF[9]; //ортогональная  матрица перехода из сисиемы координат плоскости в исходную систему координат (скоростную)
                                    // состоит из столбцов координат ортов системы координат плоскости. Первый столбец -
                                    // это ось X СКП, последнгий столбец - вектор единичной нормали к плоскости
                                    // Преобразование координат точки из СКП в исходную СК выглядит так:
                                    // S = arrF * Sплоск +  marrS0
     TLongPlane();
 // конструктор копирования
 TLongPlane (const TLongPlane &R) ;
 // оператор присваивания
 TLongPlane &operator=(const TLongPlane & R);
 // параметр констр
 TLongPlane(long  double* arrS0,long  double* arrF);

bool findIntersectingPoint_with_Line(long double *arrPosWorking, long double *arrVeloWorking
                                     , long double *pvalPointIntersectX, long double *pvalPointIntersectY) ;

void fillNormalVect(long double *arrN);

void transform_xyzSSK_to_xyzSKP(long double *arrPointINtersect_PrSK,long  double *arrPointINtersect_SKP);

void transform_xyzSKP_to_xyzSSK(long double *arrPoint_SKP,long  double *arrPoint_PrSK);


};

#endif // LONGPLANE_H
