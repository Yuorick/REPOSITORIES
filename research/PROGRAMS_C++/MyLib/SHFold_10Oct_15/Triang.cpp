


#pragma hdrstop

#include "Triang.h"
#include <stdio.h>
 #include <math.h>
//#include "UrPointXY.h"
//#include "SHPStruct.h"
#include "URPolygon.h"
#include "YrLinTransform.h"
#include "UrPointXY.h"
#include <string.h>
#include <stdlib.h>

extern const double TOLRNC ;

  TTriang ::TTriang()
{
  memset(m_pVert,0,4 * sizeof(double));
}

 // конструктор копирования
 TTriang ::TTriang (const TTriang &R)
 {
	for (int i = 0 ; i < 4; i++) m_pVert[i] = R.m_pVert[i];

 }
 // оператор присваивания
 TTriang &TTriang::operator=(const TTriang  &R)
 {

  for (int i = 0 ; i < 4; i++) m_pVert[i] = R.m_pVert[i];
  return *this ;
 }
 // парам конструктор

TTriang :: TTriang( TURPointXY  *arrPoints)
 {
   if(arrPoints != NULL)
   {
   for (int i = 0 ; i < 3; i++)
   {
	m_pVert[i].X = arrPoints[i].X;
	m_pVert[i].Y = arrPoints[i].Y;
   }
   m_pVert[3].X  = m_pVert[0].X ;
   m_pVert[3].Y  = m_pVert[0].Y ;
   }

 }
 // парам конструктор
/* TTriang :: TTriang( TURPointXY  p0,TURPointXY  p1,TURPointXY  p2)
 {
  m_pVert[0].X = p0.X ;
   m_pVert[0].Y = p0.Y ;
   m_pVert[1].X = p1.X ;
   m_pVert[1].Y = p1.Y ;
   m_pVert[2].X = p2.X ;
   m_pVert[2].Y = p2.Y ;
   m_pVert[3] = m_pVert[0] ;
 }  */
  // парам конструктор
 TTriang :: TTriang( TURPointXY  p0,TURPointXY  p1,TURPointXY  p2)
 {

   m_pVert[0] = p0 ;
   m_pVert[1] = p1 ;
   m_pVert[2] = p2 ;
   m_pVert[3] = p0 ;

 }
double TTriang ::calcSq()
{
	double S =0 ;
	for (int i =0; i < 3; i++)
	{
	   S = S + calcVectS(m_pVert[i],m_pVert[i + 1]);
	}
	return -S;
}
// площадь треугольника с центром в (0,0)
double TTriang ::calcVectS(TURPointXY P1,TURPointXY P2)
{
	return (P1.X * P2.Y - P1.Y* P2.X)/2;
}

 double TTriang ::max_d(const double x0, const double x1 )
 {
	 return (x0 > x1)? x0:x1;
 }
 double TTriang ::min_d(const double x0, const double x1 )
 {
	 return (x0 > x1)? x1:x0;
 }
// Вычитание треугольника из треугольника trT0\ trT1
//Если треугольники не пересекаются то  quantTr = -2
// Если весь 0треугольник trT0 лежит вутри trT1 и результат вычитания есть пустое множество, то  quantTr = 0


 void TTriang ::SubTwoTrsToTrs( TTriang *ptrT0, TTriang *ptrT1, TTriang *pTr,int &quantTr )
{
	// для отладки
   /*	 TURPolygon plgTr0(1, ptrT0);
	 TURPolygon plgTr1(1, ptrT1);
	 plgTr0.WriteToASCII( L"D:\\GeoProcC++\\FILES\\plgTr0.txt");
	 plgTr1.WriteToASCII( L"D:\\GeoProcC++\\FILES\\plgTr1.txt");
	*/
	///
   TURPolygon urplgArr[3];
   for (int i = 0; i < 3; i++) urplgArr[i] = TURPolygon(12);


   int  quantPlg = -1;
   quantTr = 0;
   SubTwoTrsToRings(ptrT0, ptrT1, urplgArr,quantPlg );

   // для отладки
  //	urplgArr[0].WriteToASCII(L"D:\\GeoProcC++\\FILES\\urplgArr0.txt") ;
  //	urplgArr[1].WriteToASCII(L"D:\\GeoProcC++\\FILES\\urplgArr1.txt") ;
   ///
   if( (quantPlg == -2) || (quantPlg == 0) )
   {
	quantTr = quantPlg;
	return;

   }
   for (int i = 0; i < quantPlg; i++)
   {
	 if( urplgArr[i].NumPoints == 4)
	 {
		pTr[quantTr] =  TTriang(urplgArr[i].Points);
		quantTr++;
	 }
	 else
	 {
	   urplgArr[i].TriangulateRing(&pTr[quantTr]) ;
	   quantTr += urplgArr[i].NumPoints -3;
	 }
   }
}
// Вычитание треугольника из треугольника (*trT0)\ (*ptrT1). (*ptrT0) -положительный, (*ptrT1) - отрицательный
//Если треугольники не пересекаются то  quantPlg = -2
// Если весь треугольник *trT0 лежит вутри *trT1 и результат вычитания есть пустое множество, то  quantPlg = 0
// OUTPUT: urplgArr - массив полигонов, quantPlg - к-во элементов этого массива
// перед обращением к этой функции следует создать массив urplgArr[3] с количеством вершин каждого полигона 12
// TURPolygon urplgArr[3]; for (int i = 0; i < 3; i++) urplgArr[i] = TURPolygon(12);
// результатом вычитания из одного треуг другого является максимум 3 полигона. Количество вершин у каждого полигона
// не предышает 12
void TTriang ::SubTwoTrsToRings( TTriang *ptrT0, TTriang *ptrT1, TURPolygon *urplgArr,int &quantPlg )
{

  TURPolygon plgCircle0(12),plgCircle1(12);
  memcpy( plgCircle0.Points,(*ptrT0).m_pVert,4 * sizeof(TURPointXY)) ;
  memcpy( plgCircle1.Points,(*ptrT1).m_pVert,4 * sizeof(TURPointXY)) ;
  plgCircle0.NumPoints = 4;
  plgCircle1.NumPoints = 4;
  int iarrNums0[12] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int iarrNums1[12] = { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int iarrVerNums0[4] = {0,1,2,3};
  int iarrVerNums1[4] = {0,1,2,3};
//  bool barrOut[12] = {true,true,true,false,false,false,false,false,false,false,false,false};
  int iCount =0;
  int quantCross = 0 ; //количество пересечений сторон
  int quantJoinVerts = 0; // количество общих вершин
  int quantTouchs = 0;  // количество касаний вершин и сторон
  int quantJoinLines = 0;
  int iCase = -100;

  int quantCase4 = 0 ;
  for (int i = 0; i < 3; i++)
  for (int j = 0; j < 3; j++)
  {
	 TURPointXY p0(0,0);
	 int num0,num1;
	 int ibegin0 = iarrVerNums0[i] ;
	 int ibegin1 = iarrVerNums1[j] ;
	// int numEdge0 = ibegin0;
	// int  numEdge1 = ibegin1;
	if (   IntersectTwoSegments( (*ptrT0).m_pVert[i], (*ptrT0).m_pVert[i +1]
					,(*ptrT1).m_pVert[j], (*ptrT1).m_pVert[j +1]
					,&p0 ,&iCase))
	{
	  iCount ++;
	  // Номер ребра полигона plgCircle0 на котором лежит точка p0
	  // сразу три стороны треугольника не могум пересекатьодин отрезок.
	  // Поэтому, разница между iarrVerNums0[i + 1] и iarrVerNums0[i ] не может превосходить 2
	  //  то есть на каждой  стороне треугольника может быть не более 2 дополгнительных вершин
		int numEdge0 = ibegin0;
		if (( iarrVerNums0[i + 1] - iarrVerNums0[i] ) != 1 )
		{
		  TURPointXY point0(0,0),point1(0,0);
		  point0.X = plgCircle0.Points[ibegin0 +1].X - plgCircle0.Points[ibegin0 ].X;
		  point0.Y = plgCircle0.Points[ibegin0 +1].Y - plgCircle0.Points[ibegin0 ].Y ;
		  point1.X = p0.X - plgCircle0.Points[ibegin0 +1].X ;
		  point1.Y = p0.Y - plgCircle0.Points[ibegin0 +1].Y ;
		  if ( TURPointXY::ScalMult(point1,point0)>0 ) numEdge0++;

		}

	  // Номер ребра полигона plgCircle1 на котором лежит точка p0
		int  numEdge1 = ibegin1;
		if (( iarrVerNums1[j + 1] - iarrVerNums1[j] ) != 1 )
		{
		  TURPointXY point0(0,0),point1(0,0);
		  point0.X = plgCircle1.Points[ibegin1 +1].X - plgCircle1.Points[ibegin1 ].X;
		  point0.Y = plgCircle1.Points[ibegin1 +1].Y - plgCircle1.Points[ibegin1 ].Y ;
		  point1.X = p0.X - plgCircle1.Points[ibegin1 +1].X ;
		  point1.Y = p0.Y - plgCircle1.Points[ibegin1 +1].Y ;
		if ( TURPointXY::ScalMult(point1,point0)>0 ) numEdge1++;

		}

		switch(iCase)
		{
			case 25:
				quantTouchs++;
				if (iarrNums0 [ iarrVerNums0[i + 1] % (plgCircle0.NumPoints-1)] == -1 )
				{
					plgCircle1.NumPoints ++;
					for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > numEdge1 )) iarrNums0[ii] ++;
					InsertNum(iarrNums1,numEdge1,iarrVerNums0[i + 1] % (plgCircle0.NumPoints-1));
					iarrNums0[ iarrVerNums0[i + 1] % (plgCircle0.NumPoints-1) ] = (numEdge1 + 1) % (plgCircle1.NumPoints-1) ;
					for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > numEdge1) iarrVerNums1[k] ++;
					InsertVert(&plgCircle1,numEdge1,(*ptrT0).m_pVert[i +1]);
					break;
				}
				else
				break;
			case 26:
				quantTouchs++;
				if (iarrNums0 [ iarrVerNums0[ i ]]  == -1 )
				{
					plgCircle1.NumPoints ++;
					for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > numEdge1 )) iarrNums0[ii] ++;
					InsertNum(iarrNums1,numEdge1,iarrVerNums0[i] );
					iarrNums0[ ibegin0  ] = (numEdge1 + 1) % (plgCircle1.NumPoints-1) ;
					iarrNums1[ (numEdge1 +1) % (plgCircle1.NumPoints-1) ] =  ibegin0  ;
					for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > numEdge1) iarrVerNums1[k] ++;
					InsertVert(&plgCircle1,numEdge1,(*ptrT0).m_pVert[i]);
					break;
				}
				else
				break;

			case 27:
				quantTouchs++;
				if (iarrNums1 [iarrVerNums1 [j + 1] % (plgCircle1.NumPoints-1)] == -1 )
				{
					plgCircle0.NumPoints ++;
					for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > numEdge0 )) iarrNums1[ii] ++;
					InsertNum(iarrNums0,numEdge0,iarrVerNums1 [j + 1] % (plgCircle1.NumPoints-1));
					InsertVert(&plgCircle0,numEdge0,(*ptrT1).m_pVert[j +1]);
					for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > numEdge0) iarrVerNums0[k] ++;
				   //	iarrNums1[ iarrVerNums1 [j + 1] % (plgCircle1.NumPoints-1) ] =  (ibegin0 + 1) % (plgCircle0.NumPoints-1) ;
					iarrNums1[ iarrVerNums1 [j + 1] % (plgCircle1.NumPoints-1) ] =  numEdge0 + 1;
					//InsertVert(&plgCircle0,ibegin0,(*ptrT1).m_pVert[j +1]);
					break;
				}
				else
				break;

			case 28:
				quantTouchs++;
				if (iarrNums1 [iarrVerNums1[j] ]  == -1 )
				{
					plgCircle0.NumPoints ++;
					for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > numEdge0 )) iarrNums1[ii] ++;
					InsertNum(iarrNums0,numEdge0,ibegin1 );
					iarrNums1[ ibegin1  ] =  (numEdge0 + 1) % (plgCircle0.NumPoints-1) ;
					for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > numEdge0) iarrVerNums0[k] ++;
					InsertVert(&plgCircle0,numEdge0,(*ptrT1).m_pVert[j]);
					break;
				}
				else
				break;

			case 24:
				quantCross++;
			   /*	if ((iarrVerNums0[i +1] - iarrVerNums0[i]) > 1)
				{
				if( (fabs((*p0).X - plgCircle0.Points[ibegin0 ].X ) + fabs((*p0).Y - plgCircle0.Points[ibegin0 ].Y ) )
				>  (fabs(plgCircle0.Points[ibegin0 +1 ].X - plgCircle0.Points[ibegin0 ].X )
				+ fabs(plgCircle0.Points[ibegin0 +1].Y  - plgCircle0.Points[ibegin0 ].Y ) ))
				ibegin0++;
				}


				if ((iarrVerNums1[j +1] - iarrVerNums1[j]) > 1)
				{
				if( (fabs((*p0).X - plgCircle1.Points[ibegin1 ].X ) + fabs((*p0).Y - plgCircle1.Points[ibegin1 ].Y ) )
				>  (fabs(plgCircle1.Points[ibegin1 +1 ].X - plgCircle1.Points[ibegin1 ].X ) + fabs(plgCircle1.Points[ibegin1 +1].Y  - plgCircle1.Points[ibegin1 ].Y ) ))
				ibegin1++;
				}
				plgCircle0.NumPoints++;
				plgCircle1.NumPoints++;
				InsertVert(&plgCircle0,ibegin0,*p0);
				InsertVert(&plgCircle1,ibegin1,*p0);
				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums0[ii]++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii]++;

				InsertNum(iarrNums0,ibegin0,(ibegin1 +1) % (plgCircle1.NumPoints-1));
				InsertNum(iarrNums1,ibegin1,(ibegin0 +1) % (plgCircle0.NumPoints-1));

				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k]++ ;
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k]++ ;  */
				plgCircle0.NumPoints++;
				plgCircle1.NumPoints++;
				InsertVert(&plgCircle0,numEdge0,p0);
				InsertVert(&plgCircle1,numEdge1,p0);
				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > numEdge1 )) iarrNums0[ii]++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > numEdge0 )) iarrNums1[ii]++;

				InsertNum(iarrNums0,numEdge0,(numEdge1 +1) % (plgCircle1.NumPoints-1));
				InsertNum(iarrNums1,numEdge1,(numEdge0 +1) % (plgCircle0.NumPoints-1));
				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > numEdge0) iarrVerNums0[k]++ ;
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > numEdge1) iarrVerNums1[k]++ ;



			break;
			case 29:
			quantJoinVerts++;
			iarrNums0[iarrVerNums0[i]] = iarrVerNums1[j] ;
			iarrNums1[iarrVerNums1[j]] = iarrVerNums0[i] ;
			break;
			case 30:
			quantJoinVerts++;
			iarrNums0[iarrVerNums0[i]] = iarrVerNums1[j+1] % (plgCircle1.NumPoints-1) ;
			iarrNums1[iarrVerNums1[j +1]%(plgCircle1.NumPoints-1)] = iarrVerNums0[i] ;

			break;
			case 31:
			quantJoinVerts++;
			iarrNums0[iarrVerNums0[i+1]% (plgCircle0.NumPoints-1)] = iarrVerNums1[j] ;
			iarrNums1[iarrVerNums1[j]] = iarrVerNums0[i + 1] % (plgCircle1.NumPoints-1) ;

			break;
			case 32:
			quantJoinVerts++;
			iarrNums0[iarrVerNums0[(i + 1)] % (plgCircle0.NumPoints -1)] = iarrVerNums1[(j+ 1)] %(plgCircle1.NumPoints -1) ;
			iarrNums1[iarrVerNums1[(j + 1)] % (plgCircle1.NumPoints -1)] = iarrVerNums0[(i +1)] %(plgCircle0.NumPoints -1);

			break;
			case 0:
			quantJoinLines++;
			plgCircle0.NumPoints ++;
			plgCircle1.NumPoints ++;

			for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums0[ii] ++;
			for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii] ++;

			InsertNum(iarrNums0,ibegin0,(ibegin1 +2) % (plgCircle1.NumPoints-1));
			InsertNum(iarrNums1,ibegin1,(ibegin0 +2) % (plgCircle0.NumPoints-1));
			 iarrNums0[ (ibegin0 +2) % (plgCircle0.NumPoints-1) ] = (ibegin1 + 1) % (plgCircle1.NumPoints-1) ;
			 iarrNums1[ (ibegin1 +2) % (plgCircle1.NumPoints-1) ] =  (ibegin0 + 1) % (plgCircle0.NumPoints-1) ;
			for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;
			for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;

			InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[ (ibegin1 + 1) % (plgCircle1.NumPoints-1)]);
			InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[ (ibegin0 + 2) % (plgCircle0.NumPoints-1)]);

			break;
			case 1:
			quantJoinLines++;
			plgCircle0.NumPoints ++;
			plgCircle1.NumPoints ++;

			for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums0[ii] ++;
			for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii] ++;

			InsertNum(iarrNums0,ibegin0,ibegin1);
			InsertNum(iarrNums1,ibegin1,ibegin0);

			for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;
			for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;

			InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[ibegin1 ]);
			InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[ibegin0 ]);

			break;
			case 2:
			quantJoinLines++;
			if (iarrNums1[ ibegin1] == -1)
			{
			  if (iarrNums1[( ibegin1 +1)% (plgCircle1.NumPoints-1)] == -1)
			  {
				plgCircle0.NumPoints += 2;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii] += 2;
				InsertNum(iarrNums0,ibegin0,ibegin1);
				InsertNum(iarrNums0,ibegin0,(ibegin1 + 1) % (plgCircle1.NumPoints-1));
				iarrNums1[ ibegin1] = (ibegin0 + 2) % (plgCircle0.NumPoints-1);
				iarrNums1[ (ibegin1 + 1) % (plgCircle1.NumPoints-1)] = (ibegin0 + 1) % (plgCircle0.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] += 2;

				InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[ibegin1 ]);
				InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[ (ibegin1 + 1) % (plgCircle1.NumPoints-1) ]);

				break;

			  }
			  else
			  {
				plgCircle0.NumPoints ++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > (ibegin0 +1) )) iarrNums1[ii] ++;
				InsertNum(iarrNums0,ibegin0 +1,ibegin1);
				iarrNums1[ ibegin1] = (ibegin0 + 2) % (plgCircle0.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;

				InsertVert(&plgCircle0,ibegin0 +1,plgCircle1.Points[ibegin1 ]);

				break;

			  }
			}
			else
			{
			  if (iarrNums1[( ibegin1 +1)% (plgCircle1.NumPoints-1)] == -1)
			  {
				plgCircle0.NumPoints ++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0  )) iarrNums1[ii] ++;
				InsertNum(iarrNums0,ibegin0 ,(ibegin1 +1) % (plgCircle1.NumPoints-1));
				iarrNums1[ (ibegin1 +1) % (plgCircle1.NumPoints-1)] = (ibegin0 + 1) % (plgCircle0.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;

				InsertVert(&plgCircle0,ibegin0 ,plgCircle1.Points[ibegin1 +1 ]);

				break;

			  }
			  else
			  {
				break;
			  }

			}
			case 3:
			quantJoinLines++;
			if (iarrNums0[ ibegin0] == -1)
			{
			  if (iarrNums0[( ibegin0 +1)% (plgCircle0.NumPoints-1)] == -1)
			  {
				plgCircle1.NumPoints += 2;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums1[ii] += 2;
				InsertNum(iarrNums1,ibegin1,ibegin0);
				InsertNum(iarrNums1,ibegin1,(ibegin0 + 1) % (plgCircle0.NumPoints-1));
				iarrNums0[ ibegin0] = (ibegin1 + 2) % (plgCircle1.NumPoints-1);
				iarrNums1[ (ibegin0 + 1) % (plgCircle0.NumPoints-1)] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] += 2;

				InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[ibegin0 ]);
				InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[ (ibegin0 + 1) % (plgCircle0.NumPoints-1) ]);

				break;

			  }
			  else
			  {
				plgCircle1.NumPoints ++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > (ibegin1 +1) )) iarrNums0[ii] ++;
				InsertNum(iarrNums1,ibegin1 +1,ibegin0);
				iarrNums0[ ibegin0] = (ibegin1 + 2) % (plgCircle1.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;

				InsertVert(&plgCircle1,ibegin1 +1,plgCircle0.Points[ibegin0 ]);

				break;

			  }
			}
			else
			{
			  if (iarrNums0[( ibegin0 +1)% (plgCircle0.NumPoints-1)] == -1)
			  {
				plgCircle1.NumPoints ++;
				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1  )) iarrNums0[ii] ++;
				InsertNum(iarrNums1,ibegin1 ,(ibegin0 +1) % (plgCircle0.NumPoints-1));
				iarrNums0[ (ibegin0 +1) % (plgCircle0.NumPoints-1)] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;

				InsertVert(&plgCircle1,ibegin1 ,plgCircle0.Points[ibegin0 +1 ]);

				break;

			  }
			  else
			  {
				break;
			  }

			}


			case 4:
			quantJoinLines++;
			quantCase4 ++;
			if (quantCase4 == 2)
			{
			quantPlg=0;
			return;
			}
			 iarrNums0[ ibegin0] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
			 iarrNums0[ (ibegin0 + 1) % (plgCircle0.NumPoints-1)] = ibegin1 ;
			 iarrNums1[ ibegin1] = (ibegin0 + 1) % (plgCircle0.NumPoints-1);
			 iarrNums1[ (ibegin1 + 1) % (plgCircle1.NumPoints-1)] = ibegin0;

			break;

			case 5:
			quantJoinLines++;
			if (iarrNums1[ (ibegin1 + 1)% (plgCircle1.NumPoints-1)] == -1)
			{
			plgCircle0.NumPoints ++;
			for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii] ++;

			InsertNum(iarrNums0,ibegin0,(ibegin1 +1) % (plgCircle1.NumPoints-1));
			iarrNums1[(ibegin1 + 1) % (plgCircle1.NumPoints-1)] = (ibegin0 +1) % (plgCircle0.NumPoints-1);
			for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;
			InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[(ibegin1 +1) % (plgCircle1.NumPoints-1)]);
			}
			if (iarrNums1[ ibegin1 ] == -1)
			{
			  iarrNums0[(ibegin0 + 2) % (plgCircle0.NumPoints-1)] = ibegin1;
			  iarrNums1[ibegin1 ] = (ibegin0 + 2) % (plgCircle0.NumPoints-1);
			}
			break;
			case 6:
			quantJoinLines++;
			  if (iarrNums1[ ibegin1 ] == -1)
			  {
				plgCircle0.NumPoints ++;

				for (int ii = 0; ii < 12; ii++) if ((iarrNums1[ii] > ibegin0 )) iarrNums1[ii] ++;

				InsertNum(iarrNums0,ibegin0,ibegin1 );
				iarrNums0[ibegin0 ] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				iarrNums1[ibegin1 ] = (ibegin0 + 1) % (plgCircle0.NumPoints-1);
				iarrNums1[(ibegin1 + 1) % (plgCircle1.NumPoints-1)] = ibegin0 ;
				for (int k= 0; k < 4; k++) if (iarrVerNums0[k] > ibegin0) iarrVerNums0[k] ++;
				InsertVert(&plgCircle0,ibegin0,plgCircle1.Points[ibegin1 ]);
				break;
			  }
			  else
			  {
				iarrNums0[ibegin0 ] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				iarrNums1[(ibegin1 + 1) % (plgCircle1.NumPoints-1)] = ibegin0 ;
				break;
			  }

			case 7:
			quantJoinLines++;
			  if (iarrNums1[ ibegin1 ] == -1)
			  {
				plgCircle1.NumPoints ++;

				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums0[ii] ++;

				InsertNum(iarrNums1,ibegin1,ibegin0 );
				iarrNums1[ibegin1 ] = (ibegin0 + 1) %  (plgCircle0.NumPoints-1);
				iarrNums0[ibegin0 ] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				iarrNums0[(ibegin0 + 1)% (plgCircle0.NumPoints-1)] = ibegin1 ;
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;
				InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[ibegin0 ]);
				break;
			  }
			  else
			  {
				iarrNums1[ibegin1 ] = (ibegin0 + 1) % (plgCircle0.NumPoints-1);
				iarrNums0[(ibegin0 + 1) % (plgCircle0.NumPoints-1)] = ibegin1 ;
				break;
			  }

		   //	break;
			case 8:
			quantJoinLines++;
			  if (iarrNums0[ (ibegin0 + 1) % (plgCircle0.NumPoints-1) ] == -1)
			  {
				plgCircle1.NumPoints ++;

				for (int ii = 0; ii < 12; ii++) if ((iarrNums0[ii] > ibegin1 )) iarrNums0[ii] ++;

				InsertNum(iarrNums1,ibegin1,(ibegin0 +1) % (plgCircle0.NumPoints-1) );
				iarrNums1[(ibegin1 +2) % (plgCircle1.NumPoints-1) ] = ibegin0 ;
				iarrNums0[ibegin0 ] = (ibegin1 + 2) % (plgCircle1.NumPoints-1);
				iarrNums0[(ibegin0 + 1)% (plgCircle0.NumPoints-1)] = (ibegin1 +1)% (plgCircle1.NumPoints-1) ;
				for (int k= 0; k < 4; k++) if (iarrVerNums1[k] > ibegin1) iarrVerNums1[k] ++;
				InsertVert(&plgCircle1,ibegin1,plgCircle0.Points[(ibegin0 +1)% (plgCircle1.NumPoints-1) ]);
				break;

			  }
			  else
			  {
				iarrNums0[ibegin0 ] = (ibegin1 + 1) % (plgCircle1.NumPoints-1);
				iarrNums1[(ibegin1 + 1) % (plgCircle1.NumPoints-1)] = ibegin0 ;
				break;
			  }
			case 9:
			quantPlg = -2;
			urplgArr[0] = TURPolygon (4,(*ptrT0).m_pVert);
			return;

		   //	break;
			case 10:
			quantPlg = -2;
			urplgArr[0] = TURPolygon (4,(*ptrT0).m_pVert);
             return;
		  //	break;

			default:  // треугольники не пересекаются   (iCase >= 9) && (iCase <24)
			quantPlg = -2;
			urplgArr[0] = TURPolygon (4,(*ptrT0).m_pVert);
			return;
		  //	break;


	   }

	}
 }

	//if (iCount == 0)
   if( (quantCross == 0) && (quantJoinLines == 0) && ( (quantJoinVerts == 0)|| (quantTouchs == 0))  )

	{
	  TURPointXY pp =  (*ptrT1).m_pVert[0];
	  if ( (*ptrT0).PntInTriangle(pp) == 1)
	  {  //  (*ptrT1) внутри (*ptrT0)
		SubTwoEnclosedTrToRings(  &(*ptrT0), &(*ptrT1), urplgArr );
		quantPlg = 2 ;
		return ;
	  }

	  else
	  {  //  (*ptrT0) внутри (*ptrT1)
		  TURPointXY pp =  (*ptrT0).m_pVert[0];
		  if ( (*ptrT1).PntInTriangle(pp) == 1)
		  {

			quantPlg = 0 ;
			return ;
		  }
		  else
		  {  //не пересекаются
			quantPlg = -2 ;
			return ;

          }
		//quantPlg = 0;
		//urplgArr[0] = TURPolygon (4,(*ptrT0).m_pVert);
		//return;
	  }
   }
  // пересечение треугольников равно единственной вершине(общей)
   if( (quantCross == 0) && (quantJoinLines == 0) && ( (quantJoinVerts == 4)|| (quantTouchs == 2))  )
   {
	 // проверяем существует ли вершина треуг (*ptrT0) лежащая внутри (*ptrT1)
	 // если существует, то весь (*ptrT0) лежит внутри (*ptrT1)
	 for (int i =0; i < 3; i++)
	 {
		if((*ptrT1).PntInTriangle((*ptrT0).m_pVert[i]) == 1)
		{
		  quantPlg = 0;
		  return ;
		}
	 }

   //
	 // проверяем существует ли вершина треуг (*ptrT1) лежащая вне (*ptrT1)
	 // если существует, то весь (*ptrT1) лежит вне (*ptrT0)
	 for (int i =0; i < 3; i++)
	 {
		if((*ptrT0).PntInTriangle((*ptrT1).m_pVert[i]) == 0)
		{
		quantPlg = -2;
		urplgArr[0] = TURPolygon (4,(*ptrT0).m_pVert);
		return;
		}
	 }
   }
   //

   // пересечение треугольников равно единственной точке являющейся
   // вершиной одного из треугольников и токе на стороне другого треугольника
 //  if( (quantCross == 0) && (quantJoinLines == 0) && (quantTouchs == 2) )
  // {

 //  }
	 // для отдладки   !!!!
	  //plgCircle0.

	 ///
	  bool barrUsed[3] = {false,false,false};
	   quantPlg = 0;
	  for (int i = 0; i < 3; i++)
	  {

	   TURPointXY pntarrTemp[24] ;
	  int lenarr = 0;
		bool bProceed = false;
 // построение цикла из вершины i треугольника (*ptrT0) не лежащей в треугольнике (*ptrT1)
		if( ( (*ptrT1).PntInTriangle((*ptrT0).m_pVert[i]) == 0) && (!barrUsed[i ]))
		{
			int index0  =iarrVerNums0[i];
			int index = index0;
			bool bCircle0 = true;
			barrUsed[i ] = true ;
			for (int j = 0; j< 24; j++)
			{
			if( bCircle0 )
			{

				pntarrTemp[j] = plgCircle0.Points[index] ;
				lenarr++;
					if ( ( (index%plgCircle0.NumPoints) == index0) && ( j > 0))
					{
					  urplgArr[quantPlg].NumPoints = lenarr;
					  memcpy(urplgArr[quantPlg].Points, pntarrTemp,lenarr * sizeof(TURPointXY));
					  quantPlg++ ;
					  break ;
					}
				for (int k =0; k < 3; k++)
				{
					if (iarrVerNums0 [k] == index)
					{
					barrUsed[k ] = true;
					break;
					}
				}
				if (iarrNums0[index] == -1) //|| ( bProceed == true))
				{
				index = (index +1)% (plgCircle0.NumPoints -1);
				}
				else
				{
				  if ( bProceed )
				  {
					index = (index +1)% (plgCircle0.NumPoints -1);
					bProceed = false ;

				  }
				  else
				  {
				   index = iarrNums0[index ];
					bProceed = true ;
				   bCircle0 = false ;
				   }
				}


			}
			else
			{
				pntarrTemp[j] = plgCircle1.Points[index] ;
				lenarr++;
				if (iarrNums1[index] == -1) //|| ( bProceed == true))
				{
					index = (index +1)% (plgCircle1.NumPoints -1);
				}
				else
				{
				  if (bProceed)
				  {
					index = (index +1)% ( plgCircle1.NumPoints -1);
					bProceed = false ;
				  }
				   else
				   {
				   index = iarrNums1[index ];

					bProceed = true;;
					bCircle0 = true ;
					if ( (index % (plgCircle0.NumPoints -1)) == index0)
					{
					  pntarrTemp[lenarr] = plgCircle0.Points[index0] ;
					  lenarr++;
					  urplgArr[quantPlg] = TURPolygon(lenarr, pntarrTemp) ;
					  quantPlg++ ;
					  break ;
					}
				   }
				}
			}
		  }
	   }


 }
 for (int i =0; i < quantPlg; i++)
  urplgArr[i].subtractEqualVerts() ;


}

// возвращает:
// 0 - точка вне полигона
// 1 - точка внутри полигона
// 3  - точка p0  лежит на ребре полигона
// 4  - точка  p0  является вершиной полигона

int TTriang::PntInTriangle( TURPointXY v ) // NEW 21.12.2012
{
  int iSq = Signum(calcSq());
  for (int i =0; i < 3; i++)
  {
	TTriang temrTr(m_pVert[i],m_pVert[i + 1], v);
	double tempSq =  temrTr.calcSq();
	if (fabs (tempSq) < TOLRNC )
	{
	  double d0 = dist2(&m_pVert[i],&m_pVert[i + 1]);
	  double d1 = dist2(&m_pVert[i],&v);
	  double d2 = dist2(&v,&m_pVert[i + 1]);
	  if ((d1 > (d0 + TOLRNC)) ||(d2 > (d0 + TOLRNC)) )return 0;
	  else
	   if ((d1 < (d0 - TOLRNC)) && (d2 < (d0 - TOLRNC)) )return 3;
	   else return 4 ;


	}
	else
	if (iSq * Signum (tempSq)< 0 ) return 0;

  }
  return 1;
}


// пересечение 2-х отрезков
//возвращает true если отрезки пересек и false в противном случае
//если true: если отрезки лежат на одной прямой - то  отрезок- пересечение
// определяется по коду *ipCaseType
// если не наодной прямой, то p0 - точка пересечения и *ipCaseType = 24
// опдробное описание приведено вю файле IntersectTwoSegments
 bool TTriang::IntersectTwoSegments(
					const TURPointXY p00,const  TURPointXY p01 //1-ый сегмент
					,const TURPointXY p10,const TURPointXY p11  // 2-ой сегмент
					, TURPointXY *p0  // точка пересечения
					, int *ipCaseType //
					)
{

 try
 {
		TURPointXY pa;;
	   double a0; ;
	   double alf0, b0,b1,by0,by1;

   switch (TURPolygon::IntersectTwoLine(  p00, p10,p11, p01, p0) )
   {
	   case 0: // не пересекаются параллельны

	   *ipCaseType  = -101;
	   return false ;
	 //  break;

	   case 1:  // песекпряаются в одной точке прямые не параллельны
	   if ( !(
		   ((*p0).X <=  max_d(p00.X,p01.X) + TOLRNC  )
		&& ((*p0).X >=  min_d(p00.X,p01.X) - TOLRNC )
		&& ((*p0).X <=  max_d(p10.X,p11.X) + TOLRNC )
		&& ((*p0).X >=  min_d(p10.X,p11.X) - TOLRNC )
		&& ((*p0).Y <=  max_d(p00.Y,p01.Y) + TOLRNC  )
		&& ((*p0).Y >=  min_d(p00.Y,p01.Y) - TOLRNC )
		&& ((*p0).Y <=  max_d(p10.Y,p11.Y) + TOLRNC )
		&& ((*p0).Y >=  min_d(p10.Y,p11.Y) - TOLRNC )
		))
	   {

		  return false;
	   }
	   else
	   {
	   int ip0 = 2;
	   int ip1 =2;
	   if (dist( *p0,p00) <= TOLRNC )ip0 = 0;
		else
		  if(dist(*p0,p01) <= TOLRNC )ip0 = 1;

	   if (dist( *p0,p10) <= TOLRNC )ip1= 0;
		else
		  if(dist(*p0,p11) <= TOLRNC )ip1 = 1;
		  switch(ip0)
		  {   case 0:
				switch (ip1)
				{
				case 0:
				*ipCaseType = 29;
				break;

				case 1:
				*ipCaseType = 30;

				break ;
				case 2:
				*ipCaseType = 26;

				break;

				default:
				break;	;
				}
			  break;
			  case 1:
				switch (ip1)
				{
				case 0:
				*ipCaseType = 31;

				break;
				case 1:
				*ipCaseType = 32;

				break ;
				case 2:
				*ipCaseType = 25;
				break;

				default:
				break;	;
				}

			  break ;
			  case 2:
				switch (ip1)
				{
				case 0:
				*ipCaseType = 28;

				break;
				case 1:
				*ipCaseType = 27;

				break ;
				case 2:
				*ipCaseType = 24;

				break;

				default:
				break;	;
				}

			  break;
			  default:break;
          }

		return true;
	   }
  //	   break;

	   case 2: // прямые накоторых лежат отрезки совпадают

	   pa.X = p01.X- p00.X ;
		pa.Y = p01.Y - p00.Y;
		a0 = sqrt(pa.X * pa.X + pa.Y * pa.Y) ;
	  alf0 = asin(pa.Y/a0);
	   if (pa.X < 0) alf0 = 3.1415926 - alf0 ;


	   // линейное преобразование координат
	 LinTransf(p00.X, p00.Y, // x0, y0 - координаты
										//точки (0,0) новой системы координат (СК)
				 alf0, // alf0 - угол  между осями 0X старой и новой СК
							  // исчисляется против часовой стрелки от 0X1 к 0X2
				p10.X, p10.Y// x,y - координаты точки в старой СК
				, &b0,  &by0  // xTr,  yTr  - координаты той же точки в новой СК
				);
	 LinTransf(p00.X, p00.Y, // x0, y0 - координаты
										//точки (0,0) новой системы координат (СК)
				 alf0, // alf0 - угол  между осями 0X старой и новой СК
							  // исчисляется против часовой стрелки от 0X1 к 0X2
				p11.X, p11.Y// x,y - координаты точки в старой СК
				, &b1,  &by1  // xTr,  yTr  - координаты той же точки в новой СК
				);
   *ipCaseType =  IntersectTwoSegmInLine( a0, b0, b1);
	if( (*ipCaseType%12) < 11 ) return true;
	else return false ;
   //	   break;
	  default: break;
	}
	return true ;
   }
   catch(...)
   {

   }

   return true ;
}
void TTriang :: InsertFalse(bool *barrOut,const int i0 ,const bool b)
{
	memcpy(&barrOut[i0 +1],&barrOut[i0], (9 - i0)* sizeof(bool)) ;
	barrOut[i0] = false;
}
// вставить вершину p вслед за вершиной i0 - то есть вершину i0+1
void TTriang ::InsertVert(TURPolygon *plg, const int i0 ,const  TURPointXY p)
{
  memcpy( (&(*plg).Points[i0 + 2])   ,&((*plg).Points[i0 + 1]), (8 - i0)* sizeof( TURPointXY )) ;
  (* plg).Points[i0 + 1] = p ;
}
void TTriang ::InsertNum(int *iarrNums,const int i,const int j)
{
  memcpy( &iarrNums[ i + 2 ],&iarrNums[ i + 1 ], (10 - i)* sizeof(int )) ;
  iarrNums[i + 1] = j;
}
int TTriang ::urSign(double a)
{
	int  val = 0;
	if(a > 0) val = 1;
	else
	if (a <0)val = -1;

return val;
}
// описание приве"ено в файле "описание IntersectTwoSegmInLine .docx"
int TTriang ::IntersectTwoSegmInLine(double a0,double b0,double b1)
{
  int i0 ;
  int i1 ;
  switch(urSign(b1 - b0) )
  {
	  case -1: // направления встречные
		if ((b1 > a0 + TOLRNC) || ( b0 < (-TOLRNC)))return 11;
		 i0 = PntInSegment(a0 , b0) ;
		 i1 = PntInSegment(a0 , b1) ;
		 switch(i1)
		 {
			 case 0:  // b1 внутри отрезка
			   switch(i0)
			   {
					 case 0: return 2; // b1 внутри отрезка и b0  внутри отрезка
				   //	 break;
					 case 1: return -1;  // не существ
				  //	 break;
					 case 2: return 5;  //b1 внутри отрезка b0= a0
				  //	 break;
					 case 3: return 0;  //b1 внутри отрезка b0> a0
				  //	 break;
					 default:
					 break;

			   }
			 break;
			 case 1:  // b1 = 0
			   switch(i0)
			   {
					 case 0: return 6;// b1 = 0   a0 > b0 > 0
				   //	 break;
					 case 1: return -1;  // не существ
				   //	 break;
					 case 2: return 4; // b1 = 0  b0 = a0
				  //	 break;
					 case 3: return 8;// b1 = 0  b0 > a0
				   //	 break;
					 default:
					 break;

			   }

			 break;
			 case 2:  // b1 = a0
			   switch(i0)
			   {
					 case 0: return -1;   // не существ
				   //	 break;
					 case 1: return -1;   // не существ
				   //	 break;
					 case 2: return -1;   // не существ
					// break;
					 case 3: return 9;  // b1 = a0 b0 > a0
				  //	 break;
					 default:
					 break;

			   }

			 break;
			 case 3:  // b1 не принадлежит отрезку
			   switch(i0)
			   {
					 case 0: return 1; // b1<0 не принадлежит отрезку bo принадлежит
				   //	 break;
					 case 1: return 10;// b1 <0 не принадлежит отрезку bo = 0 принадлежит
				   //	 break;
					 case 2: return 7; // b1 <0 не принадлежит отрезку bo = a0 принадлежит
					// break;
					 case 3: return 3; //  b1 <0 b2 > a0
				   //	 break;
					 default:
					 break;

			   }

			 break;
			 default:
			 break;
		 }

	  break;
	  case 1: // направления попутные   РЕДАКТИРОВАТЬ!!
		if ((b1 <(- TOLRNC)) || ( b0 > ( a0 + TOLRNC)))return 23;
		 i0 = PntInSegment(a0 , b0) ;
		 i1 = PntInSegment(a0 , b1) ;
		 switch(i1)
		 {
			 case 0: // b1 внутри сегмента
			   switch(i0)
			   {
					 case 0: return 14;// b1 внутри сегмента b0 внутри сегм
				  //	 break;
					 case 1: return 17; // не существ
				  //	 break;
					 case 2: return 17; //b1 внутри сегмента b0 = a0
				 //	 break;
					 case 3: return 12;  //b1 внутри сегмента b0 <0
				   //	 break;
					 default:
					 break;

			   }
			 break;
			 case 1:  // b1 = 0
			   switch(i0)
			   {
					 case 0: return -1; //
			//	 break;
					 case 1: return -1;
			 //		 break;
					 case 2: return -1;
			   //		 break;
					 case 3: return 21; // b1 = 0 b0<0
			   //		 break;
					 default:
					 break;

			   }

			 break;
			 case 2:   // b1 = a0
			   switch(i0)
			   {
					 case 0: return 18;
			 //		 break;
					 case 1: return 16;
			 //		 break;
					 case 2: return -1;
			 //		 break;
					 case 3: return 20;
			  //		 break;
					 default:
					 break;

			   }

			 break;
			 case 3:
			   switch(i0)
			   {
					 case 0: return 13;
			  //		 break;
					 case 1: return 19;
			   //		 break;
					 case 2: return 22;
			   //		 break;
					 case 3: return 15;
			   //		 break;
					 default:
					 break;

			   }

			 break;
			 default:
			 break;
		 }

	  break;
	  case 0:return -2;
 //	  break;
	  default:
	  break;
  }
  return -1001;
}
//точка х в отрезке [0;a]
// возвращпет
// 0 - если в сегменте
// 1 - если x == 0
// 2- если x == 1
// 3 - если не принадлежит сегменту
int TTriang ::PntInSegment(double a,double x)
{
   double b =  (x+TOLRNC)*(x-a -TOLRNC);
  if(b > 0)return 3;
  else
	if ( (x - TOLRNC)*(x-a + TOLRNC)< 0 ) return 0;
	else
	  if (fabs(x) < TOLRNC) return 1;
	  else return 2;
}
double TTriang::dist(TURPointXY p0 ,TURPointXY p1)
{
	return (double)sqrt( (p0.X - p1.X)*(p0.X - p1.X)+ (p0.Y - p1.Y)*(p0.Y - p1.Y)) ;
}
bool TTriang::IsTrianglesIntersectsLine(TTriang *ptrT0,const int i0,TTriang *ptrT1,const int i1)
{
  TTriang tTr( (*ptrT0).m_pVert[(i0 + 2) % 3],(*ptrT1).m_pVert[i1],(*ptrT1).m_pVert[i1 + 1]);
  if(tTr.calcSq()*(*ptrT0).calcSq() >0 )return false;

 return true;
}

// Вычитание треугольника из треугольника trT0\ trT1. trT1 -положительный, trT1 - отрицательный
//  trT1 лежит внутри trT0 не касаясь его границ
void TTriang ::SubTwoEnclosedTrToRings( TTriang *trT0, TTriang *trT1, TURPolygon *urplgArr )
{
	int numPlg = -1;
	int iarrNumEdge[2];
	(*trT0).CutTriangle((*trT1).m_pVert[0],(*trT1).m_pVert[1],iarrNumEdge,urplgArr);
	urplgArr[0].WriteToASCII(L"D:\\GeoProcC++\\FILES\\urplgArr0_1.txt");
	urplgArr[1].WriteToASCII(L"D:\\GeoProcC++\\FILES\\urplgArr1_1.txt");
	if (urplgArr[0].PntInPolygon((*trT1).m_pVert[2]) == 1) numPlg = 0;
	else numPlg = 1 ;
	TURPolygon plgTemp(8);
	memcpy(plgTemp.Points,urplgArr[numPlg].Points,urplgArr[numPlg].NumPoints * sizeof(TURPointXY));
	plgTemp.NumPoints = urplgArr[numPlg].NumPoints +3;
	memcpy(	&plgTemp.Points[iarrNumEdge[numPlg] + 4],&plgTemp.Points[iarrNumEdge[numPlg] + 1]
			   ,(urplgArr[numPlg].NumPoints - iarrNumEdge[numPlg] - 1)*sizeof(TURPointXY));
	if ( (fabs((*trT1).m_pVert[0].X - plgTemp.Points[iarrNumEdge[numPlg]].X)
		  + fabs((*trT1).m_pVert[0].Y - plgTemp.Points[iarrNumEdge[numPlg]].Y))
		  <
		  (fabs((*trT1).m_pVert[1].X - plgTemp.Points[iarrNumEdge[numPlg]].X)
		  + fabs((*trT1).m_pVert[1].Y - plgTemp.Points[iarrNumEdge[numPlg]].Y)) )
		  {
		  plgTemp.Points[iarrNumEdge[numPlg] + 1] = (*trT1).m_pVert[0];
		  plgTemp.Points[iarrNumEdge[numPlg] + 2] = (*trT1).m_pVert[2];
		  plgTemp.Points[iarrNumEdge[numPlg] + 3] = (*trT1).m_pVert[1];

		  }
	else
		  {
		  plgTemp.Points[iarrNumEdge[numPlg] + 1] = (*trT1).m_pVert[1];
		  plgTemp.Points[iarrNumEdge[numPlg] + 2] = (*trT1).m_pVert[2];
		  plgTemp.Points[iarrNumEdge[numPlg] + 3] = (*trT1).m_pVert[0];

		  }


	urplgArr[numPlg] = plgTemp;
}
// рассечение треугольника прямой, заданной 2 точками p1,p2, которые лежат внутри треугольника
// прямая рассекает треугольник на 2 полигона purplgRings
// в массивеи  iarrNumEdge возвращается номера общего ребра полигонов purplgRings
void TTriang ::CutTriangle(TURPointXY p0 ,TURPointXY p1, int *iarrNumEdge,TURPolygon *purplgRings )
{
  TURPointXY cutPoints[3];
  // если прямая пересекает сторону i - (i+1) то numPriznak[i] = 1
  // если не пересекает,numPriznak[i] = 0
  // если пересекает вершину i, то  numPriznak[i] = 2
  int numPriznak[3] = {-1,-1,-1};

	if(fabs(p0.X - p1.X) <= TOLRNC)
	{   // прямая вертикальная
	   for (int i = 0; i < 3; i++)
	   {
		  if(fabs(m_pVert[i].X - m_pVert[i + 1].X) <= TOLRNC)
		  {
		  numPriznak[i] = 0 ;
		  continue;
		  }
		  else
		  {
			  if (
					 ( p0.X <= ( max_d(m_pVert[i].X,m_pVert[i + 1].X ) - TOLRNC) )
				  && ( p0.X >= ( min_d(m_pVert[i].X,m_pVert[i + 1].X ) + TOLRNC) )
				  )
			  {
				double k1 = (m_pVert[i + 1].Y- m_pVert[i ].Y )/(m_pVert[i + 1].X- m_pVert[i ].X );
				cutPoints[i] =  TURPointXY(p0.X, m_pVert[i].Y + k1 * (p0.X - m_pVert[i].X));
				numPriznak[i] = 1 ;

			  }
			  else
			  {
				if (fabs( m_pVert[i].X - p0.X) <= TOLRNC  )
				{
				  numPriznak[i] = 2 ;
				  cutPoints[i]  =  m_pVert[i];
				}
				continue;
			  }
		  }
	   }
	}
	else  // прямая не вертикальная
	{

	  double  k0 = (p1.Y - p0.Y)/(p1.X - p0.X);
	 for (int i = 0; i < 3; i++)
	{


    // пересечение линии и отрезка(сегмента)
// возвращает:
// 0 - не пересекает
// 1 - пересекает
// 2  - прямая проходит через вершину сегмента
// 3  - точка p0,    лежит внутри сегмента
// 4  - точка  p0  является вершиной сегмента
// 5  - сегмент лежит на линии
		switch( TURPolygon::LineCutSegment(m_pVert[i], m_pVert[i + 1]// сегмент
							, p0, k0 // прямая k0 - тангенс наклона
							,cutPoints[i] // точка пересечения
								)  )
      {
		  case 0:
		   numPriznak[i] = 0;
		  continue;
   //		  break;
		  case 1:
		  numPriznak[i] = 1;
   //		  break;
		  case 2:
		   numPriznak[i] = 2;
  //		  break;
		  default:
		  continue;
   //		  break;
      }


	 }
	}

	//
	purplgRings[0] = TURPolygon(5);
	purplgRings[1] = TURPolygon(5);
	switch(numPriznak[ 0 ])
	{
		case 0:
			switch(numPriznak[ 1 ])
			{
			//case 0:
		   //	break;
			case 1:
			purplgRings[0].Points[0] = cutPoints[1] ;
			purplgRings[0].Points[1] = cutPoints[2] ;
			purplgRings[0].Points[2] = m_pVert[0] ;
			purplgRings[0].Points[3] = m_pVert[1] ;
			purplgRings[0].Points[4] = cutPoints[1] ;
			iarrNumEdge[0] = 0 ;
			purplgRings[1].Points[0] = cutPoints[2] ;
			purplgRings[1].Points[1] = cutPoints[1] ;
			purplgRings[1].Points[2] = m_pVert[2] ;
			purplgRings[1].Points[3] = cutPoints[2] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;
			break;
			case 2:
			purplgRings[0].Points[0] = m_pVert[1] ;
			purplgRings[0].Points[1] = cutPoints[2] ;
			purplgRings[0].Points[2] = m_pVert[0] ;
			purplgRings[0].Points[3] = m_pVert[1 ];
			purplgRings[0].NumPoints = 4;
			iarrNumEdge[0] = 0 ;
			purplgRings[1].Points[0] = cutPoints[2] ;
			purplgRings[1].Points[1] = m_pVert[1] ;
			purplgRings[1].Points[2] = m_pVert[2] ;
			purplgRings[1].Points[3] = cutPoints[2] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;

			break;
			default:
			break;
			}

		break;
		case 1:
			switch(numPriznak[ 2 ])
			{
			case 0:
			purplgRings[0].Points[0] = cutPoints[0] ;
			purplgRings[0].Points[1] = cutPoints[1] ;
			purplgRings[0].Points[2] = m_pVert[2] ;
			purplgRings[0].Points[3] = m_pVert[0] ;
			purplgRings[0].Points[4] = cutPoints[0] ;
			iarrNumEdge[0] = 0 ;
			purplgRings[1].Points[0] = cutPoints[1] ;
			purplgRings[1].Points[1] = cutPoints[0] ;
			purplgRings[1].Points[2] = m_pVert[1] ;
			purplgRings[1].Points[3] = cutPoints[1] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;

			break;
			case 1:
			purplgRings[0].Points[0] = cutPoints[2] ;
			purplgRings[0].Points[1] = cutPoints[0] ;
			purplgRings[0].Points[2] = m_pVert[1] ;
			purplgRings[0].Points[3] = m_pVert[2] ;
			purplgRings[0].Points[4] = cutPoints[2] ;
			iarrNumEdge[0] = 0 ;
			purplgRings[1].Points[0] = cutPoints[0] ;
			purplgRings[1].Points[1] = cutPoints[2] ;
			purplgRings[1].Points[2] = m_pVert[0] ;
			purplgRings[1].Points[3] = cutPoints[0] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;

			break;
			case 2:
			purplgRings[0].Points[0] = cutPoints[0] ;
			purplgRings[0].Points[1] = m_pVert[2] ;
			purplgRings[0].Points[2] = m_pVert[1] ;
			purplgRings[0].Points[3] = cutPoints[0] ;
			iarrNumEdge[0] = 0;
			purplgRings[0].NumPoints = 4;
			purplgRings[1].Points[0] = m_pVert[2] ;
			purplgRings[1].Points[1] = cutPoints[0] ;
			purplgRings[1].Points[2] = m_pVert[1] ;
			purplgRings[1].Points[3] = m_pVert[2] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;

			break;
			default:
			break;
			}

		break;
		case 2:
			purplgRings[0].Points[0] = cutPoints[1] ;
			purplgRings[0].Points[1] = m_pVert[0] ;
			purplgRings[0].Points[2] = m_pVert[1] ;
			purplgRings[0].Points[3] = cutPoints[1] ;
			iarrNumEdge[0] = 0;
			purplgRings[0].NumPoints = 4;
			purplgRings[1].Points[0] = m_pVert[0] ;
			purplgRings[1].Points[1] = cutPoints[1] ;
			purplgRings[1].Points[2] = m_pVert[2] ;
			purplgRings[1].Points[3] = m_pVert[0] ;
			iarrNumEdge[1] = 0 ;
			purplgRings[1].NumPoints = 4;


		break;
		default:
		break;
	}

}

void TTriang::ChangeDir()
{
	TURPointXY pointT = m_pVert[2];
	m_pVert[2] = m_pVert[1];
	m_pVert[1] = pointT;
}

// если 2 треугольника пересекаются по стороне,
// то каждый из них надо разбить на 2 треугольника
//То есть, треугольники из массива *ppTr могут пересекаться или
// только по единственной вершине или по стороне
// напнример треугольники (0,0),(2,2),(4,0),(0,0)  и  (2,0)(6,0),(4,-2),(2,0)
//пересекаются по отрезку (4,0),(6,0)
// каждый из этих треугольн7иков надо разбить на 2.
// тка первый треуг надор разбить так: треуг 1-1 :
// (0,0),(2,2),(2,0),(0,0)
// треуг 1-2:   (2,0),(2,2),(4,0),(2,0)

//  *lenTrArr - на входе начальная длина массива *ppTr
// на выходе - длина  массива *ppTr
void TTriang::FinishDisassembleTriangleArr(TTriang **ppTr,int *lenTrArr)
{
  int lenTemp = *lenTrArr  ;
  int lenMemory  = *lenTrArr  ;
  int iCounter = 0;
  while(1)
  {
	  iCounter++;
	  bool bBreak = false;
	  for (int i = 0; i < lenTemp; i++)
	  {


	  for (int j = 0; j < lenTemp; j++)
	  {
	   //	if( (iCounter == 4) && (i == 4)&&(j == 5))
	   //	{
		//  iCounter = iCounter;
		//}
		if (j == i) continue ;
		TURPointXY urpntP;
		int i0 = -1;
		int i1 = -1;
		int iType =  TypeOfTrianglesInersection( (&(*ppTr)[i]), (&(*ppTr)[j]),&urpntP,&i0,&i1);
	   /*	if  (iType == 0 )
		{
			 if (lenTemp+ 2 > lenMemory )
			 {
			  lenMemory += 100;
			  *ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			 }
			 TTriang trT00((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 + 1]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);
			 TTriang trT01((*ppTr)[j].m_pVert[i1 + 1],(*ppTr)[i].m_pVert[i0 +1 ]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );
			 (*ppTr)[i] = trT00;
			 (*ppTr)[lenTemp] = trT01;

			 TTriang trT10((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 + 1]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 TTriang trT11((*ppTr)[i].m_pVert[i0 + 1],(*ppTr)[j].m_pVert[i1 + 1 ]
					,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 (*ppTr)[j] = trT10;
			 (*ppTr)[lenTemp + 1] = trT11;
			 lenTemp = lenTemp+ 2 ;
			bBreak = true;
		}
		if  (iType == 1 )
		{
			 if (lenTemp+ 2 > lenMemory )
			 {
			  lenMemory += 100;
			  *ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			 }
			 TTriang trT00((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 ]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);
			 TTriang trT01((*ppTr)[j].m_pVert[i1 ],(*ppTr)[i].m_pVert[i0 +1]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );

			 TTriang trT10((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 ]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 TTriang trT11((*ppTr)[i].m_pVert[i0 ],(*ppTr)[j].m_pVert[i1 + 1 ]
					,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 (*ppTr)[j] = trT10;
			 (*ppTr)[lenTemp + 1] = trT11;
			 (*ppTr)[i] = trT00;
			 (*ppTr)[lenTemp] = trT01;

			 lenTemp = lenTemp+ 2 ;
			bBreak = true;
		}   */
		TTriang trT00,trT01,trT02,trT10,trT11,trT12 ;
		 switch (iType)
		 {
		  case 0:
			 if (lenTemp+ 2 > lenMemory )
			 {
			  lenMemory += 100;
			  *ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			 }
			 trT00 = TTriang((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 + 1]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);
			 trT01= TTriang((*ppTr)[j].m_pVert[i1 + 1],(*ppTr)[i].m_pVert[i0 +1 ]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );
			 (*ppTr)[i] = trT00;
			 (*ppTr)[lenTemp] = trT01;

			 trT10= TTriang((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 + 1]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			  trT11= TTriang((*ppTr)[i].m_pVert[i0 + 1],(*ppTr)[j].m_pVert[i1 + 1 ]
					,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 (*ppTr)[j] = trT10;
			 (*ppTr)[lenTemp + 1] = trT11;
			 lenTemp = lenTemp+ 2 ;
			bBreak = true;

		  break;

		  case 1:
			 if (lenTemp+ 2 > lenMemory )
			 {
			  lenMemory += 100;
			  *ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			 }
			  trT00= TTriang((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 ]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);
			  trT01= TTriang((*ppTr)[j].m_pVert[i1 ],(*ppTr)[i].m_pVert[i0 +1]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );

			  trT10= TTriang((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 ]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			  trT11= TTriang((*ppTr)[i].m_pVert[i0 ],(*ppTr)[j].m_pVert[i1 + 1 ]
					,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);
			 (*ppTr)[j] = trT10;
			 (*ppTr)[lenTemp + 1] = trT11;
			 (*ppTr)[i] = trT00;
			 (*ppTr)[lenTemp] = trT01;

			 lenTemp = lenTemp+ 2 ;
			bBreak = true;

		  break;
		  case 2:
			if (lenTemp+ 2 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT00= TTriang((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 + 1 ]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);

			  trT01= TTriang((*ppTr)[j].m_pVert[i1 + 1 ],(*ppTr)[j].m_pVert[i1 ]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );
			  trT02= TTriang((*ppTr)[i].m_pVert[(i0 -1 + 3)%3],(*ppTr)[j].m_pVert[i1 ]
			 ,(*ppTr)[i].m_pVert[i0 + 1]);
			 (*ppTr)[i] = trT00;
			  (*ppTr)[lenTemp] = trT01;
			  (*ppTr)[lenTemp + 1] = trT02;
			  lenTemp = lenTemp+ 2 ;
			bBreak = true;
		  break;
		  case 3:
			if (lenTemp+ 2 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT10= TTriang((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 + 1 ]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);

			  trT11= TTriang((*ppTr)[i].m_pVert[i0 + 1 ],(*ppTr)[i].m_pVert[i0  ]
			 ,(*ppTr)[j].m_pVert[(i1 - 1 + 3) % 3] );
			   trT12= TTriang((*ppTr)[j].m_pVert[i1 + 1],(*ppTr)[j].m_pVert[(i1 - 1 + 3) % 3]
			 ,(*ppTr)[i].m_pVert[i0 ]);
			 (*ppTr)[j] = trT10;
			  (*ppTr)[lenTemp] = trT11;
			  (*ppTr)[lenTemp + 1] = trT12;
			  lenTemp = lenTemp+ 2 ;
			bBreak = true;

		  break;
		  case 5:
			if (lenTemp+ 1 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT00= TTriang((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1 + 1 ]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);

			  trT01= TTriang((*ppTr)[j].m_pVert[i1 + 1 ],(*ppTr)[i].m_pVert[i0 + 1 ]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );
			 (*ppTr)[i] = trT00;
			  (*ppTr)[lenTemp] = trT01;

			  lenTemp = lenTemp+ 1 ;
			bBreak = true;

		  break;
		  case 6:
			if (lenTemp+ 1 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT00= TTriang((*ppTr)[i].m_pVert[i0],(*ppTr)[j].m_pVert[i1  ]
			 ,(*ppTr)[i].m_pVert[(i0 -1 + 3)%3]);

			  trT01= TTriang((*ppTr)[j].m_pVert[i1 ],(*ppTr)[i].m_pVert[i0 + 1 ]
			 ,(*ppTr)[i].m_pVert[(i0 - 1 + 3) % 3] );
			 (*ppTr)[i] = trT00;
			  (*ppTr)[lenTemp] = trT01;

			  lenTemp = lenTemp+ 1 ;
			bBreak = true;
		   break;
		  case 7:
			if (lenTemp+ 1 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT10= TTriang((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 ]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);

			  trT11= TTriang((*ppTr)[i].m_pVert[i0 ],(*ppTr)[j].m_pVert[i1 + 1 ]
			 ,(*ppTr)[j].m_pVert[(i1 - 1 + 3) % 3] );
			 (*ppTr)[j] = trT10;
			  (*ppTr)[lenTemp] = trT11;

			  lenTemp = lenTemp+ 1 ;
			bBreak = true;

		  break;
		  case 8:
			if (lenTemp+ 1 > lenMemory )
			{
			lenMemory += 100;
			*ppTr = (TTriang *)realloc(*ppTr,lenMemory  * sizeof(TTriang));
			}
			trT10= TTriang((*ppTr)[j].m_pVert[i1],(*ppTr)[i].m_pVert[i0 +1 ]
			 ,(*ppTr)[j].m_pVert[(i1 -1 + 3)%3]);

			  trT11= TTriang((*ppTr)[i].m_pVert[i0 +1],(*ppTr)[j].m_pVert[i1 + 1 ]
			 ,(*ppTr)[j].m_pVert[(i1 - 1 + 3) % 3] );
			 (*ppTr)[j] = trT10;
			  (*ppTr)[lenTemp] = trT11;

			  lenTemp = lenTemp+ 1 ;
			bBreak = true;

		  break;
		 default:
		 break ;
            }

		if (bBreak)break;
	  }
	if (bBreak)break;
   }
   if(!bBreak) break;
  }
   *ppTr = (TTriang *)realloc(*ppTr,lenTemp  * sizeof(TTriang));
   *lenTrArr = lenTemp;
}
// 2 положительных треугольникка не имеют общих
// внутренних точек, но могут пересекаться по стороне.
// В этом случае стороны имеют противоположное направление
// если пересечение есть, то в случае 0 пересечения сегментов возвращается 0
// в случае 1 возвращается 1
//  *urpntP - точка пересечения
// в остальнх лучаях возвращается -1
int TTriang::TypeOfTrianglesInersection( TTriang *pTr0, TTriang *pTr1
		,TURPointXY *urpntP, int *i0,int * i1)
{
	int iReturn = -1;
	int ipCaseType =-1;
	for (int i =0; i < 3; i++)
	for (int j =0; j < 3; j++)
	{
	   if(IntersectTwoSegments(
				(*pTr0).m_pVert[i],(*pTr0).m_pVert[i + 1] //1-ый сегмент
				,(*pTr1).m_pVert[j],(*pTr1).m_pVert[j + 1]  // 2-ой сегмент
					, urpntP  // точка пересечения
					, &ipCaseType //
					) )
	   {
		  /* if ((ipCaseType  == 0) || (ipCaseType  == 1))
		   {
			*i0 = i;
			*i1 = j;
			return ipCaseType;
		   }
			  */
		 switch (ipCaseType)
		 {
		  case 0:

		  case 1:

		  case 2:

		  case 3:

		  case 5:

		  case 6:

		  case 7:

		  case 8:
		  *i0 = i;
		  *i1 = j;
		  return ipCaseType;
	  //	  break;
		 default:
		 break ;
		 }

	   }
	}
	return iReturn;
}
#pragma package(smart_init)
