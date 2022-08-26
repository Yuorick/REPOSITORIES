


#pragma hdrstop
#include <vcl.h>
#include "ShipTraj.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "MatrixProccess.h"
#include <time.h>
#include "UrPointXY.h"
#include "HomingSituation.h"
#include "URPolygon.h"

// ДЛЯ ДАТЧИКА СЛУЧ ЧИСЕЛ
#define DEBUG



//---------------------------------------------------------------------------
TShipTraj::TShipTraj()
{
	miType    = 0;
	mTCur = 5;
	mTBegin = 0;
	memset(marrParams,0,25 * sizeof(long double  )) ;
	marrParams[24] = 0.5 ;
	mSigW = 33. ;

	memset(marrVectSostGSK, 0, 9*sizeof(long double  )) ;
	memset(marrVectDeviations, 0, 9*sizeof(long double  )) ;
	memset(marrVectSostGSK_Begin, 0, 9*sizeof(long double  )) ;




}

//---------------------------------------------------------------------------


// конструктор копирования
 TShipTraj ::TShipTraj (const TShipTraj &R)
 {
	mTCur = R.mTCur ;
	miType = R.miType;
	mTBegin = R.mTBegin;
	mT = R.mT;
	mSigW = R.mSigW;
	memcpy(marrParams,R.marrParams, 25 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK,R.marrVectSostGSK, 9 * sizeof(long double  )) ;
	memcpy(marrVectDeviations,R.marrVectDeviations, 9 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK_Begin,R.marrVectSostGSK_Begin, 9 * sizeof(long double  )) ;


 }
 // оператор присваивания
 TShipTraj TShipTraj::operator=(TShipTraj  R)
 {

	mTCur = R.mTCur ;
	miType = R.miType;
	mTBegin = R.mTBegin;
	mT = R.mT;
	mSigW = R.mSigW;
	memcpy(marrParams,R.marrParams, 25 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK,R.marrVectSostGSK, 9 * sizeof(long double  )) ;
	memcpy(marrVectDeviations,R.marrVectDeviations, 9 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK_Begin,R.marrVectSostGSK_Begin, 9 * sizeof(long double  )) ;



	return *this ;
 }

  // парам конструктор
 TShipTraj::TShipTraj (const long double   TCur, const int  iType , const long double   	SigW
		  ,long double   *arrParams 	,long double   *arrVectDeviations,long double   *arrVectSostGSK
		  ,long double   *arrVectSostGSK_Begin )
 {


	mTCur = TCur;
	miType = iType ;
	mSigW = SigW ;
	memcpy(marrParams,arrParams, 25 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK,arrVectSostGSK, 9 * sizeof(long double  )) ;
	memcpy(marrVectDeviations,arrVectDeviations, 9 * sizeof(long double  )) ;
	memcpy(marrVectSostGSK_Begin,arrVectSostGSK_Begin, 9 * sizeof(long double  )) ;

 }
  // парам конструктор
  // PartData - член класса, содержащий параметры, описывающие траекторию   цели
  // arrVSTarget -  вектор состояния цели, 9 компонентов
  // TCur -  текущее время

 TShipTraj::TShipTraj( TPartData PartData,long double   *arrVSTarget, const long double   TCur)
{
   miType = PartData.miTypePart;
   mTCur = TCur ;
   mTBegin = TCur;
   mT = PartData.mTimePart;
   mSigW = PartData.mSigW ;
 // memset(marrVS_OSKPSb, 0, 6*sizeof(long double  )) ;
  memset(marrParams, 0, 25 *sizeof(long double  )) ;
  marrParams[24] = 0.5;// первый участок траектории нужен только для горки

  memcpy(marrVectSostGSK, arrVSTarget, 9 * sizeof(long double  )) ;
  memset(marrVectDeviations, 0, 9 * sizeof(long double  )) ;  // исправить !!!
  memcpy(marrVectSostGSK_Begin, arrVSTarget, 9 * sizeof(long double  )) ;
// формирование массива  arrParams для трпаекторий различных типов
  if(miType == 1)
  {

  }

  if(miType == 2)   //"кабрирование"
  {
// формирование массива  marrParams для трпаектории типа "Кабрирование":
//marrParams[0],marrParams[1],marrParams[2],marrParams[3] ,marrParams[4] - координаты центра  окружномсти , радиус, продолжительность
//marrParams[5],marrParams[6],marrParams[7],marrParams[8],marrParams[9] - нули
//marrParams[10],marrParams[11],marrParams[12],marrParams[13],marrParams[4] - нули
// marrParams[15] - marrParams[23] - матрица перехода из СКПМ в ГСК
// СКПМi - система координат плоскости маневра. Центр находится в центре  окружности
// Ось X направлена по вектору и направленному из центра окружности к начальной точке
// Ось Y направлена по начальной скорости
// Ось Z направлена по нормали к плоскости
// marrParams[24 ] - нуль
// В частности arrParams[21],arrParams[22],arrParams[23] - вектор ед нормали к плоскости менвра
// описание  TPartData
// PartData.marrData[0]  - угол плоскости маневра между нормальью к плоскости и вертикалью

// PartData.marrData[1],PartData.marrData[2]- время, ускорение
// PartData.marrData[3],PartData.marrData[4]- нули
// PartData.marrData[5],PartData.marrData[6]- нули
		long double   v2 =   norm(&arrVSTarget[3])* norm(&arrVSTarget[3]); // v2 = v*v
		long double   valW = PartData.marrData[ 2 ] ;
		long double   valR =  v2 / valW ;             // радиус окружности
	   //	long double   valTime  = PartData.marrData[ 1 ] ;  // время участка - чтобы удобней писать текст
		long double   angGam = PartData.marrData[0] ;    // угол наклона траектории (нормали плоскости маневра и вертикали

		// вектор ед нормали к плосксти маневра
		long double   arrN[3] ={0} ;
	if (!FindNormVect(arrVSTarget, angGam, arrN) )
	{
	 ShowMessage(L"Ошибка в конструкторе траектории кабрирования ") ;
	}

		// нахождение центра окружности arrCentre

		FindCircleCentre(arrN, arrVSTarget,  valR, marrParams) ;



		CreateMtrxTansf_SKPM_TO_GSK(arrVSTarget, marrParams,arrN, &marrParams[15]); // матрица преобразования заполнена
		marrParams[ 3 ]  =  valR;
		marrParams[ 4]  =  PartData.mTimePart;

  }

  if(miType == 3)   //"равноускоренное криволинейное движение"
  {
   // формирование массива  marrParams для трпаектории типа "равноускоренное криволинейное движение":
		//marrParams[0],marrParams[1],marrParams[2],- ускорения по осям ГСК
		 // описание  TPartData
		 // PartData.marrData[0],PartData.marrData[1],PartData.marrData[2]  - ускорения по осям ГСК

		// огстальные нули
		memset (marrParams, 0, 25 * sizeof(long double  ) ) ;
		memcpy(marrParams, PartData.marrData, 3 * sizeof(long double  )) ;
  }


}

 void  TShipTraj::CreateMtrxTansf_SKPM_TO_GSK(long double   *arrS, long double   *arrCentre,long double   *arrN,long double   *arrOut)
 {
	 // ось x  SKPM
	 long double   arrT[3] = {0},arrT0[3] = {0};
	 MtrxMinusMatrx(arrS, arrCentre,3, 1, arrT) ;
	 long double   valn = 1./ norm(arrT);
	 MatrxMultScalar(arrT, 3, 1, valn, arrT0) ;
	 for (int i = 0; i < 3 ; i++) arrOut [i * 3] = arrT0[i] ;
	// Ось Y  SKPM
	  long double   valv = 1./ norm(&arrS[3]);
	  MatrxMultScalar(&arrS[3], 3, 1, valv, arrT0) ;
	  for (int i = 0; i < 3 ; i++) arrOut [i * 3 + 1] = arrT0[i] ;
	 // Ось Z
	 for (int i = 0; i < 3 ; i++) arrOut [i * 3 + 2] = arrN[i] ;
 }

void TShipTraj::FindFinalPoint_SKPM(const long double   valV,const long double   valR, long double   valT, long double   *arrSFin)
{
  long double   valOm =  valV/valR ;
  arrSFin[0] =  valR* cos(valOm * valT) ;
  arrSFin[1] =  valR* sin(valOm * valT) ;
  arrSFin[2] =   0;
  arrSFin[3] =   -valV * sin(valOm * valT) ;
  arrSFin[4] =   valV * cos(valOm * valT) ;
  arrSFin[5] =   0 ;
}
//
bool TShipTraj::FindNormVect(long double   *arrS,const long double   angGam, long double   *arrN)
{

 long double   *parrV =  &arrS[3] ;


 // нахождение нормали к плоскости

 if ( fabs(parrV[1]) > 0.00001 )
 {
   long double   sq2 =  (1. + parrV[0] *parrV[0]/ parrV[1]/ parrV[1]) * sin(angGam)* sin(angGam)
		- parrV[2] *parrV[2]/ parrV[1]/ parrV[1] * cos(angGam)* cos(angGam) ;
	if ( sq2 < 0)
	{
	  ShowMessage(L"Ошибка в задании траектории Горка void TShipTraj::FindNormVect 1") ;
	  return false;
	}
	arrN[0] = ( -parrV[0] * parrV[2] * cos (angGam)/ parrV[1]/ parrV[1] + Sign(angGam)* sqrt(sq2) )/(1. + parrV[0] *parrV[0]/ parrV[1]/ parrV[1]);
	arrN[1] =  -parrV[2]/ parrV[1]* cos(angGam)- parrV[0]/ parrV[1]* arrN[0] ;
	arrN[2] = cos (angGam) ;

}
else
{
		 if (fabs(parrV[0]) > 0.00001)
		 {
		   arrN[2] = cos (angGam) ;
		   arrN[0] = - parrV[2] * cos(angGam) / parrV[0];
		   long double   sq12 = 1. - arrN[2] * arrN[2] -  arrN[0] * arrN[0];
		   if (fabs(sq12) < 0.000001)
		   {
			sq12 = 0.;

		   }
		   else
		   {
		   if ( sq12 < 0)
			{
			  ShowMessage(L"Ошибка в задании траектории Горка void TShipTraj::FindNormVect 2") ;
			  return false;
			}
			}
		   arrN[1] = Sign(angGam)* sqrt (sq12);
		 }
		 else
		 {
		  if ((fabs (angGam ) - M_PI/2.) > 0.00001)
		  {
			  ShowMessage(L"Ошибка в задании траектории Горка void TShipTraj::FindNormVect 3") ;
			  return false;
		  }
		  else
		  {
		  arrN[0] = 1.;
		  arrN[1] = 0.;
		  arrN[2] = 0. ;
		  return false;
		  }
		 }
}
return true ;
}

int TShipTraj::Sign(const long double   a)
{
	if (fabs ( a) < 0.00000001) return 0;


	return (a >= 0)?1:-1;
}

// нахождение центра вращения
// INPUT:
// arrN - ед нормаль к плоскости
// arrS - положение и скорость точки
//  valR - радиус
// idir - направление вращения idir = 1, если положительное и idir = -1, если отрицат
//OutPut
//  arrOut - центр вращения

void  TShipTraj::FindCircleCentre(long double   *arrN,long double   *arrS, const long double   valR, /*int idir,*/ long double   *arrOut)
{
  long double   arrT [3] ={0} ,arrT0 [3] ={0} ;
  VectMult(arrN, &arrS[3], arrT) ;
  long double   valD = valR/norm(arrT)/**idir*/ ;
  MatrxMultScalar(arrT, 3, 1, valD,arrT0) ;
  MtrxSumMatrx(arrS, arrT0,3, 1, arrOut) ;

}

long double    TShipTraj::norm(long double   *arr)
{
	return sqrt(arr[0]* arr[0] + arr[1]* arr[1] +arr[2]* arr[2]) ;
}

void  TShipTraj::VectMult(long double   *px, long double   *py, long double   *arrRez)
{
  arrRez [0] = px[1] * py[2] - px[2] * py[1];
  arrRez [1] =  -px[0] * py[2] + px[2] * py[0];
  arrRez [2] = px[0] * py[1] - px[1] * py[0];
}

//
bool TShipTraj::recalcTrajPoint(const long double   tNext)
{
	switch(miType)
	{
		case 0:
		recalcTrajPoint_0( tNext) ; // равногмерное
		break ;

		break;
		case 2:
		recalcTrajPoint_2( tNext) ;  // кабрирование
		break;
		case 3:
		recalcTrajPoint_3( tNext) ;  // равноускоренное криволин движение
		break;
		default:
		return false ;
	}
	return true ;
}
// случай равномерного движения
void TShipTraj::recalcTrajPoint_0(const long double   tNext)
{
  long double   h = tNext - mTCur;
  if ( h < 0)
  {
	 ShowMessage(L"Error.TShipTraj::recalcTrajPoint_0. h<0") ;
	 return ;
  }
  long double   arrTemp[9] = {0} ;
  for (int i = 0; i < 3; i++) arrTemp[6 + i] = getGauss(0, mSigW) ;
  for (int i = 0; i < 3; i++)
  {
	 arrTemp[3 + i] = marrVectSostGSK  [3 + i]+ h *arrTemp[6 + i];
	 arrTemp[    i] = marrVectSostGSK  [ i] + marrVectSostGSK  [3 + i]* h
				+ arrTemp[6 + i]* arrTemp[6 + i]*h*h/2. ;;
  }
  memcpy( marrVectSostGSK, arrTemp, 9 * sizeof(long double  )) ;
  mTCur = tNext;


}


// случай маневра Кабрирование
void TShipTraj::recalcTrajPoint_2(const long double   tNext)
{
  long double   h = tNext - mTCur;
  if ( h < 0)
  {
	 ShowMessage(L"Error.TShipTraj::recalcTrajPoint_0. h<0") ;
	 return ;
  }


	// время от начачла текущего участка
	 long double   tPart = tNext - mTBegin;

	 //

	 // расчет точки вектора состояния в СКПМi
	 // arrVectSostGSK_Ideal - вектор состояния без добавления возмущений
	 // marrVectSostGSK - вектор состояния с добавлеными возмущениями
	 long double   valV = norm(&marrVectSostGSK_Begin[3]);
	 long double   valR = marrParams[3 ];
	 long double   valOm =  valV/valR ;
	 long double   arrSSKPM[9] = {0};
	 long double   arrMtrxTrans[4] = {0}, arrX[2] = {0}, arrV[2] = {0}, arrW[2]={0} ;
	 arrX[0] = valR ;  arrX[1] = 0 ;
	 arrV[0] = 0; arrV[1] = valV;
	 arrW[0] = -valV*valV/valR ;  arrW[1] = 0 ;
	 createTransfMtrx(valOm * tPart, arrMtrxTrans) ;

	 MtrxMultMatrx(arrMtrxTrans,2, 2, arrX,1, arrSSKPM);
	 MtrxMultMatrx(arrMtrxTrans,2, 2, arrV,1, &arrSSKPM[3]);
	 MtrxMultMatrx(arrMtrxTrans,2, 2, arrW,1, &arrSSKPM[6]);
	arrSSKPM[2] =  0 ;
	arrSSKPM[5] =  0;
	arrSSKPM[8] =  0 ;



	 // сюда вставить добавление возмущений к  arrSSKPM
	 long double   arrDevTemp[9] ={0} ; // вектор возмущений в СК связанной с целью
	 calcAddDeviations_Gorka(h, arrDevTemp);
	 //
	 long double   arrDevSKPM[9] ={0} ; // вектор возмущений пересчитанный в CКПМ
	 memcpy(arrDevSKPM, arrDevTemp, 9 * sizeof(long double  )) ;
	 MtrxMultMatrx(arrMtrxTrans,2, 2, arrDevTemp,1, arrDevSKPM);
	 MtrxMultMatrx(arrMtrxTrans,2, 2, &arrDevTemp[3],1, &arrDevSKPM[3]);
	 MtrxMultMatrx(arrMtrxTrans,2, 2, &arrDevTemp[6],1, &arrDevSKPM[6]);
	 // добавление возмущений к arrSSKPM
	 long double   arrSTemp[9] = {0};
	 MtrxSumMatrx(arrDevSKPM, arrSSKPM,9, 1, arrSTemp) ;
	 memcpy(arrSSKPM,arrSTemp , 9 * sizeof(long double  )) ;

	 ///
	 // персесчет  вектора состояния в КГСК
	  TrasformSSKPM_TO_GSK(0,arrSSKPM, marrVectSostGSK);

	 mTCur =  tNext ;
}

// случай маневра Равноускоренное криволинейное движение
void TShipTraj::recalcTrajPoint_3(const long double   tNext)
{
  long double   h = tNext - mTCur;
  long double   delT =  tNext - mTBegin ;
  if ( h < 0)
  {
	 ShowMessage(L"Error.TShipTraj::recalcTrajPoint_3. h<0") ;
	 return ;
  }


	// время от начачла текущего участка
	 long double   tPart = tNext - mTBegin;

	 //  Нахождение возмущений по 3- оясям ГСК
	  long double   arrDevTemp[9] ={0} ; // вектор возмущений в СК связанной с целью

	 memcpy(marrVectDeviations, arrDevTemp, 9 * sizeof(long double  )) ;

	 long double   arrVectSost[9] ={0.} ; // вектор состояния идеальный
	 for (int i = 0; i <3; i++)
	 {
	   arrVectSost[ i] = marrVectSostGSK_Begin[  i] + delT * marrVectSostGSK_Begin[ 3 + i]
			 + delT *delT *marrParams[i]/2;
	   arrVectSost[ 3 + i] = marrVectSostGSK_Begin[ 3 + i] +  delT *marrParams[i];
	   arrVectSost[ 6 + i] =  marrParams[i];
	 }
	 MtrxSumMatrx(arrVectSost, marrVectDeviations,9, 1, marrVectSostGSK) ;

	 mTCur =  tNext ;
}


 void TShipTraj::createTransfMtrx(const long double   alf, long double   *arrMtrxOut)
 {
  memset(arrMtrxOut,0,4 * sizeof(long double  )) ;
  arrMtrxOut[0] = cos (alf) ;
  arrMtrxOut[1] = -sin (alf) ;
  arrMtrxOut[2] = -arrMtrxOut[1];
  arrMtrxOut[3] =  arrMtrxOut[0];

 }

 void TShipTraj::calcAddDeviations_Gorka(const long double   h, long double   *arrOut)
 {
   // al, be - параметры динам системы 2-го порядка, формирующей возмущение по ускорению
   // x1(t+h) =  x1(t) + h * x2(t) + t*t/2 * w
   // x2(t+h) = -al * h *x1(t) + (1 + 2 * be * h) * x2(t) + (h + h*h * be)* w
   const long double   al = 6;
   const long double   be = -2;
   // arrA - матрица этой системы
  long double   arrA[4] = {0};
  arrA[0] = 1.;
  arrA[1] = h ;
  arrA[2] = -al * h ;
  arrA[3] =  1. + 2 * be * h ;
  // arrB - матрица возмущения
  long double   arrB[2] = {0} ;
  arrB [0] =  h * h /2. ;
  arrB [1] =  h + be * h * h ;
  //
  long double   arrNoise[3] = {0.};
  arrNoise[0] = getGauss( 0, 5.5) ;
  arrNoise[1] = -0.21 * arrNoise[0] ;
  arrNoise[2] = getGauss( 0, mSigW) ;
  //
   for (int i =0; i < 3; i++)
   {
	 recalcSys2(arrA,arrB,arrNoise[i],&marrVectDeviations[i * 3]);
   }

   reformVectDev(marrVectDeviations,arrOut);


 }

 void TShipTraj::reformVectDev(long double   *arrInp,long double   *arrOut)
 {
  arrOut[0] = arrInp[0] ;
  arrOut[1] = arrInp[3] ;
  arrOut[2] = arrInp[6] ;
  arrOut[3] = arrInp[1] ;
  arrOut[4] = arrInp[4] ;
  arrOut[5] = arrInp[7] ;
  arrOut[6] = arrInp[2] ;
  arrOut[7] = arrInp[5] ;
  arrOut[8] = arrInp[8] ;
 }
 void TShipTraj::recalcSys2(long double   *arrA, long double   *arrB, const long double   valNoise,long double   *arrVectDeviations)
 {
	long double   arrT0[2]= {0},arrT1[2]= {0};
	MtrxMultMatrx(arrA, 2, 2, arrVectDeviations,1, arrT0) ;
	MatrxMultScalar(arrB, 2, 1, valNoise,arrT1) ;
	MtrxSumMatrx(arrT0, arrT1,2, 1, arrVectDeviations) ;
	arrVectDeviations[2] = valNoise ;

 }
// пересчет вектора состояния из СКПМ в КГСК
 void TShipTraj::TrasformSSKPM_TO_GSK(const int numPart,long double   *arrSSKPM, long double   *arrSGSK)
 {
	for (int i = 0; i < 3; i++)
	{
	MtrxMultMatrx(&marrParams[15],3, 3, &arrSSKPM[3*i],1, &arrSGSK[ 3 * i]) ;
	}

	long double   arrT0[3] = {0},arrT1[3] = {0},arrT2[3] = {0},arrT3[3] = {0} ;
	MtrxSumMatrx(arrSGSK, &marrParams[5 *numPart],3, 1,  arrT0);
	memcpy(arrSGSK,arrT0, 3 * sizeof(long double  )) ;
 }
// вычисление номера участка траектории для момента воемени valT
int TShipTraj::getPartNum(const long double   valT)
{
  if (valT <= mTBegin ) return -2;
  long double   tTemp = mTBegin ;
  for (int i = 0; i < 3; i++)
  {
	tTemp += marrParams[4 + 5*i]  ;
	if (valT <= (tTemp - 0.00001)) return i;

  }

  return -1 ;
}


long double   TShipTraj::getGauss(const long double   a, const long double   sig )

{
  long double   val = 0;
  int num = 12;
 // for (int i = 0; i < num; i++) val += ( (long double  )rand())/((long double  )RAND_MAX) - 0.5;
  for (int i = 0; i < num; i++) val += ( Rand_())/((long double  )RAND_MAX) - 0.5;


  val = val * sig + a ;

  return val ;
}

// получение случайного числа с равномерным на отрезке [0;1]  распределением
long double   TShipTraj::getRand01( )
{

  return (Rand_())/ ((long double  )RAND_MAX);
}
//
long double   TShipTraj::Rand_()
{
  #ifndef DEBUG
  static int rand_init = 0;
  if( 0  == rand_init )
  {
	  rand_init = time( NULL );
	  srand((unsigned) time(NULL));
  }

  #endif

 // return (long double  )rand();
 // long double   breturn =  (long double  )rand(); // для отладки
  return (long double  )rand();                  // для отладки

}
 
#pragma package(smart_init)
