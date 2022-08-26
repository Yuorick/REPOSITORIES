//---------------------------------------------------------------------------


#pragma hdrstop
 #include <string.h>
 #include <stdlib.h>
#include "AntFar.h"
#include "URPolygon.h"
#include "EmitterImage.h"
#include "YrWriteShapeFile.h"
// INPUT:
// wchFoldName - путь к папке с файлами
// Emitter:
// ValLength0 - длина горизонтальной части эмиттера
// ValLength1 - длина наклонной части эмиттера
// ValFi - угол между горизонтальной и наклонной частями эмиттера
//  AM:
// QuantEmit - к-во эмиттеров
// DistEmit - расстояние мжде эмиттерами
// ValLengthBoardAM - расстояние между границами АМ
// SideWidthAM - толщина боковых частей полигона кожуха эмиттеров
// BackWidthAM - толщина задней части  полигона кожуха эмиттеров
// LenghPlnOutAM  - длина линии выхода с АМ
// FAR:
// QuantAM - к-во АМ
//
//

void createFarPict(wchar_t *wchFoldName, const int QuantEmit,const double DistEmit
  , const int QuantAM, const double ValLengthBoardAM
  ,const double SideWidthAM, const double BackWidthAM , const double LenghPlnOutAM
  ,const double ValLength0 ,const double ValLength1,const double ValFi)

{    // рассточние между центьрами АМ
	  double valDistAM =  ((double)QuantEmit)*DistEmit  + ValLengthBoardAM * 2.;
	// mсуммарная высота ФАР
	double valHeight =  (((double)QuantEmit)*DistEmit  + ValLengthBoardAM* 2.) * ((double)QuantAM);
	double valWidth = ValLength0 + ValLength1 * cos(ValFi) * 3. + BackWidthAM;// ширина
	double valRad = valWidth/ 2.;  // радиус закругления
	// оьолочка ФАР
	TURPolygon plgFarHood1 = TURPolygon::createFarHood( valHeight, valWidth, valRad);
	TURPointXY pntHoodSdvig(valWidth/2. - (ValLength0 + ValLength1 * cos(ValFi)  + BackWidthAM) ,0.);
	TURPolygon plgFarHood = plgFarHood1.SdvigTransform(pntHoodSdvig );
	///
	// построение  АМ

		 // эмиттер
	TURPolyLine plnarrEmit[3];
	TURPointXY pnt[4];
	pnt[0] = TURPointXY( -ValLength0 - ValLength1* cos(ValFi) ,0.);
	pnt[1] = TURPointXY( - ValLength1* cos(ValFi),0.);
	pnt[2] = TURPointXY( 0.,ValLength1* sin(ValFi));
	pnt[3] = TURPointXY( 0., -ValLength1* sin(ValFi));
	plnarrEmit[0] =  TURPolyLine(pnt[0], pnt[1]);
	plnarrEmit[1] =  TURPolyLine(pnt[1], pnt[2]);
	plnarrEmit[2] =  TURPolyLine(pnt[1], pnt[3]);

	 // центр эмиттера
	 TURPointXY pntEmitCentre(0.,0.);
	   // кожух АМ
	TURPolygon plgAMHood = TURPolygon ( 9) ;

	double valLenEmit = plnarrEmit[0].Points[1].X- plnarrEmit[0].Points[0].X ;

	double valLen = ((double)(QuantEmit -1)) *  DistEmit/2.; // высота АМ

	plgAMHood.Points[0] =  TURPointXY (plnarrEmit[0].Points[1].X- 0.25 * valLenEmit,valLen +DistEmit/2. );
	plgAMHood.Points[8] = plgAMHood.Points[0] ;
	plgAMHood.Points[3] = TURPointXY( plgAMHood.Points[0].X ,- plgAMHood.Points[0].Y);
	plgAMHood.Points[1] =  TURPointXY (plnarrEmit[0].Points[0].X, plgAMHood.Points[0].Y);
	plgAMHood.Points[2] = TURPointXY( plgAMHood.Points[1].X, -plgAMHood.Points[1].Y);
	plgAMHood.Points[7] = TURPointXY (plgAMHood.Points[0].X, plgAMHood.Points[0].Y + SideWidthAM);
	plgAMHood.Points[4] = TURPointXY (plgAMHood.Points[7].X,-plgAMHood.Points[7].Y);
	plgAMHood.Points[6] = TURPointXY (plgAMHood.Points[1].X -BackWidthAM   ,plgAMHood.Points[7].Y);
	plgAMHood.Points[5] = TURPointXY (plgAMHood.Points[6].X,-plgAMHood.Points[6].Y);

   // выход из АМ
	TURPointXY  pnt1(plgAMHood.Points[5].X, 0.);
	TURPointXY  pnt2 (pnt1.X - LenghPlnOutAM, 0.);
	// выход АМ
   //	TURPolyLine plnOutAM(   pnt1,  pnt2) ;
	TURPolyLine plnOutAM = TURPolyLine::fncCreateArrow( pnt1,  pnt2, LenghPlnOutAM/10., 20./ 180. * M_PI);
	///
	// ИЗЛУЧАТЕЛИ в одном модуле   с центрами
	TURPointXY *pPntEmitCentresAM = (TURPointXY *) malloc( 3 * QuantEmit * sizeof(TURPointXY));
	TURPolyLine *pplnEmitAM = (TURPolyLine *) malloc( 3 * QuantEmit * sizeof(TURPolyLine));
	for (int i =0; i < QuantEmit; i++)
	{
	 TURPointXY pntSdvig (0., valLen -  ((double)i) *DistEmit);
	// TURPointXY  pnt00 = Emit.marrPlnSegm.    (, pnt01, pnt02, pnt03;
	 pplnEmitAM [3 * i]    =  plnarrEmit[0].SdvigTransform(pntSdvig ) ;
	 pplnEmitAM [3 * i +1] =  plnarrEmit[1].SdvigTransform(pntSdvig ) ;
	 pplnEmitAM [3 * i +2] =  plnarrEmit[2].SdvigTransform(pntSdvig ) ;
	 pPntEmitCentresAM[i].X =  pntEmitCentre.X + pntSdvig.X;
	 pPntEmitCentresAM[i].Y =  pntEmitCentre.Y + pntSdvig.Y;
	}
	  // центр АМ
	TURPointXY pntCentreАМ (0.,0.);
	///

	// формование ФАР
		// подготовкам массивов
	  // массивы ФАР
	   // массив излучателей
		  TURPolyLine *pPlnArrEmitFAR = (TURPolyLine *) malloc( 3 * QuantEmit * QuantAM * sizeof(TURPolyLine));
	   // массив точек центорв излучателей
		  TURPointXY *pPntArrEmitCentres = (TURPointXY *) malloc( 3 * QuantEmit * QuantAM * sizeof(TURPointXY));
	   ///
	   // массив центров АМ
		  TURPointXY *pPntArrАМCentres = (TURPointXY *) malloc( 3 *  QuantAM * sizeof(TURPointXY));
	   // массив кожухов АМ
		 TURPolygon *pPlgArrHoodsAM = (TURPolygon *) malloc( 3 *  QuantAM * sizeof(TURPolygon));
	   // массив выходов АМ
		TURPolyLine *pPlnArrOutLineAM = (TURPolyLine *) malloc( 3 *  QuantAM * sizeof(TURPolyLine));
		///

		//    заполнение массивов
		for (int i = 0; i < QuantAM; i++)
		{
		   TURPointXY pntSdvig(0.,valDistAM * ((double)QuantAM -1)/2. - ( ((double)i) * valDistAM ));
		   pPntArrАМCentres[i] = pntCentreАМ.SdvigTransform(pntSdvig ) ;
		   pPlgArrHoodsAM[i] =  plgAMHood.SdvigTransform(pntSdvig ) ;
		   pPlnArrOutLineAM [i]=  plnOutAM.SdvigTransform(pntSdvig ) ;
		   for (int j = 0; j < QuantEmit * 3; j++)
		   {
			 pPlnArrEmitFAR[i * QuantEmit * 3+ j] = pplnEmitAM[j].SdvigTransform(pntSdvig ) ;
		   }

		   for (int j = 0; j < QuantEmit ; j++)
		   {
			 pPntArrEmitCentres[i * QuantEmit + j] = pPntEmitCentresAM[j].SdvigTransform(pntSdvig ) ;
		   }

		}




	wchar_t wchFileName[300]= {0};
	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\HoodFAR.shp");
	plgFarHood.WriteSetSHPFiles(wchFileName,&plgFarHood, 1) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\EmittersAM.shp");
	TURPolyLine::WriteSetSHPFiles(wchFileName,pplnEmitAM , 3 * QuantEmit) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\AMCentres.shp");
	TURPointXY::WriteSetSHPFiles(wchFileName,pPntArrАМCentres , QuantAM) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\AMHoods.shp");
	TURPolygon::WriteSetSHPFiles(wchFileName,pPlgArrHoodsAM , QuantAM) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\AMOutLines.shp");
	TURPolyLine::WriteSetSHPFiles(wchFileName,pPlnArrOutLineAM , QuantAM) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\FAREmitters.shp");
	TURPolyLine::WriteSetSHPFiles(wchFileName,pPlnArrEmitFAR , 3 * QuantEmit * QuantAM) ;

	wcscpy( wchFileName,  wchFoldName);
	wcscat(wchFileName, L"\\FAREmitCentres.shp");
	TURPointXY::WriteSetSHPFiles(wchFileName,pPntArrEmitCentres ,  QuantEmit * QuantAM) ;

   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ///////   создание файлов, связанных с волной /////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// угол цели
  double valAl = 10./ 180. * M_PI;
 // линия фазового фронта проходящая через фазовый центр первой строки
 TURPointXY pntPhaseCntr = pPntArrАМCentres[0];
 double valxfront = (valHeight + 10.)/cos(valAl);
 TURPointXY pntFront0( pntPhaseCntr.X - valxfront*sin(valAl), pntPhaseCntr.Y + valxfront*cos(valAl) );
 TURPointXY pntFront1 ( pntPhaseCntr.X + valxfront*sin(valAl), pntPhaseCntr.Y - valxfront*cos(valAl) );
 TURPolyLine plnFront(pntFront0,pntFront1);
 // точки пересечения лучей с линией фазового фронта
 TURPointXY *arrPointFront =  new TURPointXY[1 + QuantAM];
 for (int i = 0 ; i < (QuantAM ); i++)
 {
  double vald = pntPhaseCntr.Y - pPntArrАМCentres[i].Y;
  arrPointFront[i] = TURPointXY(vald * cos(valAl) * sin(valAl)
	  ,pntPhaseCntr.Y - vald * cos(valAl) * cos(valAl));
 }
 arrPointFront[QuantAM ] = pntPhaseCntr;

 // линии с сигналом от цели проходят через точки arrPoints
 TURPolyLine *arrPlyline   = new TURPolyLine[1 + QuantAM];
 TURPolyLine *arrSegmSdvig = new TURPolyLine[1 + QuantAM] ;

 double valx =100.;
 for (int i = 0; i < (1 + QuantAM); i++)
 {
   TURPointXY  pnt2(valx,  valx* tan(valAl) + pPntArrАМCentres[i].Y);
   arrPlyline[i] = TURPolyLine(arrPointFront[i],   pnt2) ;
   arrSegmSdvig[i] = TURPolyLine(pPntArrАМCentres[i], arrPointFront[i]);
 }

 // нормаль к ФАР
 TURPointXY pntNorm(100.,pntPhaseCntr.Y);
 TURPointXY pnt00(0.,pntPhaseCntr.Y);
 TURPolyLine plnNorm(pntNorm,pnt00) ;
///

// линия фазовых центров излучателей первого модуля
  TURPointXY pntTemp(0.,pPntArrEmitCentres[0].Y + 20.);
  TURPolyLine plnAM0(pntTemp,pPntArrEmitCentres[QuantEmit-1]) ;

   // файлы


 wchar_t   pwcharrFile[300] ={0};
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Arrays.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,arrPlyline, QuantAM) ;

	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\plnNOrm.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,&plnNorm, 1) ;



  // обозначение угла между нормалью и лучом
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect0.shp");
 // TURPointXY  pnt0(0.,0.);
  double Dist0 = 12., Dist1 = 13.;
  TYrWriteShapeFile::CreateAngleMarks(pwcharrFile, pnt00
	 , pntNorm, arrPlyline[0].Points[1],Dist0, Dist1 ) ;
	 ///

 // обозначение угла между фронтом и ФАР
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect1.shp");


  TYrWriteShapeFile::CreateAngleMarks(pwcharrFile, pnt00
	 , pPntArrEmitCentres[0], pntFront0,Dist0+ 6, Dist1+ 6 ) ;
	 ///

	 // обозначение угла между фронтом и ФАР
  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Sect2.shp");


  TYrWriteShapeFile::CreateAngleMarks(pwcharrFile, pnt00
	 , pPntArrEmitCentres[QuantEmit-1], pntFront1,Dist0, Dist1 ) ;
	 ///


  wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\plnAM0.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,& plnAM0, 1) ;

  // линия фронта волны
	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\plnFront.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,&plnFront, 1) ;
  // отрезки фазоовых сдвигов

	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\sgmSdvig.shp");
  TURPolyLine::WriteSetSHPFiles(pwcharrFile,arrSegmSdvig, QuantAM) ;
  // фазовый центр
  	wcscpy(pwcharrFile, wchFoldName);
  wcscat(pwcharrFile, L"\\Centre.shp");
  TURPointXY::WriteSetSHPFiles(pwcharrFile,&pntPhaseCntr, 1) ;
  ///


  // Создание картинок строковых диаграмм
   	// Создание диаграммы с центром в точке  pointCentre
	// радиусом  valR с шагом по углу iNUM
	// valWidthMilrad -  ширина диаграммы в милирадианах
	// valRotAng - угол поворота относительно горизонтали
	// положит направление против час стрелки
  /*	double valWidth = 4./180. * M_PI ;
	const double valR = 12000.;
	const double valRotAng = 0.;//-2./180. * M_PI ;
	const int iNUM = 3000;
	TURPolygon Diagr = TURPolygon::fncCreateDiagr(pntCentreNew,  valR,valRotAng , valWidth, iNUM) ;
	wcscpy(pwcharrFileVess, pwcharrPath);
	wcscat(pwcharrFileVess, L"\\Diagr.shp");
	Diagr.WriteSetSHPFiles(pwcharrFileVess, &Diagr ,1); */

  ////
	delete [] arrPointFront;
	delete []arrPlyline ;
	delete []arrSegmSdvig;


	free (pplnEmitAM);
	free (pPlnArrEmitFAR );
	free(pPntArrEmitCentres);
	free(pPntArrАМCentres);
	free(pPlgArrHoodsAM);
	free(pPlnArrOutLineAM);
	free(pPntEmitCentresAM);

}

#pragma package(smart_init)
