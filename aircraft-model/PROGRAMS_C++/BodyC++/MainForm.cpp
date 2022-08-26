﻿// ---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "MainForm.h"
#include "SimpleBody_3D.h"
// ---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
#include "Plane.h"
#include "CircleCylinder.h"
#include "TruncatedConeCircle.h"

#include "Complicated_Body.h"

#include "URPolygon.h"
#include "URPointXY.h"
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "MatrixProccess.h"
#include <stdlib.h>
#include "ArcParab.h"
#include "URPolyLine.h"
#include "Circle.h"
#include "YrRastr.h"
#include "URPolyLine.h"
#include "YrWriteShapeFile.h"

TForm1 *Form1;

// ---------------------------------------------------------------------------
__fastcall TForm1::TForm1(TComponent* Owner) : TForm(Owner) {

}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button1Click(TObject *Sender)
{
  // задание данных для апроксимации тела вертолета
  // в начальной сиситеме координат простейшими телами (для расчета центра масс и тензора инерции)
	const double VAlM = 8500;
	const double VAlH1 = 6.864;
	const double VAlR1 = 0.9813;
	const double VAlH2 = 3.076;
	const double VAlH3 = 2.017;
	const double VAlR3 = 1.145 / 2.;
	const double VAlH4 = 1.89;
	const double VAlX5 = 8.391;
	const double VAlR5 = 0.178 / 2.;
	const double VAlH5 = 2.324;
	///

	// создание апроксимации тела вертолета
	TComplicated_Body HeliBody(VAlM, VAlH1, VAlR1, VAlH2, VAlH3, VAlR3, VAlH4,
		VAlX5, VAlR5, VAlH5);
	///

	// вычисление координат центра масс в начальной сиситеме координат
	double arrCentreGrav[3] = {
		0.
	};
	HeliBody.calcCentreOfGravity(arrCentreGrav);
	///

	// формирование шейп файла с точкой центра масс в плоскости XY
	TURPointXY pntGravityCentreXY(arrCentreGrav[0] * 1000.,
		arrCentreGrav[1] * 1000.);
	pntGravityCentreXY.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\pntGravityCentreXY.shp", &pntGravityCentreXY, 1);
	///

	// формирование токи центра масс в 3 плоскостях начальной сиситемы координат
	TURPointXY pntGravCentreXY(-arrCentreGrav[0] * 1000., -arrCentreGrav[1]
		* 1000.);
	TURPointXY pntGravCentreXZ(-arrCentreGrav[0] * 1000., 0.);
	TURPointXY pntGravCentreYZ(0., -arrCentreGrav[1] * 1000.);
	 ///

	 // ценнтрирование тела вертолета относительно центра масс
	HeliBody.doCentreUp();
	///

	// вычисление матрицы моментов иннерции относительно центра масс
	double arrMtrxInertia[9] = {
		0.
	};
	HeliBody.calcInertiaMtrx(arrMtrxInertia);
	///


  // центрирование полигонов элементов планера вертолета относительно центра масс
  // вычисление центра масс каждого полигона в собственной связанной сиситеме координат
	// чтение shp файлов, сдвиг в ЦТ, нахождение центра масс полигона

	// 1. вид спереди  фюзеляж
	wchar_t wcharrTEmp[300] = {
		0
	};
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezFace.shp");
	int quantPlg = 1;

	TURPolygon *pPlg = (TURPolygon*)malloc(quantPlg*sizeof(TURPolygon));
	TURPolygon **ppPlg = &pPlg;
	// закачиваем файл с верт проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredFace = (*ppPlg)[0];
	plgCentredFace = plgCentredFace.SdvigTransform(pntGravCentreYZ);

	wchar_t wcharrTEmp1[300] = {
		0
	};
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredFace.shp");
	plgCentredFace.WriteSetSHPFiles(wcharrTEmp1, &plgCentredFace, 1);

	double cs = 5.;
	TURPointXY pntCentreFace = plgCentredFace.calcCentreMass(cs);  //   !!pntCentreFace
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreFace.shp");
	pntCentreFace.WriteSetSHPFiles(wcharrTEmp1, &pntCentreFace, 1);

	///

	// 2. вид сверху фюзеляж
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezFuz.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredFuz = (*ppPlg)[0];
	plgCentredFuz = plgCentredFuz.SdvigTransform(pntGravCentreXZ);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredFuz.shp");
	plgCentredFuz.WriteSetSHPFiles(wcharrTEmp1, &plgCentredFuz, 1);
	TURPointXY pntCentreFuz = plgCentredFuz.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreFuz.shp");
	pntCentreFuz.WriteSetSHPFiles(wcharrTEmp1, &pntCentreFuz, 1);

	///

	// 2. вид сверху  правый стаб
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezRightStab.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredRighttStab = (*ppPlg)[0];
	plgCentredRighttStab = plgCentredRighttStab.SdvigTransform(pntGravCentreXZ);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredRightStab.shp");
	plgCentredRighttStab.WriteSetSHPFiles
		(wcharrTEmp1, &plgCentredRighttStab, 1);
	TURPointXY pntCentreRighttStab = plgCentredRighttStab.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreRightStab.shp");
	pntCentreRighttStab.WriteSetSHPFiles(wcharrTEmp1, &pntCentreRighttStab, 1);

	///

	// 2. вид сверху  левый стаб
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezLeftStab.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredLefttStab = (*ppPlg)[0];
	plgCentredLefttStab = plgCentredLefttStab.SdvigTransform(pntGravCentreXZ);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredLeftStab.shp");
	plgCentredLefttStab.WriteSetSHPFiles(wcharrTEmp1, &plgCentredLefttStab, 1);
	TURPointXY pntCentreLefttStab = plgCentredLefttStab.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreLeftStab.shp");
	pntCentreLefttStab.WriteSetSHPFiles(wcharrTEmp1, &pntCentreLefttStab, 1);

	///

	// 2. вид сверху  правое крыло
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezRightWing.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredRighttWing = (*ppPlg)[0];
	plgCentredRighttWing = plgCentredRighttWing.SdvigTransform(pntGravCentreXZ);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredRightWing.shp");
	plgCentredRighttWing.WriteSetSHPFiles
		(wcharrTEmp1, &plgCentredRighttWing, 1);
	TURPointXY pntCentreRighttWing = plgCentredRighttWing.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreRightWing.shp");
	pntCentreRighttWing.WriteSetSHPFiles(wcharrTEmp1, &pntCentreRighttWing, 1);

	///

	// 2. вид сверху  левое крыло
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezLeftWing.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredLefttWing = (*ppPlg)[0];
	plgCentredLefttWing = plgCentredLefttWing.SdvigTransform(pntGravCentreXZ);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredLeftWing.shp");
	plgCentredLefttWing.WriteSetSHPFiles(wcharrTEmp1, &plgCentredLefttWing, 1);
	TURPointXY pntCentreLefttWing = plgCentredLefttWing.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreLeftWing.shp");
	pntCentreLefttWing.WriteSetSHPFiles(wcharrTEmp1, &pntCentreLefttWing, 1);

	///

	// 2. вид сверху фюзеляж
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezFuzVert.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredFuzVert = (*ppPlg)[0];
	plgCentredFuzVert = plgCentredFuzVert.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredFuzVert.shp");
	plgCentredFuzVert.WriteSetSHPFiles(wcharrTEmp1, &plgCentredFuzVert, 1);
	TURPointXY pntCentreFuzVert = plgCentredFuzVert.calcCentreMass(cs);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreFuzVert.shp");
	pntCentreFuzVert.WriteSetSHPFiles(wcharrTEmp1, &pntCentreFuzVert, 1);

	///
	// 2. вид сверху фюзеляж
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezRule.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredRule = (*ppPlg)[0];
	plgCentredRule = plgCentredRule.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredRule.shp");
	plgCentredRule.WriteSetSHPFiles(wcharrTEmp1, &plgCentredRule, 1);
	TURPointXY pntCentreRule = plgCentredRule.calcCentreMass(5);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreRule.shp");
	pntCentreRule.WriteSetSHPFiles(wcharrTEmp1, &pntCentreRule, 1);

	///

	// 2. вид сверху фюзеляж
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezShaft.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredShaft = (*ppPlg)[0];
	plgCentredShaft = plgCentredShaft.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredShaft.shp");
	plgCentredShaft.WriteSetSHPFiles(wcharrTEmp1, &plgCentredShaft, 1);
	TURPointXY pntCentreShaft = plgCentredShaft.calcCentreMass(5);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreShaft.shp");
	pntCentreShaft.WriteSetSHPFiles(wcharrTEmp1, &pntCentreShaft, 1);

	// 2. вид сверху фюзеляж
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plgRezSharnirs.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolygon::ReadSHPFile(wcharrTEmp, ppPlg, &quantPlg);
	TURPolygon plgCentredSharnirs = (*ppPlg)[0];
	plgCentredSharnirs = plgCentredSharnirs.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plgRezCentredSharnirs.shp");
	plgCentredSharnirs.WriteSetSHPFiles(wcharrTEmp1, &plgCentredSharnirs, 1);
	TURPointXY pntCentreSharnirs = plgCentredSharnirs.calcCentreMass(5);
	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\CEntreSharnirs.shp");
	pntCentreSharnirs.WriteSetSHPFiles(wcharrTEmp1, &pntCentreSharnirs, 1);

	///
	int quantPln = 1;

	TURPolyLine *pPln = (TURPolyLine*)malloc(quantPln*sizeof(TURPolyLine));
	TURPolyLine **ppPln = &pPln;
	// 2. вид сбоку линия стабилиатора
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plnRezStabLine.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolyLine::ReadSHPFile(wcharrTEmp, ppPln, &quantPln);
	TURPolyLine plnCentredStabLine = (*ppPln)[0];
	plnCentredStabLine = plnCentredStabLine.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plnCentredStabLine.shp");
	plnCentredStabLine.WriteSetSHPFiles(wcharrTEmp1, &plnCentredStabLine, 1);

	///

	// 2. вид сбоку линия крыла
	wcscpy(wcharrTEmp, mwchFoldInp);
	wcscat(wcharrTEmp, L"\\plnRezWingLine.shp");

	// закачиваем файл с передн  проекцией фюзеляжа
	TURPolyLine::ReadSHPFile(wcharrTEmp, ppPln, &quantPln);
	TURPolyLine plnCentredWingLine = (*ppPln)[0];
	plnCentredWingLine = plnCentredWingLine.SdvigTransform(pntGravCentreXY);

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\plnCentredWingLine.shp");
	plnCentredWingLine.WriteSetSHPFiles(wcharrTEmp1, &plnCentredWingLine, 1);

	///

	wcscpy(wcharrTEmp1, mwchFoldOut);
	wcscat(wcharrTEmp1, L"\\Axes.shp");
	TYrWriteShapeFile::CreateShpAxes(wcharrTEmp1, -20000., 20000., -20000.,
		20000.);
	free(pPlg);
}
// ------------------------------                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ---------------------------------------------

void __fastcall TForm1::Button5Click(TObject *Sender) {
	OpenDialog1->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog1->Execute()) {
		mpwchSHP_HorProj = (OpenDialog1->FileName).w_str();

	}
	Edit6->Text = mpwchSHP_HorProj;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button2Click(TObject *Sender) {
	int quantPlg = 1;

	TURPolygon *pPlg = (TURPolygon*)malloc(quantPlg*sizeof(TURPolygon));
	TURPolygon **ppPlg = &pPlg;
	TURPolygon::ReadSHPFile(mpwchSHP_HorProj, ppPlg, &quantPlg);
	TURPolygon plg0 = (*ppPlg)[0];
	int iii = 0;

	// линейное преобразование полигона
	// INPUT:
	// valAng - угол поворота
	// pntCentre - точка центра полигона
	// pntSdvig - точка куда перемещается центр полигона
	// valRastigenie - коэффициент растяжения
	// OUTPUT:
	// возвращает преобразованный полигон
	const TURPointXY pntSdvig(0., 0.);
	const double valRastigenieHor = 13870. / 2628.; // 7835./1463.;
	const double valAngHor = -atan
		((plg0.Points[plg0.NumPoints - 1].Y - plg0.Points[plg0.NumPoints - 2].Y) / (plg0.Points[plg0.NumPoints - 1]
			.X - plg0.Points[plg0.NumPoints - 2].X));
	TURPolygon plgRightFuz = plg0.LinTransform
		(valAngHor, plg0.Points[plg0.NumPoints - 2], pntSdvig, valRastigenieHor);
	///
	// plgRightFuz.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plg1.shp",&plgRightFuz, 1 )  ;

	TURPolygon plgLeftFuz = plgRightFuz.SimOtragenieTransform();
	// plgLeftFuz.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plg2.shp",&plgLeftFuz, 1 )  ;
	// plgRightFuz.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgRightFuz.shp",&plgRightFuz, 1 )  ;

	// перетаскивание крыла
	TURPolygon::ReadSHPFile(mpwchSHP_Wing, ppPlg, &quantPlg);
	TURPolygon plg3 = (*ppPlg)[0];
	TURPolygon plgRightWing = plg3.LinTransform
		(valAngHor, plg0.Points[plg0.NumPoints - 2], pntSdvig, valRastigenieHor);
	// plgRightWing.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgRightWing.shp",&plgRightWing, 1 )  ;

	TURPolygon plgLeftWing = plgRightWing.SimOtragenieTransform();
	// plgLeftWing.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgLeftWing.shp",&plgLeftWing, 1 )  ;

	// перетаскивание стабилизатора
	TURPolygon::ReadSHPFile(mpwchSHP_Stab, ppPlg, &quantPlg);
	plg3 = (*ppPlg)[0];
	TURPolygon plgRightStab = plg3.LinTransform
		(valAngHor, plg0.Points[plg0.NumPoints - 2], pntSdvig, valRastigenieHor);
	// plgRightStab.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgRightStab.shp",&plgRightStab, 1 )  ;

	TURPolygon plgLeftStab = plgRightStab.SimOtragenieTransform();
	// plgLeftStab.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgLeftStab.shp",&plgLeftStab, 1 )  ;

	// фыормирование левого фюзеляжа

	// симметр отражение относительно оси OY
	// симметричное отражение полигона относительно оси OY

	TURPolygon plgRezRightFuz = plgLeftFuz.SimOtragenieY_and_flip();
	TURPolygon plgRezLeftFuz = plgRightFuz.SimOtragenieY_and_flip();

	TURPolygon plgRezLeftWing = plgRightWing.SimOtragenieY_and_flip();
	TURPolygon plgRezRightWing = plgLeftWing.SimOtragenieY_and_flip();

	TURPolygon plgRezLeftStab = plgRightStab.SimOtragenieY_and_flip();
	TURPolygon plgRezRightStab = plgLeftStab.SimOtragenieY_and_flip();
	plgRezRightFuz.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezRightFuz.shp", &plgRezRightFuz, 1);
	plgRezLeftFuz.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezLeftFuz.shp", &plgRezLeftFuz, 1);

	plgRezLeftWing.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezLeftWing.shp", &plgRezLeftWing, 1);
	plgRezRightWing.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezRightWing.shp", &plgRezRightWing, 1);

	plgRezLeftStab.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezLeftStab.shp", &plgRezLeftStab, 1);
	plgRezRightStab.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezRightStab.shp", &plgRezRightStab, 1);

	free(pPlg);

}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button3Click(TObject *Sender) {
	OpenDialog2->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog2->Execute()) {
		mpwchSHP_Wing = (OpenDialog2->FileName).w_str();

	}
	Edit1->Text = mpwchSHP_Wing;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button4Click(TObject *Sender) {
	OpenDialog3->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog3->Execute()) {
		mpwchSHP_Stab = (OpenDialog3->FileName).w_str();

	}
	Edit2->Text = mpwchSHP_Stab;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button6Click(TObject *Sender) {
	OpenDialog4->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog4->Execute()) {
		mpwchSHP_VertProj = (OpenDialog4->FileName).w_str();

	}
	Edit3->Text = mpwchSHP_VertProj;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button8Click(TObject *Sender) {
	OpenDialog5->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog5->Execute()) {
		mpwchSHP_Rule = (OpenDialog5->FileName).w_str();

	}
	Edit4->Text = mpwchSHP_Rule;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button9Click(TObject *Sender) {
	OpenDialog6->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog6->Execute()) {
		mpwchSHP_Shaft = (OpenDialog6->FileName).w_str();

	}
	Edit5->Text = mpwchSHP_Shaft;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button10Click(TObject *Sender) {
	OpenDialog7->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog7->Execute()) {
		mpwchSHP_Sharnirs = (OpenDialog7->FileName).w_str();

	}
	Edit7->Text = mpwchSHP_Sharnirs;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button7Click(TObject *Sender) {
	int quantPlg = 1;

	TURPolygon *pPlg = (TURPolygon*)malloc(quantPlg*sizeof(TURPolygon));
	TURPolygon **ppPlg = &pPlg;
	// закачиваем файл с верт проекцией фюзеляжа
	TURPolygon::ReadSHPFile(mpwchSHP_VertProj, ppPlg, &quantPlg);
	TURPolygon plgFuzVert = (*ppPlg)[0];
	int iii = 0;
	///
	// закачиваем файл с рулем
	TURPolygon::ReadSHPFile(mpwchSHP_Rule, ppPlg, &quantPlg);
	TURPolygon plgRule = (*ppPlg)[0];
	//

	// закачивам фосьюайл с осью
	TURPolygon::ReadSHPFile(mpwchSHP_Shaft, ppPlg, &quantPlg);
	TURPolygon plgShaft = (*ppPlg)[0];
	///

	// закачивам фосьюайл с шарнирома
	int quantPlgShar = 1;
	TURPolygon::ReadSHPFile(mpwchSHP_Sharnirs, ppPlg, &quantPlgShar);
	TURPolygon plgSharnirs = TURPolygon(*ppPlg, quantPlgShar); // (*ppPlg)[0];
	///

	// закачиваем файл с осью СГФ

	TURPolyLine *pPln = (TURPolyLine*)malloc(quantPlg*sizeof(TURPolyLine));
	TURPolyLine **ppPln = &pPln;
	TURPolyLine::ReadSHPFile(mpwchSHP_SGF_Line, ppPln, &quantPlg);
	TURPolyLine plnSGF = (*ppPln)[0];

	// линейное преобразование полигона
	// INPUT:
	// valAng - угол поворота
	// pntCentre - точка центра полигона
	// pntSdvig - точка куда перемещается центр полигона
	// valRastigenie - коэффициент растяжения
	// OUTPUT:
	// возвращает преобразованный полигон
	const TURPointXY pntSdvig(0., 0.);
	const double valRastigenieHor = 13870. / 4134.;
	const double valAngHor = -atan((plnSGF.Points[1].Y - plnSGF.Points[0].Y) /
		(plnSGF.Points[1].X - plnSGF.Points[0].X));
	TURPolygon plgFuz1 = plgFuzVert.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolygon plgRezFuzVert = plgFuz1.SimOtragenieY_and_flip();
	plgRezFuzVert.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezFuzVert.shp", &plgRezFuzVert, 1);

	///

	TURPolygon plgRule1 = plgRule.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolygon plgRezRule = plgRule1.SimOtragenieY_and_flip();
	plgRezRule.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezRule.shp", &plgRezRule, 1);

	///

	///

	TURPolygon plgShaft1 = plgShaft.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolygon plgRezShaft = plgShaft1.SimOtragenieY_and_flip();
	plgRezShaft.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgRezShaft.shp", &plgRezShaft, 1);
	///

	TURPolygon plgSharnirs1 = plgSharnirs.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolygon plgRezSharnirs = plgSharnirs1.SimOtragenieY_and_flip();
	plgRezSharnirs.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\ plgRezSharnirs.shp", &plgRezSharnirs, plgRezSharnirs.NumParts);
	///

	// // закачивам фосьюайл с линией крыла

	TURPolyLine::ReadSHPFile(mpwchSHP_Stab_Line, ppPln, &quantPlg);
	TURPolyLine pln_Stab_Line = (*ppPln)[0];
	TURPolyLine pln_Stab_Line1 = pln_Stab_Line.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolyLine plnRezStabLine = pln_Stab_Line1.SimOtragenieY_and_flip();
	plnRezStabLine.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plnRezStabLine.shp", &plnRezStabLine, 1);

	///

	// закачивам фосьюайл с линией крыла

	TURPolyLine::ReadSHPFile(mpwchSHP_Wing_Line, ppPln, &quantPlg);
	TURPolyLine pln_Wing_Line = (*ppPln)[0];
	TURPolyLine pln_Wing_Line1 = pln_Wing_Line.LinTransform
		(valAngHor, plnSGF.Points[0], pntSdvig, valRastigenieHor);
	TURPolyLine plnRezWingLine = pln_Wing_Line1.SimOtragenieY_and_flip();
	plnRezWingLine.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plnRezWingLine.shp", &plnRezWingLine, 1);

	///

	//
	TYrWriteShapeFile::CreateShpAxes(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\Axes.shp", -20000., 20000., -20000., 20000.);

	free(pPln);
	free(pPlg);

}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button11Click(TObject *Sender) {
	OpenDialog8->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog8->Execute()) {
		mpwchSHP_SGF_Line = (OpenDialog8->FileName).w_str();

	}
	Edit8->Text = mpwchSHP_SGF_Line;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button12Click(TObject *Sender) {
	OpenDialog9->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog9->Execute()) {
		mpwchSHP_FaceProj = (OpenDialog9->FileName).w_str();

	}
	Edit9->Text = mpwchSHP_FaceProj;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button14Click(TObject *Sender) {
	OpenDialog10->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog10->Execute()) {
		mpwchSHP_Horiz_Line = (OpenDialog10->FileName).w_str();

	}
	Edit10->Text = mpwchSHP_Horiz_Line;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button13Click(TObject *Sender) {
	int quantPlg = 1;

	TURPolygon *pPlg = (TURPolygon*)malloc(quantPlg*sizeof(TURPolygon));
	TURPolygon **ppPlg = &pPlg;
	TURPolygon::ReadSHPFile(mpwchSHP_FaceProj, ppPlg, &quantPlg);
	TURPolygon plgFace = (*ppPlg)[0];
	///

	// закачиваем файл с осью

	TURPolyLine *pPln = (TURPolyLine*)malloc(quantPlg*sizeof(TURPolyLine));
	TURPolyLine **ppPln = &pPln;
	TURPolyLine::ReadSHPFile(mpwchSHP_Horiz_Line, ppPln, &quantPlg);
	TURPolyLine plnSGF = (*ppPln)[0];
	///

	// линейное преобразование полигона
	// INPUT:
	// valAng - угол поворота
	// pntCentre - точка центра полигона
	// pntSdvig - точка куда перемещается центр полигона
	// valRastigenie - коэффициент растяжения
	// OUTPUT:
	// возвращает преобразованный полигон
	const TURPointXY pntSdvig(0., 0.);
	const double valRastigenieHor = 2670. / 126.63; // 7835./1463.;
	const double valAngHor = -atan((plnSGF.Points[1].Y - plnSGF.Points[0].Y) /
		(plnSGF.Points[1].X - plnSGF.Points[0].X));
	TURPointXY pnt0(198.205, -207.191);
	TURPointXY pnt1(0., -1003.216);
	TURPointXY pnt2(pnt1.X - pnt0.X, pnt1.Y - pnt0.Y);
	// TURPolygon  plgFaceRight =plgFace.LinTransform(valAngHor , pnt2,valRastigenieHor ) ;

	TURPolygon plgFaceRight = plgFace.LinTransform(valAngHor, pnt0, pnt1,
		valRastigenieHor);
	///
	// plgFaceRight.WriteSetSHPFiles(L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\plgFaceRight.shp",&plgFaceRight, 1 )  ;

	TURPolygon plgFaceLeft = plgFaceRight.SimOtragenieY_and_flip();
	plgFaceRight.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgFaceRight.shp", &plgFaceRight, 1);
	plgFaceLeft.WriteSetSHPFiles(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\plgFaceLeft.shp", &plgFaceLeft, 1);

	TYrWriteShapeFile::CreateShpAxes(
		L"D:\\REPOSITORIES\\aircraft-model\\K-52_GRAPHS\\ScaledProject\\!_RezSobstvSK\\Axes.shp", -20000., 20000., -20000., 20000.);

	// TURPolygon plgRezLeftStab =  plgRightStab.SimOtragenieY_and_flip() ;

}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button15Click(TObject *Sender) {
	OpenDialog11->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog11->Execute()) {
		mpwchFileInp = (OpenDialog11->FileName).w_str();

	}
	wcscpy(mwchFoldInp, mpwchFileInp);
	wchar_t *pwchr = wcsrchr(mwchFoldInp, L'\\');
	pwchr[0] = 0;
	Edit11->Text = mwchFoldInp;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button16Click(TObject *Sender) {
	SaveDialog1->Filter = L"SHP файлы (*.shp)|*.shp";

	if (SaveDialog1->Execute()) {
		// ShowMessage( (SaveDialog1->FileName).w_str()) ;
		mpwchFileOut = (SaveDialog1->FileName).w_str();

	}
	wcscpy(mwchFoldOut, mpwchFileOut);
	wchar_t *pwchr = wcsrchr(mwchFoldOut, L'\\');
	pwchr[0] = 0;
	Edit12->Text = mwchFoldOut;

}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button17Click(TObject *Sender) {
	OpenDialog12->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog12->Execute()) {
		mpwchSHP_Wing_Line = (OpenDialog12->FileName).w_str();

	}
	Edit13->Text = mpwchSHP_Wing_Line;
}
// ---------------------------------------------------------------------------

void __fastcall TForm1::Button18Click(TObject *Sender) {
	OpenDialog13->Filter = L"файлы с графиками (*.shp)|*.shp";

	if (OpenDialog13->Execute()) {
		mpwchSHP_Stab_Line = (OpenDialog13->FileName).w_str();

	}
	Edit14->Text = mpwchSHP_Stab_Line;
} // fgets
// ---------------------------------------------------------------------------
