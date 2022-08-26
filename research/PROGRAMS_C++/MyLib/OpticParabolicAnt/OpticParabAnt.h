//---------------------------------------------------------------------------

#ifndef OpticParabAntH
#define OpticParabAntH
class TRectWaveGuide;
class TURPolyLine;
class TComp;
class TURPointXY;
class TOpticParabAnt
{
public:
	// фокальный параметр зеркала антенны (ml/2 = фокусное расстояние)
	double ml;
	// диаметр зеркала антенны
	double mParabDiam;
   // к-во волноводов
   int mNumWaveGuide;
  // расстояние между центрами волноводов
   double mWaveGuideDist;
   //высота излучателя
   double mWaveGuideHeight;
   // массив волноводов
   TRectWaveGuide *mpArrWaveGuide;

	__fastcall ~TOpticParabAnt() ;

  TOpticParabAnt();
  // парам констр
   TOpticParabAnt(const double l,const double ParabDiam,const int NumWaveGuide
  ,const double WaveGuideDist,const double WaveGuideHeight,const TRectWaveGuide WaveGuide);
 // парам констр
 TOpticParabAnt(const double l,const double ParabDiam,const int NumWaveGuide
  ,const double WaveGuideDist,const double WaveGuideHeight,const double b);
 // конструктор копирования
 TOpticParabAnt (const TOpticParabAnt &R) ;
 // оператор присваивания
 TOpticParabAnt operator=(TOpticParabAnt  R);

 bool createPictWithRay(wchar_t *wchFoldPict, const double VAlContrReflAng
  , const double VAlTargAng, const double VAlXf);

 bool  BuildRayWayPln(const double VAlContrReflAng
  , const double VAlTargAng, const double VAlYf, TURPolyLine &plnRez);

TURPolyLine calcFront0 (const double VAlContrReflAng , const double VAlTargAng);

int   calcRaysQuant (const double VAlContrReflAng , const double VAlTargAng
  , double valY0, double valY1,double valStep);

bool  createPhaseFrontPict(wchar_t *wchFoldPict, const double VAlContrReflAng
  , const double VAlTargAng , double valY0, double valY1);

double totalLengthWaveGuide();

void createRaysArray(const double VAlContrReflAng , const double VAlTargAng
  , double valY0, double valY1,double valStep, int quantRays, TURPolyLine *plnarr) ;

void  calcVeerDiagrams(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const int ILenDiagrArr, TComp *cmparrVeerDiagrs );

void  calcDiagrArrayFromDirect(const double  VAlLambda
		 ,const double VAlContrReflAng, const double VAlTargAngCur, TComp *pcmparrDiagrCur ) ;

void  calcDiagrArrayFromDirect__(const double  VAlLambda
,const double VAlContrReflAng, const double VAlTargAngCur, TComp *pcmparrDiagrCur ) ;


TURPolyLine createWaveGuideSegm(const int INumWaveGuide  );

TURPointXY createWaveGuideCentre(const int INumWaveGuide  )  ;

void  createVeerDiagrams(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const int ILenDiagrArr, TComp *cmparrVeerDiagrs );

void  createDiagrArrayFromDirect(const double  VAlLambda
		 ,const double VAlContrReflAng, const double VAlTargAngCur, TComp *pcmparrDiagrCur );

TComp calcW( double valYTemp, double VAlLambda, const double VAlContrReflAlf, const double VAlTargTetta)  ;

void  createVeerDiagramsNew(wchar_t *wchFoldPict, const double  VAlLambda,  const double VAlContrReflAng
	 , const double VAlAngDiagrStep, const double VAlMaxTargTetta); // , TComp *cmparrVeerDiagrs );

void createWeightFuncArr(const double  VAlLambda,const double VAlContrReflAng,const int LEnInpCurrentArr
, const TURPointXY PNtWaveGuideCentre, TComp *pcmparrQ ) ;

void createInpCurrentArray(const double  VAlLambda,double valTargEps, const int LEnInpCurrentArr
  , TComp *pcmparrFGr, TComp *pcmparrFGr1) ;

double calTay(const int LEnInpCurrentArr, const int INum);

TComp calcInpCurrentInPoint(const double  VAlLambda,double valTargEps, double valTay
   ,const double VAlSegmLength, const double STepIntegr );

void createWeightFuncArrNew(const double  VAlLambda,  const double VAlContrReflAng,const int LEnInpCurrentArr
, const TURPointXY PNtWaveGuideCentre, TComp *pcmparrQ );



};

TComp calcEMFieldPlaneWave_(const double  VAlLambda,double valA0
   ,const double VAlSegmLength, const double STepIntegr );

TComp calcEMFieldPlaneWave(const double  VAlLambda,double valTay, double valWaveX0,
   double valAlfaWave,const double VAlSegmLength, const double STepIntegr );

#endif
