// ---------------------------------------------------------------------------

#ifndef PeaceVessH
#define PeaceVessH



#include "Environment.h"
#include "PeaceSins.h"
#include "Platform.h"
#include "Gps.h"
class QTrueMeasParams;



class QPeaceSins;
class QZamer;
class QAbstractMeasure;
class TTable_1D;

class QPlatform;
class QPeaceVess
{
public:
    QPeaceSins mSins;

    QPlatform mPlatform;

    QGps mGps;
    //первые 3 координаты - векьтор параллакса
    // далее - betta, eps, alf
	// параметры корабля
	double mWidth; // ширина(м)
	double mLength; // длина(м)

	double mMaxQ; // максимальный угол отклонения от генерального курса(амплитуда угла рыскания)
	double mT_Q; // период рыскания
	double mMaxPsi; // максимальный угол килевой качки(амплитуда)
	double mT_Psi; // период килевой качки
	double mMaxTet; // максимальный угол боротовой качки(амплитуда)
	double mT_Tet; // период бортовой качки
	double mMaxVert; // амплитуда вертикальной качки

	// парамеитры движения
	double mQ0; // генеральный курс
	double mVVess; // скорость корабля своего

	double mTVess; // время привязки траекторной информации

	// protected:
	// парметры углов палубы в центре качания
	double mQ; // угол курса
	double mPsi; // угол килевой качки
	double mTet; // угол бортовой качки
	double mVQ; // скорость изменения угла курса
	double mVPsi; // скорость изменения угла килевой качки
    double mVTet; // вкорость изменения угла бортовой качки

    // углы начальных фаз колебания угла курса, килевой качки, бортовой качки, высоты
	double marrDelt[4];
    // вектор состояния корабля в ГСК истинный
    double marrVectSost[9];

    // вектор состояния корабля в ГСК начальный
   // double marrVectSost0[9];


	// конструктор по умолчанию
    QPeaceVess();
	// конструктор копирования
    QPeaceVess(const QPeaceVess &R);
	// оператор присваивания
    QPeaceVess &operator=(const QPeaceVess  &R);

    // парам конструктор 1
    QPeaceVess (const TEnvironment Environment,const double Width,const double Length
                ,const  double MaxQ ,const  double T_Q
                ,const double MaxPsi,const double T_Psi ,const  double MaxTet
                ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
                , const QPeaceSins Sins
                ,const  QPlatform Platform, const double TVess, const QGps Gps);

     // парам конструктор 2
    QPeaceVess (const TEnvironment Environment,const double Width,const double Length
                     ,const  double MaxQ ,const  double T_Q
                     ,const double MaxPsi,const double T_Psi ,const  double MaxTet
                     ,const double T_Tet,const double MaxVert, const double Q0,const double VVess
                     , const QPeaceSins Sins,const  QPlatform Platform, const double TVess
                     , const QGps Gps, const double *arrS_Vess0);




    void tuneCurrentSinsMeasures();

    void extrapolateTrueVS_GSK(const double VAlTExtr,
        double *arrVessExtrapVS_GSK);

    void calcCentreDeckAngles(const double valT);

	void calcDeckAngles(const double valT, double *arrPointPositionPSK,
		double *pvalQ, double *pvalVQ, double *pvalPsi, double *pvalVPsi, double *pvalTet,
        double *pvalVTet);

    bool recalcVess(const double valT);

    void Move(const double valT);



    bool createTrueMeasureParams( TTable_1D &tblPrfl,const double TObrabotki
                                  , double *arrTrueBeaconPos,QTrueMeasParams &TrueMeasParams  );

    bool createTrueMeasureParams( TTable_1D &tblPrfl,const double TObrabotki
             , double *arrTrueBeaconPos, double *arrTrueBeaconVelo,QTrueMeasParams &trueMeasParams  );


     //static void  recalcZamerFrom_SphSK_to_SobSK_Gdg(QZamer InpSphPSKZamer, QZamer *pOutPSKZamer ); // +



    // void  integrate_arr_dEilers_po_dt(double *arr_dEilers_po_dt); //+



    static void  RecalcVect_KGSK_INTO_GdgPSK(double *arrKGSK  //+
            ,double *arrEilerCntrKP, double *arrVectOmegaPSK0
           , double *arrGdgPosParams,double *arrPSK,const int lenarrKGSK );

     static void recalcVect_KGSK_INTO_GdgSobSK( double *arrKGSK  //+
           , double *arrEilerCntrKP, double *arrVectOmegaPSK0
           ,  double *arrGdgPosParams,double *arrSobSK,const int lenarrKGSK );

     static void recalcPositionFromKGSK_to_GdgSphericalSK(double *arrEilerCntrKP //+
                 , double *arrKGSK,double *arrGdgPosParams,double *arrVGdg);

     static void recalcPositionFromGSK_to_GdgSphericalSK(double *arrSVess_gsk,double *arrEilerCntrKP
                , double *arrSObj_gsk,double *arrGdgPosParams,double *arrVGdg);

     static void recalcPositionFromGdgSphericalSK_to_KGSK(double *arrVGdg,double *arrEilerCntrKP//+
                               ,double *arrGdgPosParams, double *arrKGSK);

     //static void recalcZamerFrom_GdgSphericalSK_to_KGSK(QZamer *pZamerGdgSh,double *arrEilerCntrKP //+
         //      ,double *arrGdgPosParams, QZamer**ppZamerKGSK);

    // static void recalcZamerFrom_GdgSphericalSK_to_PSK(QZamer *pZamerGdgSh,double *arrEilerCntrKP
            //                                                     ,double *arrGdgPosParams, QZamer **ppZamerPSK);



     void integrateEilersVect (double *arrMu);

     void  integrate_arr_dEilers_po_dt(double *arr_dMu);



     static void RecalcVect_KGSK_INTO_GdgPSK_differentiation( double *arrKGSK
                                           , double *arrEilerCntrKP, double *arrVectOmegaPSK0,  double *arrGdgPosParams
                                            ,double *arrPSK,const int lenarrKGSK );


     int createInputDataReport(wchar_t*FileName, const bool bHeader);

     void createTrajectReport(const double VAlBeginT,const double VAlEndT
                              ,const double VAlStepT,double* parrData , int *piQuantRows );

     static void  RecalcVect_PSK_INTO_GdgPSK( double *arrPSK_CT, double *arrGdgPosParams
                                       ,double *arrGdgPSK,const int lenarrKGSK );

     static void recalcVect_PSK_INTO_GdgSobSK( double *arrPSK_CT
                ,  double *arrGdgPosParams,double *arrSobSK,const int lenarrKGSK );



     static void  RecalcVect_KGSK_INTO_PSK_CT( double *arrKGSK
         , double *arrEilers0, double *arr_dEilers0_po_dt
          ,double *arrPSK_CT,const int lenarrKGSK );

     static void  RecalcVect_PSK_CT_INTO_KGSK( double *arr_PSK_CT, double *arrEilers0
           , double *arr_dEilers0_po_dt,double *arr_KGSK,const int lenarrKGSK );




     bool calcOtvetT( TTable_1D &tblPrfl,const double TObrabotki,double *arrTrueBeaconPos, double &valOtvetT);

     bool calc_tBeam_and_dtBeam_po_dT( TTable_1D &tblPrfl,const double *arrTrueBeaconPos
                                      , double &val_tBeam, double &val_dtBeam_po_dT);

     void  RecalcVect_PSK_CT_INTO_GSK( double *arr_PSK_CT, double *arrEilers0
                 , double *arr_dEilers0_po_dt,double *arr_GSK,const int lenarrKGSK );


     void calc_dYAnt_po_dT(double *arr_dYAnt_po_dT);

     static double calcCourseAngle(const double x, const double y);
     

     
     bool calcOtvetT_( TTable_1D &tblPrfl,const double TObrabotki
                       ,double *arrBeaconPos,double *arrBeaconVelo, double &valOtvetT
                       , double &val_tOtv, double &val_tZapr, double *arrBeaconPosInp
                       , double *arrBeaconPosOut);

};
#endif
