#ifndef BALLANCECALC_H
#define BALLANCECALC_H

class THelic;
class TEnvironment;
class TBallanceCalc
{
public:
    TBallanceCalc();

static bool calcBallParamsForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
                                              ,const long double VAlY,const long double VAlVx,const long double VAlPsi
                                              ,long double *arrX0, long double *arrXRez );

static void calcFgr_and_JacFgr_ForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
        , long double VAlY, long double VAlVx, long double VAlPsi
        ,long double *arrXCur, long double *arrFGr, long double *arrJacFgr);

static void calcFgr_ForSteadyLineMoving(const THelic Helic,const TEnvironment Environment
                                                       , long double VAlY, long double VAlVx, long double VAlPsi
                                                       ,long double *arrXCur, long double *arrFGr);

static bool calcBallParamsForTurn(const THelic Helic,const TEnvironment Environment
                                                 ,const long double VAlY,const long double VAlVx,const long double VAlPsi
                                                 ,const long double VAlRad,long double *arrX0, long double *arrXRez );

static void calcFgr_and_JacFgr_ForTurn(const THelic Helic,const TEnvironment Environment
                                                      , long double VAlY, long double VAlVx, long double VAlPsi, const long double VAlRad
                                                      ,long double *arrX, long double *arrFGr, long double *arrJacFgr);

static void calcFgr_ForTurn(const THelic Helic,const TEnvironment Environment
                                           , long double VAlY, long double VAlVx, long double VAlPsi, const long double VAlRad
                                           ,long double *arrX, long double *arrFGr);

static bool calcBallParamsForRotation(const THelic Helic,const TEnvironment Environment
                                                     ,const long double VAlY,const long double VAl_dPsi_po_dt
                                                     ,long double *arrX0, long double *arrXRez );

static void calcFgr_and_JacFgr_ForRotation(const THelic Helic,const TEnvironment Environment
                                                          , long double VAlY, const long double VAl_dPsi_po_dt
                                                          ,long double *arrX, long double *arrFGr, long double *arrJacFgr);

static void calcFgr_ForRotation(const THelic Helic,const TEnvironment Environment
                                    , long double VAlY,const long double VAl_dPsi_po_dt
                                    ,long double *arrXCur, long double *arrFGr);

static long double calcAlfaZakl(const long double valNu,const long double valXT,const long double valDist);
};

#endif // BALLANCECALC_H
