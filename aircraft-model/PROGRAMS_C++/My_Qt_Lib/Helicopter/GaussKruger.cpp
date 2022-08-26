#include <math.h>
#include "GaussKruger.h"

// в вычислениях единицы измерения:
// координаты - в метрах
// времена - в секундах
// скорости - в м/с
// углы - в радианах
// температура - в градусах Цельсия
// направление и курс отсчитывается от севера по часовой стрелке
// направление ветра показывает откуда дует ветер

#define  kfMhPi_    M_PI
#define  kfMh2Pi_   M_PI*2.0

#define GradToRad_           ( M_PI/180.0 )
#define speed_east(V,Q)_     ( V * cos(Q*GradToRad_) )
#define We(V,Q)_     ( V * sin(Q) )
#define speed_north_(V,Q)    ( V * sin(Q*GradToRad_) )
#define Wn_(V,Q)    ( V * cos(Q) )
#define GKx_(A,X,Y) (cos(A)*X-sin(A)*Y)
#define GKy_(A,X,Y) (cos(A)*Y+sin(A)*X)

// INPUT:
// L - долгота в рад
// B - широта в рад
void    mhGeo2Gauss_( double B, double L, double *X, double *Y )
{
    double  sin2B, sinB2, sinB4, sinB6, cosB;
    double  PX, PY, l, l2, Lgr;
    int     n;
    sin2B = sin(2.*B);
    cosB  = cos(B);
    sinB2 = sin(B)*sin(B);
    sinB4 = sinB2*sinB2;
    sinB6 = sinB4*sinB2;
    // долгота в градусах в диапозоне от 0 до 360
    L = mhMod2Pi_(L);
    Lgr = (L*180.)/kfMhPi_;
    n   = static_cast<int>((Lgr + 6.)/6.);
    l = (Lgr-(3+6*(n-1)))/( 180./kfMhPi_ );
    l2 = l * l;
    PY = l2*(278194.-830174.*sinB2 + 572434.*sinB4-16010.*sinB6
             +l2*(109500.-574700.*sinB2+863700.*sinB4-398600.*sinB6));
    PX = l2*(270806.-1523417.*sinB2 +1327645.*sinB4 - 21701.*sinB6
             +l2*(79690.-866190.*sinB2 + 1730360.*sinB4 -945460.*sinB6));
    *Y = 6367558.4968*B
            - sin2B*(16002.89 + 66.9607*sinB2 + 0.3515*sinB4
                     -l2*(1594561.25 + 5336.535*sinB2 + 26.79*sinB4+0.149*sinB6
                          +l2*(672483.4-811219.9*sinB2 + 5420.*sinB4 -10.6*sinB6+ PY)));

    *X = (5+10*n)*100000.
            +l*cosB*(6378245. + 21346.1415*sinB2 +107.159*sinB4 +0.5977*sinB6
                     +l2*(1070204.16 - 2136826.66*sinB2 + 17.98*sinB4 -11.99*sinB6 +PX));

    return;
}

// INPUT:
// X - на восток
// Y - на север
void    mhGauss2Geo_( double X, double Y, double *B, double *L )
{
    double  sin2b, sinb2, sinb4;
    double  sin2Bo, sinBo2, sinBo4, sinBo6, cosBo;
    double  PdB, Pl;
    double  b, zo, zo2, l, Bo, dB;
    int     n;
    n = static_cast<int>( X * 0.000001 );
    b = Y / 6367558.4968;
    sin2b = sin(2.0*b);
    sinb2 = sin(b) * sin(b);
    sinb4 = sinb2 * sinb2;
    Bo = b + sin2b*(0.00252588685 - 0.00001491860*sinb2 + 0.00000011904*sinb4);
    sin2Bo = sin(2.0*Bo);
    sinBo2 = sin(Bo) * sin(Bo);
    sinBo4 = sinBo2 * sinBo2;
    sinBo6 = sinBo4 * sinBo2;
    cosBo  = cos(Bo);
    zo = ( X - (10.0 * n + 5.0) * 100000.0 ) / ( 6378245.0 * cosBo );
    zo2 = zo * zo;
    PdB = 0.042858 - 0.025318*sinBo2 + 0.014346*sinBo4 - 0.001264*sinBo6
            -zo2*( 0.01672 - 0.0063*sinBo2 + 0.01188*sinBo4 - 0.00328*sinBo4 );
    dB = -zo2*sin2Bo*( 0.251684631 - 0.003369263*sinBo2 + 0.000011276*sinBo4
                       -zo2*( 0.10500614 - 0.04559916*sinBo2 + 0.00228901*sinBo4 - 0.00002987*sinBo6
                              -zo2*PdB ) );
    Pl = 0.01225 + 0.09477*sinBo2 + 0.03282*sinBo4 - 0.00034*sinBo6
            -zo2*( 0.0038 + 0.0524*sinBo2 + 0.0482*sinBo4 + 0.0032*sinBo6 );
    l =  zo*( 1.0 - 0.0033467108*sinBo2 - 0.0000056002*sinBo4 - 0.0000000187*sinBo6
              -zo2*( 0.16778975 + 0.16273586*sinBo2 - 0.00052490*sinBo4 - 0.00000846*sinBo6
                     -zo2*( 0.0420025 + 0.1487407*sinBo2 + 0.005942*sinBo4 - 0.000015*sinBo6
                            -zo2*Pl ) ) );
    *B = Bo + dB;
    *L = 6.0*(n - 0.5)/(180.0 / kfMhPi_) + l;
    return;
}

double  mhMod2Pi_( double a_fVal )
{
    double  fAngle = fmod( a_fVal, kfMh2Pi_ );
    if ( fAngle<0 )
        fAngle += kfMh2Pi_;
    return  fAngle;
}
