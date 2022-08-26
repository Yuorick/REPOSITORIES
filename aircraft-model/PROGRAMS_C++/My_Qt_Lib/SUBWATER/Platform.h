#ifndef PlatformORM_H
#define PlatformORM_H
#include "HidroRLS.h"


class QHidroRLS;
class TEnvironment;

class QPlatform
{
public:
    QHidroRLS *mpHidroRLS;

    double marrPosParams[6];

    QPlatform();

    QPlatform (const  QPlatform &R);

    QPlatform  &operator=( const QPlatform  &R);

    QPlatform ( QHidroRLS *pHidroRLS, const double *arrPosParams);

    virtual int createInputDataReport(wchar_t*FileName, const bool bHeader);

    virtual int createInheritedInputDataReport(wchar_t*FileName, const bool bHeader);






};

#endif // Platform_H
