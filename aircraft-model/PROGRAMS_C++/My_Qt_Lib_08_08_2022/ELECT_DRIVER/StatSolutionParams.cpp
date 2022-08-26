#include "StatSolutionParams.h"
#include <string.h>

QStatSolutionParams::QStatSolutionParams()
{
    memset(marrStatPhVect,0, QVARS * sizeof(double));
    memset(marrStatU,0, 2 * sizeof(double));
    memset(marrGears,0, 2 * QVARS * sizeof(double));

}
// конструктор копирования
 QStatSolutionParams ::QStatSolutionParams (const QStatSolutionParams &R)
 {
     memcpy(marrStatPhVect, R.marrStatPhVect, QVARS * sizeof(double));
     memcpy(marrStatU, R.marrStatU, 2 * sizeof(double));
     memcpy(marrGears, R.marrGears, 2 * QVARS * sizeof(double));
 }
 // оператор присваивания
 QStatSolutionParams &QStatSolutionParams::operator=(const QStatSolutionParams  &R)
 {
     memcpy(marrStatPhVect, R.marrStatPhVect, QVARS * sizeof(double));
     memcpy(marrStatU, R.marrStatU, 2 * sizeof(double));
     memcpy(marrGears, R.marrGears, 2 * QVARS * sizeof(double));
   return *this ;
 }
// парам констр
QStatSolutionParams :: QStatSolutionParams(  const double *arrStatPhVect,const double* arrStatU,const double*arrGears)

 {
    memcpy(marrStatPhVect, arrStatPhVect, QVARS * sizeof(double));
    memcpy(marrStatU, arrStatU, 2 * sizeof(double));
    memcpy(marrGears, arrGears, 2 * QVARS * sizeof(double));
}
