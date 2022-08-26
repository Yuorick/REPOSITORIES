#ifndef ACDFILE_H
#define ACDFILE_H

#include <QString>
#include "BigMeasure.h"

/*bool writeExchangeDataFile(QString wchReportFile,
              const int QuantMeas// к-во измерений
              ,QBigMeasure* parrMeas//
              , double* arrAntPosParams
              ,double* arrGPSPosParams
              ,double* arrBeaconPos
              ,const double TObrabotki
              );*/

bool writeExchangeDataFile_(QString wchReportFile,
              const int QuantMeas// к-во измерений
              ,QBigMeasure* parrMeas//
              , double* arrAntPosParams
              ,double* arrGPSPosParams
              ,double* arrBeaconPos
              ,const double TObrabotki
              );
void readFile_acd(QString TargFileName,double *arrP_Ant
                  , double *arrP_Gps, double *arrS_Beacon
                  , double &valTBeacon, QVector<QBigMeasure> *pvctBigMesure);

void readQuantMeasures_From_File_acd(QString DataFileName,int &iQuantMeasures);

void readBigMeasuresOnly_acd(QString TargFileName,const double *arrP_Ant
                             , double *arrP_Gps, double *arrS_Beacon,const  double valTBeacon
                             , QVector<QBigMeasure> *pvctBigMesure);

#endif // ACDFILE_H
