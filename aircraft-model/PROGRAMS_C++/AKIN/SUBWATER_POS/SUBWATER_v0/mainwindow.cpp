#include "mainwindow.h"
#include "ui_mainwindow.h"

#include "SubWaterBeam.h"
#include "math.h"
#include "Table_1D.h"
#include "URPolyLine.h"
#include "UrPointXY.h"
#include <QFileDialog>
#include "YrWriteShapeFile.h"

extern bool BEZ_SHUMOV = true;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->ui->lineEdit_3->setText("D:\\AKIN\\PROFILES\\02102021//FILE48.000");
    this->ui->lineEdit->setText("D:\\AKIN\\ОТЧЕТ_ПОЗИЦИОНИРОВАНИЕ");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    // 0.
    QString strFold = this->ui->lineEdit_3->text();
    wchar_t array [400] = {0};
    strFold.toWCharArray(array);
    array[strFold.length()] = 0;
    wcscpy(mpwchPrflFIle, array);

    strFold = this->ui->lineEdit->text();
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;
    // 1. скачивание файла профиля
  int i =   QSubWaterBeam::calcColCountFrom000(mpwchPrflFIle);

  int iNumRows = 0;
  double parrData[1000] = {0.}, parrData1[1000] = {0.};
  QSubWaterBeam::ReadDataFrom000(mpwchPrflFIle, &iNumRows, parrData);
   //1!

  // 2.сортировка по высоте
  QSubWaterBeam::sortProfile(parrData, iNumRows);
  //2!

  //3. надо добавить самую глубокую точку
  parrData[iNumRows * 2]= parrData[(iNumRows-1) * 2] +5.;
  parrData[iNumRows * 2 + 1]= parrData[(iNumRows-1) * 2 + 1] - 0.01;
  ++iNumRows;

  int iNumRowsNew = -1;
  memcpy(parrData1, parrData,1000 * sizeof(double) );
  QSubWaterBeam::adjustProfile(parrData1,  iNumRows, iNumRowsNew);
 // 3!

  // 4.создание полилинии профиля
  TURPointXY pPoints[1000];

  for(int i = 0;i< iNumRowsNew; ++i)
  {

     pPoints[i].X = parrData1[2 * i +1]-1300.;
     pPoints[i].Y = -parrData1[2 * i];
  }
  TURPolyLine plnPrfl( pPoints,iNumRowsNew) ;
  // 4!

  //TTable_1D tblPrfl_Adjust(parrData1,iNumRowsNew);
  //создание профиля скорости звука
  TTable_1D tblPrfl_Adjust;
  QSubWaterBeam::createProfileTbl(mpwchPrflFIle, &tblPrfl_Adjust, VAR0);
 // 3!

  double valTetta = -1.;
  double VAlza= 8.,VAlzm =parrData1 [(iNumRowsNew -2) * 2],VAlxm = 100.;
  bool br=calcTetta( VAlza,VAlzm,VAlxm
                 ,tblPrfl_Adjust, valTetta);
  const double VAlCosTetta = cos(valTetta);
  double valTetCrit = calcTettaCrit(  VAlza, VAlzm,tblPrfl_Adjust) + 0.017;
  //cos(tetta +0.017)
  double tetStep = 10./180.*M_PI;
  int iNumParts = M_PI/2./tetStep;

  double step = 1.;
  int nc = (VAlzm - VAlza)/step ;

  const int iNumPoints = iNumParts * nc;
  double *arr = new double[iNumPoints *2];

  int *iarrParts = new int [iNumParts];
  memset(iarrParts, 0, iNumParts * sizeof(int));

  TURPolyLine pln1(  iNumParts,iNumPoints,iarrParts
                     ,arr);

  TURPolyLine pln1_Adjust(  iNumParts,iNumPoints,iarrParts
                     ,arr);


  TURPolyLine plnCurveLength(  1,iNumParts,iarrParts
                     ,arr);
  plnCurveLength.Parts[0] = 0;
  TURPolyLine plnDist(  1,iNumParts,iarrParts
                     ,arr);
  plnDist.Parts[0] = 0;

  int iCur = 0;
  for(int i = 0; i < iNumParts; ++i)
  {
      pln1.Parts[i]= i * nc;
      pln1_Adjust.Parts[i]= i * nc;
  }
  for(int i = 0; i < iNumParts; ++i)
  {
   double tetta0 = valTetCrit + ((double)i)*tetStep;
   double x = 0., z =0.;
  for(int j = 0; j < nc; ++j)
  {

     z = VAlza  +((double)j) *  step;

     x = calcXHoriz(VAlza, z
                          ,tblPrfl_Adjust,cos(tetta0));
    pln1.Points[iCur].Y = -z;
    pln1.Points[iCur].X = x ;

    ++iCur;
  }
  plnCurveLength.Points[i].X = tetta0* 180. / M_PI *10.;
  plnCurveLength.Points[i].Y = calc_CurveLength( VAlza,z,tetta0,tblPrfl_Adjust);

 // plnCurveLength.Points[i].Y = pln1.calcPartLeng(i);
  double lenCur = pln1.calcPartLeng(i);

  plnDist.Points[i].X = tetta0* 180. / M_PI *10.;
  plnDist.Points[i].Y = sqrt((z - VAlza) * (z - VAlza) + x * x);
  int yy=0;
  }


  delete []arr ;
  delete []iarrParts;


//TURPolyLine pln(pPoints,qpoints) ;
  wchar_t wchAxesFileName0[300] ={0};
  wcscpy(  wchAxesFileName0,  mwchOutPutFold);
  wcscat(wchAxesFileName0, L"\\AxesArr.shp");
  TYrWriteShapeFile::CreateShpArrowedAxes(wchAxesFileName0,-5000., 5000.
  ,-10000., 10000.,30.) ;

  wcscpy(  wchAxesFileName0,  mwchOutPutFold);
  wcscat(wchAxesFileName0, L"\\plnBeamsTraj.shp");
  pln1.WriteSetSHPFiles(wchAxesFileName0,&pln1,1);

  wcscpy(  wchAxesFileName0,  mwchOutPutFold);
  wcscat(wchAxesFileName0, L"\\plnCurveLength1.shp");
  plnCurveLength.WriteSetSHPFiles(wchAxesFileName0,&plnCurveLength,1);

  wcscpy(  wchAxesFileName0,  mwchOutPutFold);
  wcscat(wchAxesFileName0, L"\\plnDist.shp");
  plnDist.WriteSetSHPFiles(wchAxesFileName0,&plnDist,1);


  wcscpy(  wchAxesFileName0,  mwchOutPutFold);
  wcscat(wchAxesFileName0, L"\\plnPrfl.shp");
  plnPrfl.WriteSetSHPFiles(wchAxesFileName0,&plnPrfl,1);
  int ui = 0;
}
//-----------------------------------
double F(const double z, const double VAlza
         ,TTable_1D &tblPrfl,const double VAlCosTetta)
{
    const double VAlC0 = tblPrfl.calcValue(VAlza);

    int in0 = tblPrfl.getSegmentNum(z);
    const double VAlC1 = tblPrfl.mparrVal[in0];
            const double VAlC2= tblPrfl.mparrVal[in0 +1];
            const double VAlZn1 =tblPrfl.mparrArg[in0];
            const double VAlZn2 =tblPrfl.mparrArg[in0 +1];

    double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    double a= VAlC1 - b* VAlZn1;
    double valn = VAlC0/ (a + b * z);
    return 1./ sqrt(valn* valn -VAlCosTetta * VAlCosTetta);
}
double f(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    double a= VAlC1 - b* VAlZn1;
    double valn = VAlC0/ (a + b * z);
    return 1./ sqrt(valn* valn -VAlCosTetta * VAlCosTetta);
}
double f1(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    double a= VAlC1 - b* VAlZn1;
    double valn = VAlC0/ (a + b * z);
    return valn * valn/ pow(valn* valn -VAlCosTetta * VAlCosTetta, 1.5);
}
double f2(const double z, const double VAlC0,const double VAlZn1
          ,const double VAlZn2,const double VAlC1,const double VAlC2,const double VAlCosTetta)
{
    double b = (VAlC2 -VAlC1)/(VAlZn2 -VAlZn1);
    double a= VAlC1 - b* VAlZn1;
    double valn = VAlC0/ (a + b * z);
    return valn * valn/ sqrt(valn* valn -VAlCosTetta * VAlCosTetta);
}

void MainWindow::on_pushButton_5_clicked()
{
    QString strFold = QFileDialog::getOpenFileName(0,"Выбор .000 файла с профилем скорости","D:\\AKIN\\PROFILES" ,"*.000");
    this->ui->lineEdit_3->setText(strFold);

}

void MainWindow::on_pushButton_2_clicked()
{
    QString strFold = QFileDialog::getExistingDirectory(0,"Выбор_папки_для_хранения_графиков", "D:\\AKIN\\ОТЧЕТ_ПОЗИЦИОНИРОВАНИЕ");
    this->ui->lineEdit->setText(strFold);
    strFold.toWCharArray(mwchOutPutFold);
    mwchOutPutFold[strFold.size()] = 0;
}
