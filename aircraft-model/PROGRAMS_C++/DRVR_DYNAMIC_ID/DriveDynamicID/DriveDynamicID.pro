
QT       += core gui widgets

TARGET = ElDr_RESEARCH
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11 Q0

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPointXY.cpp \
    ../../My_Qt_Lib/MatrixProccess.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URFigure.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/DBFTab.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/DBFTable.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/Triang.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URMultiPoint.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolyLine.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolyLineZ.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrLinTransform.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrRastr.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrRead.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrWrite.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPointZ.cpp \
    ../../My_Qt_Lib/Equations.cpp \
    ../../My_Qt_Lib/Comp.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolygon.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/ArcParab.cpp \
    ../../My_Qt_Lib/YrWriteShapeFile.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/ElectDriver.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/ElectMotor.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/Load.cpp \   
    ../../My_Qt_Lib/ELECT_DRIVER/LinOptCtrlSyst.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/Line.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/LinOptCtrlSyst_dim21.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/Brls.cpp \
    ../../My_Qt_Lib/Helicopter/Environment.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/Segment.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/Plane.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/Lnsgm.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlVelo.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlVeloQuiq.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/StatSolutionParams.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/Ctrl.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlPos.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlFollow.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/MeasurmentImitator.cpp \
    ../../My_Qt_Lib/Gauss.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/Filtr.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/DriveMoveImit.cpp \
    ../../My_Qt_Lib/LinDiffEq.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/ParamsID.cpp \
    ../../My_Qt_Lib/MinSquare.cpp \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlTrad.cpp

HEADERS += \
        mainwindow.h \   
    ../../My_Qt_Lib/SHFold_10Oct_15/URPointXY.h \
    ../../My_Qt_Lib/MatrixProccess.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URFigure.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/DBFTab.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/DBFTable.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/Triang.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URMultiPoint.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolyLine.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrLinTransform.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrRastr.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrRead.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/YrWrite.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPointZ.h \
    ../../My_Qt_Lib/Equations.h \
    ../../My_Qt_Lib/Comp.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolygon.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/ArcParab.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/URPolyLineZ.h \
    ../../My_Qt_Lib/CalcCorMatrx.h \
    ../../My_Qt_Lib/YrWriteShapeFile.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.h \
    ../../My_Qt_Lib/ELECT_DRIVER/ElectDriver.h \
    ../../My_Qt_Lib/ELECT_DRIVER/ElectMotor.h \
    ../../My_Qt_Lib/ELECT_DRIVER/Load.h \
    ../../My_Qt_Lib/LinDiffEq.h \
    ../../My_Qt_Lib/ELECT_DRIVER/LinOptCtrlSyst.h \
    ../../My_Qt_Lib/ELECT_DRIVER/LinOptCtrlSyst_dim21.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/Segment.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/Lnsgm.h \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlVelo.h \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlVeloQuiq.h \
    ../../My_Qt_Lib/ELECT_DRIVER/StatSolutionParams.h \
    ../../My_Qt_Lib/ELECT_DRIVER/Ctrl.h \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlPos.h \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlFollow.h \
    ../../My_Qt_Lib/ELECT_DRIVER/MeasurmentImitator.h \
    ../../My_Qt_Lib/ELECT_DRIVER/Filtr.h \
    ../../My_Qt_Lib/ELECT_DRIVER/DriveMoveImit.h \
    ../../My_Qt_Lib/ELECT_DRIVER/ParamsID.h \
    ../../My_Qt_Lib/MinSquare.h \
    ../../My_Qt_Lib/ELECT_DRIVER/CtrlTrad.h
    ../../My_Qt_Lib/ELECT_DRIVER/Load.h

FORMS += \
        mainwindow.ui
# подключение папок с файлами !!!! 03.09.2018
INCLUDEPATH += ../../My_Qt_Lib

INCLUDEPATH += ../../My_Qt_Lib/SHFold_10Oct_15
INCLUDEPATH += ../../My_Qt_Lib/ELECT_DRIVER
INCLUDEPATH += ../../My_Qt_Lib/Helicopter




# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target


