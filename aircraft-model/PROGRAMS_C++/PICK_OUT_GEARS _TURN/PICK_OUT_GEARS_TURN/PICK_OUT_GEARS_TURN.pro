#-------------------------------------------------
#
# Project created by QtCreator 2018-09-01T10:44:25
#
#-------------------------------------------------

QT       += core gui widgets

TARGET = Helic_Operation_v0
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
    ../../My_Qt_Lib/Helicopter/Blade.cpp \
    ../../My_Qt_Lib/Helicopter/Helic.cpp \
    ../../My_Qt_Lib/Helicopter/Rotor.cpp \
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
    ../../My_Qt_Lib/Helicopter/Environment.cpp \
    ../../My_Qt_Lib/Helicopter/AbstractGliderElement.cpp \
    ../../My_Qt_Lib/Helicopter/CommonGliderElement.cpp \
    ../../My_Qt_Lib/Helicopter/RuleGliderElementt.cpp \
    ../../My_Qt_Lib/YrWriteShapeFile.cpp \
    ../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.cpp \
    ../../My_Qt_Lib/liprog.cpp \
    ../../My_Qt_Lib/Helicopter/BallanceCalc.cpp \
    ../../My_Qt_Lib/Helicopter/PartHelicTraj.cpp \
    ../../My_Qt_Lib/Helicopter/HelicConstants.cpp \
    ../../My_Qt_Lib/Helicopter/TurnMove.cpp \
    ../../My_Qt_Lib/Helicopter/LineMove.cpp \
    ../../My_Qt_Lib/Helicopter/Rotating.cpp \
    ../../My_Qt_Lib/Helicopter/Hover.cpp

HEADERS += \
        mainwindow.h \
    ../../My_Qt_Lib/Helicopter/Blade.h \
    ../../My_Qt_Lib/Helicopter/Helic.h \
    ../../My_Qt_Lib/Helicopter/Rotor.h \
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
    ../../My_Qt_Lib/Helicopter/Environment.h \
    ../../My_Qt_Lib/Helicopter/AbstractGliderElement.h \
    ../../My_Qt_Lib/Helicopter/CommonGliderElement.h \
    ../../My_Qt_Lib/Helicopter/RuleGliderElement.h \
    ../../My_Qt_Lib/CalcCorMatrx.h \
    ../../My_Qt_Lib/YrWriteShapeFile.h \
    ../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.h \
    liprog.h \
    ../../My_Qt_Lib/Helicopter/BallanceCalc.h \
    ../../My_Qt_Lib/Helicopter/PartHelicTraj.h \
    ../../My_Qt_Lib/Helicopter/HelicConstants.h \
    ../../My_Qt_Lib/Helicopter/TurnMove.h \
    ../../My_Qt_Lib/Helicopter/LineMove.h \
    Rotating.h \
    ../../My_Qt_Lib/Helicopter/Rotating.h \
    ../../My_Qt_Lib/Helicopter/Hover.h

FORMS += \
        mainwindow.ui
# подключение папок с файлами !!!! 03.09.2018
INCLUDEPATH += ../../My_Qt_Lib

INCLUDEPATH += ../../My_Qt_Lib/SHFold_10Oct_15
INCLUDEPATH += ../../My_Qt_Lib/Helicopter



INCLUDEPATH += "C:/Program Files/MATLAB/R2018b/extern/include"
LIBS += -L"C:/Program Files/MATLAB/R2018b/extern/lib/win64/mingw64" -llibeng -llibmx

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
