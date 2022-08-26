
QT       += core gui widgets

TARGET = MEASURE_TRANSFORM_v1
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
    ../../../My_Qt_Lib/MatrixProccess.cpp \
    ../../../My_Qt_Lib/Equations.cpp \
    ../../../My_Qt_Lib/Comp.cpp \    
    ../../../My_Qt_Lib/Gauss.cpp \   
    ../../../My_Qt_Lib/Table_1D.cpp \    
    ../../../My_Qt_Lib/SUBWATER/Gps.cpp \
    ../../../My_Qt_Lib/SUBWATER/HidroRLS.cpp \
    ../../../My_Qt_Lib/SUBWATER/ImitMod.cpp \
    ../../../My_Qt_Lib/SUBWATER/PeaceSins.cpp \
    ../../../My_Qt_Lib/SUBWATER/PeaceVess.cpp \
    ../../../My_Qt_Lib/SUBWATER/Platform.cpp \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl2D.cpp \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl3D.cpp \
    ../../../My_Qt_Lib/SUBWATER/SubWaterBeam.cpp \
    ../../../My_Qt_Lib/SUBWATER/TrueMeasParams.cpp \
    ../../../My_Qt_Lib/CoordSystTrsf.cpp \    
    ../../../My_Qt_Lib/SUBWATER/BigMeasure.cpp \
    ../../../My_Qt_Lib/Environment.cpp




HEADERS += \
        mainwindow.h \
    ../../../My_Qt_Lib/MatrixProccess.h \
    ../../../My_Qt_Lib/Equations.h \
    ../../../My_Qt_Lib/Comp.h \
    ../../../My_Qt_Lib/Table_1D.h \   
    ../../../My_Qt_Lib/SUBWATER/Gps.h \
    ../../../My_Qt_Lib/SUBWATER/HidroRLS.h \
    ../../../My_Qt_Lib/SUBWATER/ImitMod.h \
    ../../../My_Qt_Lib/SUBWATER/PeaceSins.h \
    ../../../My_Qt_Lib/SUBWATER/PeaceVess.h \
    ../../../My_Qt_Lib/SUBWATER/Platform.h \    
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl2D.h \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl3D.h \
    ../../../My_Qt_Lib/SUBWATER/SubWaterBeam.h \
    ../../../My_Qt_Lib/SUBWATER/TrueMeasParams.h \
    ../../../My_Qt_Lib/CoordSystTrsf.h \    
    ../../../My_Qt_Lib/SUBWATER/BigMeasure.h \
    ../../../My_Qt_Lib/Environment.h



FORMS += \
        mainwindow.ui

QMAKE_CXXFLAGS += -fopenmp
LIBS += -fopenmp
# подключение папок с файлами !!!! 03.09.2018
INCLUDEPATH += ../../../My_Qt_Lib
INCLUDEPATH += ../../../My_Qt_Lib/SUBWATER
INCLUDEPATH += ../../../My_Qt_Lib/Helicopter

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target



