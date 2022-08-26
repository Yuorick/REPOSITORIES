#-------------------------------------------------
#
# Project created by QtCreator 2020-05-03T08:03:20
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = SYSTEM_ERRS_CALCULATOR
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

CONFIG += c++11

SOURCES += \
        main.cpp \
        mainwindow.cpp \
    ../../My_Qt_Lib/CalcCorMatrx.cpp \
    ../../My_Qt_Lib/Comp.cpp \
    ../../My_Qt_Lib/Equations.cpp \
    ../../My_Qt_Lib/MatrixProccess.cpp \
    ../../My_Qt_Lib/MinSquare.cpp \
    ../../My_Qt_Lib/Gauss.cpp \
    ../../My_Qt_Lib/SPATIAL_PARAMS_COORDINATION/Adjustment.cpp

HEADERS += \
        mainwindow.h \
    ../../My_Qt_Lib/CalcCorMatrx.h \
    ../../My_Qt_Lib/Comp.h \
    ../../My_Qt_Lib/Equations.h \
    ../../My_Qt_Lib/MatrixProccess.h \    
    ../../My_Qt_Lib/SPATIAL_PARAMS_COORDINATION/Adjustment.h




FORMS += \
        mainwindow.ui

INCLUDEPATH += ../../My_Qt_Lib
INCLUDEPATH += ../../My_Qt_Lib/SPATIAL_PARAMS_COORDINATION

