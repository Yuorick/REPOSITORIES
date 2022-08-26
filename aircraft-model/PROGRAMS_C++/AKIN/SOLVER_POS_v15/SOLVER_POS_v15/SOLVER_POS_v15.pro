
QT       += core gui widgets

TARGET = SOLVER_POS_v15
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
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPointXY.cpp \
    ../../../My_Qt_Lib/MatrixProccess.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URFigure.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/DBFTab.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/DBFTable.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Triang.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URMultiPoint.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolyLine.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolyLineZ.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrLinTransform.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrRastr.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrRead.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrWrite.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPointZ.cpp \
    ../../../My_Qt_Lib/Equations.cpp \
    ../../../My_Qt_Lib/Comp.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolygon.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/ArcParab.cpp \
    ../../../My_Qt_Lib/YrWriteShapeFile.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Line.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Segment.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Plane.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Lnsgm.cpp \
    ../../../My_Qt_Lib/Gauss.cpp \   
    ../../../My_Qt_Lib/Table_1D.cpp \    
    ../../../My_Qt_Lib/SUBWATER/Gps.cpp \
    ../../../My_Qt_Lib/SUBWATER/HidroRLS.cpp \    
    ../../../My_Qt_Lib/SUBWATER/PeaceSins.cpp \
    ../../../My_Qt_Lib/SUBWATER/PeaceVess.cpp \
    ../../../My_Qt_Lib/SUBWATER/Platform.cpp \
    ../../../My_Qt_Lib/SUBWATER/PosSolver.cpp \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl2D.cpp \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl3D.cpp \
    ../../../My_Qt_Lib/SUBWATER/SubWaterBeam.cpp \
    ../../../My_Qt_Lib/SUBWATER/TrueMeasParams.cpp \
    ../../../My_Qt_Lib/CoordSystTrsf.cpp \
    ../../../My_Qt_Lib/SUBWATER/BigMeasure.cpp \
    ../../../My_Qt_Lib/SUBWATER/LblSolver.cpp \    
    ../../../My_Qt_Lib/Environment.cpp \
    mythread.cpp \
    ../../../My_Qt_Lib/SUBWATER/Solver2D_Angs.cpp \
    ../../../My_Qt_Lib/SUBWATER/Solver3D_Angs.cpp \
    ../../../My_Qt_Lib/SUBWATER/acdFile.cpp \
    ../../../My_Qt_Lib/SHFold_10Oct_15/ArcEllipse.cpp \
    ../../../My_Qt_Lib/SUBWATER/PntCloud.cpp





HEADERS += \
        mainwindow.h \   
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPointXY.h \
    ../../../My_Qt_Lib/MatrixProccess.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URFigure.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/DBFTab.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/DBFTable.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Triang.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URMultiPoint.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolyLine.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrLinTransform.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrRastr.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrRead.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/YrWrite.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPointZ.h \
    ../../../My_Qt_Lib/Equations.h \
    ../../../My_Qt_Lib/Comp.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolygon.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/ArcParab.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/URPolyLineZ.h \
    ../../../My_Qt_Lib/YrWriteShapeFile.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/LongPlane.h \   
    ../../../My_Qt_Lib/SHFold_10Oct_15/Segment.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/Lnsgm.h \    
    ../../../My_Qt_Lib/Table_1D.h \   
    ../../../My_Qt_Lib/SUBWATER/Gps.h \
    ../../../My_Qt_Lib/SUBWATER/HidroRLS.h \    
    ../../../My_Qt_Lib/SUBWATER/PeaceSins.h \
    ../../../My_Qt_Lib/SUBWATER/PeaceVess.h \
    ../../../My_Qt_Lib/SUBWATER/Platform.h \
    ../../../My_Qt_Lib/SUBWATER/PosSolver.h \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl2D.h \
    ../../../My_Qt_Lib/SUBWATER/Rls_Usbl3D.h \
    ../../../My_Qt_Lib/SUBWATER/SubWaterBeam.h \
    ../../../My_Qt_Lib/SUBWATER/TrueMeasParams.h \
    ../../../My_Qt_Lib/CoordSystTrsf.h \
    ../../../My_Qt_Lib/SUBWATER/BigMeasure.h \
    ../../../My_Qt_Lib/SUBWATER/LblSolver.h \    
    ../../../My_Qt_Lib/Environment.h \
    mythread.h \
    ../../../My_Qt_Lib/SUBWATER/Solver2D_Angs.h \
    ../../../My_Qt_Lib/SUBWATER/Solver3D_Angs.h \
    ../../../My_Qt_Lib/SUBWATER/acdFile.h \
    ../../../My_Qt_Lib/SHFold_10Oct_15/ArcEllipse.h \
    ../../../My_Qt_Lib/SUBWATER/PntCloud.h




FORMS += \
        mainwindow.ui
# подключение папок с файлами !!!! 03.09.2018
INCLUDEPATH += ../../../My_Qt_Lib
INCLUDEPATH += ../../../My_Qt_Lib/SHFold_10Oct_15
INCLUDEPATH += ../../../My_Qt_Lib/SUBWATER


# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target



