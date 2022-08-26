#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    // QApplication a(argc, argv);
    // MainWindow w;
    // w.show();

    // return a.exec();
     QCoreApplication::addLibraryPath(".");
     QCoreApplication::addLibraryPath("plugins//imageformats");
     QCoreApplication::addLibraryPath("imageformats");
     QCoreApplication::addLibraryPath("plugins//platforms");
     QCoreApplication::addLibraryPath("platforms");

      QApplication a(argc, argv);
      MainWindow w;
      w.show();

      return a.exec();
}
