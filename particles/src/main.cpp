#include <QApplication>
#include "mainwindow.h"

int framenumber;

int main(int argc, char *argv[]) {
    framenumber = atoi(argv[1]);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    return a.exec();
}
