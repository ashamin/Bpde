#ifndef BPDEMAINWINDOW_H
#define BPDEMAINWINDOW_H

#include <vector>

#include <QMainWindow>
#include <QSvgWidget>
#include <QHBoxLayout>
#include <QFile>
#include <QTextStream>

#include <RInside.h>
#include "BBuilder.h"

class BpdeMainWindow : public QMainWindow
{
    Q_OBJECT
public:
    BpdeMainWindow(RInside& R, QString sourceFile);
private:
    void plot();
    void filterFile();

    RInside& R;
    QSvgWidget *svg;
    QString tempfile;
    QString svgfile;
    QString sourceFile;

    std::vector<double> x, y, H;
};

#endif // BPDEMAINWINDOW_H

