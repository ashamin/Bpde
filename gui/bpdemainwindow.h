#ifndef BPDEMAINWINDOW_H
#define BPDEMAINWINDOW_H

#include <vector>
#include <iostream>
#include <fstream>

#include <QWidget>
#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QGroupBox>
#include <QButtonGroup>
#include <QMainWindow>
#include <QHBoxLayout>
#include <QSlider>
#include <QSpinBox>
#include <QLabel>
#include <QTemporaryFile>
#include <QSvgWidget>
#include <QComboBox>
#include <QPushButton>
#include <QTextStream>
#include <QFileDialog>
#include <QSvgRenderer>
#include <QPainter>
#include <QImage>
#include <QMessageBox>

#include <RInside.h>
#include "BBuilder.h"

class BpdeMainWindow : public QMainWindow
{
    Q_OBJECT
public:
    BpdeMainWindow(RInside& R, const QString& sourceFile, QObject * parent = NULL);

private:
    void plot();
    void filterFile();
    void assignAreaToR(const Bpde::BArea& area);
    void reAssignH(double* HFunc);
    void scanDevices();
    std::string getSource(std::string filename);

    Bpde::BSolver *solver;

    RInside& R;
    QSvgWidget *svg;
    QString tempfile;
    QString svgfile;
    QString sourceFile;

    std::vector<double> x, y, H;

    QRadioButton *openMPEnabled;
    QRadioButton *openClEnabled;
    QPushButton *runButton, *exportImage, *export3D;

    QLineEdit *sourceFileEdit, *iterationsEdit, *stepEdit, *threadsLineEdit;

    QComboBox *deviceComboBox;

    std::vector<cl::Device> devices;

public slots:
    void solve();
    void loadSource();
    void selectSourceFile();
    void export3DModel();
    void exportIsoterms();

};

#endif // BPDEMAINWINDOW_H

