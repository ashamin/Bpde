#include "bpdemainwindow.h"

#include <iostream>
#include <QDebug>

BpdeMainWindow::BpdeMainWindow(RInside& R, const QString& sourceFile, QObject *parent)
    : R(R),
    sourceFile(sourceFile), x(), y(), H(), solver(NULL)
{
    tempfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("tfile <- tempfile()")));
    svgfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("sfile <- tempfile()")));


//    R.parseEval("library(\"MASS\");"
//                "library(\"lattice\");"
//                "library(\"plyr\");"
//                "library(\"emdbook\");"
//                "library(\"rgl\");"
//                "library(\"fields\");"
//                "open3d();"
//                "bg3d(\"white\");"
//                "material3d(col=\"black\");"
//                "persp3d(x, y, H, aspect=c(1, 1, 0.5), col = \"lightblue\","
//                        "xlab = \"X\", ylab = \"Y\", zlab = \"Sinc( r )\");"
//                );

    QWidget *window = new QWidget;
    window->setWindowTitle("BpdeGUI");
    setCentralWidget(window);

    QGroupBox *runParameters = new QGroupBox("Run parameters");
    openMPEnabled = new QRadioButton("&OpenMP");
    openClEnabled = new QRadioButton("&OpenCL");

    openMPEnabled->setChecked(true);

    connect(openMPEnabled, SIGNAL(clicked()), this, SLOT(loadSource()));
    connect(openClEnabled, SIGNAL(clicked()), this, SLOT(loadSource()));

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(openMPEnabled);
    vbox->addWidget(openClEnabled);
//    runParameters->setMinimumSize(260,140);
//    runParameters->setMaximumSize(260,140);




    QLabel *threadsLabel = new QLabel("Threads");
    threadsLineEdit = new QLineEdit("4");
    QHBoxLayout *threadNumber = new QHBoxLayout;
    threadNumber->addWidget(threadsLabel);
    threadNumber->addWidget(threadsLineEdit);

    QHBoxLayout *deviceLayout = new QHBoxLayout;
    QLabel *deviceLabel = new QLabel("Device");
    deviceComboBox = new QComboBox();
    deviceComboBox->addItem("CPU Intel Core i5 1.7 Mhz");
    deviceLayout->addWidget(deviceLabel);
    deviceLayout->addWidget(deviceComboBox);

    QHBoxLayout* runLayout = new QHBoxLayout;
    runButton = new QPushButton("Run computations", this);
    qDebug() << "Connect : " <<
    connect(runButton, SIGNAL(clicked()), this, SLOT(solve()));
    runLayout->addWidget(runButton);

    QVBoxLayout* ulLayout = new QVBoxLayout;
    ulLayout->addLayout(vbox);
    ulLayout->addLayout(threadNumber);
    ulLayout->addLayout(deviceLayout);
    ulLayout->addLayout(runLayout);


    runParameters->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    runParameters->setLayout(ulLayout);




    QButtonGroup *kernelGroup = new QButtonGroup;
    kernelGroup->addButton(openMPEnabled, 0);
    kernelGroup->addButton(openClEnabled, 1);

    QGroupBox *solveParamBox = new QGroupBox("Density estimation bandwidth (scaled by 100)");

    QHBoxLayout *iterationsLayout = new QHBoxLayout;
    QHBoxLayout *stepLayout = new QHBoxLayout;
    QLabel *sourceFileLabel = new QLabel("SourceFile");
    sourceFileEdit = new QLineEdit(sourceFile);
//    connect(sourceFileEdit, SIGNAL(clicked()), this, SLOT(selectSourceFile()));
    iterationsEdit = new QLineEdit("10000");
    stepEdit = new QLineEdit("3600");
    QLabel *iterationsLabel = new QLabel("Iterations");
    QLabel *stepLabel = new QLabel("Step");
    QHBoxLayout *exportLayout = new QHBoxLayout;
    exportImage = new QPushButton("Export Isoterms");
    export3D = new QPushButton("Export 3D model");
    connect(export3D, SIGNAL(clicked()), this, SLOT(export3DModel()));

    iterationsLayout->addWidget(iterationsLabel);
    iterationsLayout->addWidget(iterationsEdit);
    stepLayout->addWidget(stepLabel);
    stepLayout->addWidget(stepEdit);


    exportLayout->addWidget(exportImage);
    exportLayout->addWidget(export3D);

    svg = new QSvgWidget();
    loadSource();

    QVBoxLayout *solveParamLayout = new QVBoxLayout;
    solveParamLayout->addWidget(sourceFileLabel);
    solveParamLayout->addWidget(sourceFileEdit);
    solveParamLayout->addLayout(iterationsLayout);
    solveParamLayout->addLayout(stepLayout);
    solveParamLayout->addLayout(exportLayout);

    solveParamBox->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    solveParamBox->setLayout(solveParamLayout);


    QHBoxLayout *upperlayout = new QHBoxLayout;
    upperlayout->addWidget(runParameters);
    upperlayout->addWidget(solveParamBox);

    QHBoxLayout *lowerlayout = new QHBoxLayout;
    lowerlayout->addWidget(svg);

    QVBoxLayout *outer = new QVBoxLayout;
    outer->addLayout(upperlayout);
    outer->addLayout(lowerlayout);
    window->setLayout(outer);
//    window->show();
}

void BpdeMainWindow::assignAreaToR(const Bpde::BArea &area)
{
    x.clear();
    y.clear();
    for (int i = 0; i<area.I; i++)
        x.push_back(area.x[i]);
    for (int j = 0; j<area.J; j++)
        y.push_back(area.y[j]);

    R.assign(area.I, "I");
    R.assign(area.J, "J");
    R.assign(x, "x");
    R.assign(y, "y");


    // wtf with Y??????
    R.parseEval("x = x[2:(length(x)-1)]");
    R.parseEval("y = y[2:(length(y)-1)]");

    reAssignH(area.H);
}

void BpdeMainWindow::reAssignH(double *Hfunc)
{
    H.clear();
    Bpde::BArea area(sourceFileEdit->text().toStdString());

    for (int i = 0; i<area.I * area.J; i++)
        H.push_back(Hfunc[i]);
    R.assign(H, "Htmp");
    R.parseEval("H = matrix(0, I-2, J-2)");
    R.parseEval("for (j in 2:(J-1)){ for (i in 2:(I-1)){ H[i-1,j-1] = Htmp[(j-1)*I+(i-1) + 1]} }");

//    R.parseEval("print(H)");
}

void BpdeMainWindow::loadSource()
{
    Bpde::BArea area(sourceFileEdit->text().toStdString());
    delete solver;
    if (openMPEnabled->isChecked()){
        solver = Bpde::BSolverBuilder::getInstance()->getSolver(
                sourceFileEdit->text().toStdString(), Bpde::ParallelizationMethod::OPENMP,
                threadsLineEdit->text().toInt());
    }
    else {
        solver = Bpde::BSolverBuilder::getInstance()->getSolver(
                sourceFileEdit->text().toStdString(), Bpde::ParallelizationMethod::OPENCL);
    }
    solver->setTimeStep(0);
    solver->addExtraIterations(-area.T);

    assignAreaToR(area);
    plot();
}

void BpdeMainWindow::selectSourceFile()
{
    QString file = "";
    file = QFileDialog::getOpenFileName(
                            NULL,
                            "Select file",
                            "/home",
                            "Bpde files (*.bde)");
    if (file != "") {
        sourceFileEdit->setText(file);
        loadSource();
    }
    else QMessageBox::information(NULL, "Error", "Maybe u have wrong file format");

}

void BpdeMainWindow::export3DModel()
{
    QString file = "";
    try {
        file = QFileDialog::getSaveFileName(
                NULL,
                "Export of 3D model",
                "/home",
                "Csvfiles (*.csv)");
    }catch(...)
    {
    }
    if (file != "") {
        R.assign(file.toStdString(), "csvFilePath");
        R.parseEval("fx = c(x, rep(-1, length(H)-length(x)));"
            "fy = c(y, rep(-1, length(H)-length(y)));"
            "fH = c(H);"
            "csv = data.frame(fx, fy, fH);"
            "write.csv2(csv, file=csvFilePath);"
        );
    }
    else {
        QMessageBox::information(this, "Error", "Cannot save 3D model");
    }
}

void BpdeMainWindow::plot()
{
    std::string cmd0 = "svg(width=6,height=6,pointsize=10,filename=tfile); ";
    std::string cmd = cmd0 + "library(\"fields\");image.plot(x, y, H);"
            "contour(x, y, H, add = TRUE);dev.off()";
    R.parseEvalQ(cmd);
    filterFile();
    svg->load(svgfile);
}

void BpdeMainWindow::solve()
{
    qDebug() << "start solve";
    setWindowTitle("Equation solving. Please wait ...");
    solver->addExtraIterations(iterationsEdit->text().toInt());
    solver->setTimeStep(stepEdit->text().toInt());
    reAssignH(solver->solve());
    plot();
//    QMessageBox::information(this, "Comp ended", QString("Comp time %1").
//                             arg(solver->exec_time()));
    setWindowTitle("BpdeGUI");
}

void BpdeMainWindow::filterFile() {
    QFile infile(tempfile);
    infile.open(QFile::ReadOnly);
    QFile outfile(svgfile);
    outfile.open(QFile::WriteOnly | QFile::Truncate);

    QTextStream in(&infile);
    QTextStream out(&outfile);
    QRegExp rx1("<symbol");
    QRegExp rx2("</symbol");
    while (!in.atEnd()) {
        QString line = in.readLine();
        line.replace(rx1, "<g"); // so '<symbol' becomes '<g ...'
        line.replace(rx2, "</g");// and '</symbol becomes '</g'
        out << line << "\n";
    }
    infile.close();
    outfile.close();
}
