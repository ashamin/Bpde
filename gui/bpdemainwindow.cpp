#include "bpdemainwindow.h"

#include <iostream>
#include <QDebug>

BpdeMainWindow::BpdeMainWindow(RInside& R, const QString& sourceFile, QObject *parent)
    : R(R),
    sourceFile(sourceFile), x(), y(), H(), solver(NULL)
{
    tempfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("tfile <- tempfile()")));
    svgfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("sfile <- tempfile()")));

    QWidget *window = new QWidget;
    window->setWindowTitle("BpdeGUI");
    setCentralWidget(window);

    QGroupBox *runParameters = new QGroupBox("Параметры запуска");
    openMPEnabled = new QRadioButton("&OpenMP");
    openClEnabled = new QRadioButton("&OpenCL");

    openMPEnabled->setChecked(true);

    connect(openMPEnabled, SIGNAL(clicked()), this, SLOT(loadSource()));
    connect(openClEnabled, SIGNAL(clicked()), this, SLOT(loadSource()));

    QVBoxLayout *vbox = new QVBoxLayout;
    vbox->addWidget(openMPEnabled);
    vbox->addWidget(openClEnabled);


    QLabel *threadsLabel = new QLabel("Количество потоков");
    threadsLineEdit = new QLineEdit("4");
    QHBoxLayout *threadNumber = new QHBoxLayout;
    threadNumber->addWidget(threadsLabel);
    threadNumber->addWidget(threadsLineEdit);

    QHBoxLayout *deviceLayout = new QHBoxLayout;
    QLabel *deviceLabel = new QLabel("Устройство");
    deviceComboBox = new QComboBox();

    scanDevices();
    for (std::vector<cl::Device>::iterator it = devices.begin(); it != devices.end(); it++)
        deviceComboBox->addItem((*it).getInfo<CL_DEVICE_NAME>().c_str());

    connect(deviceComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(loadSource()));

    deviceLayout->addWidget(deviceLabel);
    deviceLayout->addWidget(deviceComboBox);

    QHBoxLayout* runLayout = new QHBoxLayout;
    runButton = new QPushButton("Начать вычисления", this);
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

    QGroupBox *solveParamBox = new QGroupBox("Настройки вычислительного метода");

    QHBoxLayout *iterationsLayout = new QHBoxLayout;
    QHBoxLayout *stepLayout = new QHBoxLayout;
    QLabel *sourceFileLabel = new QLabel("SourceFile");
    sourceFileEdit = new QLineEdit(sourceFile);
    iterationsEdit = new QLineEdit("10000");
    stepEdit = new QLineEdit("3600");
    QLabel *iterationsLabel = new QLabel("Итерации");
    QLabel *stepLabel = new QLabel("Шаг            ");
    QHBoxLayout *exportLayout = new QHBoxLayout;
    exportImage = new QPushButton("Экспорт изотерм");
    export3D = new QPushButton("Экспорт 3D модели");
    connect(exportImage, SIGNAL(clicked()), this, SLOT(exportIsoterms()));
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
}

void BpdeMainWindow::scanDevices()
{
    std::vector<cl::Platform> platforms;
    std::vector<cl::Device> tmp;
    cl::Platform::get(&platforms);
    for (std::vector<cl::Platform>::iterator it = platforms.begin();
         it!=platforms.end(); it++)
    {
        (*it).getDevices(CL_DEVICE_TYPE_ALL, &tmp);
        std::copy(tmp.begin(), tmp.end(), std::back_inserter(devices));
    }
}

std::string BpdeMainWindow::getSource(std::string filename)
{
    std::ifstream file(filename.c_str());
    std::string ret(std::istreambuf_iterator<char>(file), (std::istreambuf_iterator<char>()));
    return ret;
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
        std::vector<cl::Device> dev;
        dev.push_back(devices[deviceComboBox->currentIndex()]);
        solver = Bpde::BSolverBuilder::getInstance()->getSolver(
                    sourceFileEdit->text().toStdString(), Bpde::ParallelizationMethod::OPENCL,
                    dev);
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
                            "Выберите файл",
                            "/home",
                            "Bpde files (*.bde)");
    if (file != "") {
        sourceFileEdit->setText(file);
        loadSource();
    }
    else QMessageBox::information(NULL, "Ошибка", "Файл поврежден или имеет неверный формат");
}

void BpdeMainWindow::export3DModel()
{
    QString file = "";
    try {
        file = QFileDialog::getSaveFileName(
                NULL,
                "Экспорт 3D модели (csv формат)",
                "/home",
                "Csvfiles (*.csv)");
    }catch(...)
    {
    }
    if (file != "") {
        R.assign(file.toStdString(), "csvFilePath");
        std::string cmd = getSource("/home/ashamin/diploma/src/Bpde/gui/R/export3dmodel.R");
        R.parseEval(cmd);
    }
}

void BpdeMainWindow::exportIsoterms()
{
    QString file = "";
    try {
        file = QFileDialog::getSaveFileName(
                NULL,
                "Экпорт изотерм (формат jpg)",
                "/home",
                "Images (*.jpg)");
    }catch(...)
    {
    }
    if (file != "") {
        QSvgRenderer renderer(svgfile);
        QImage image(x.size()*40+(int)(x.size()*4), y.size()*40, QImage::Format_ARGB32);
        image.fill(0xaaA08080);

        QPainter painter(&image);
        renderer.render(&painter);

        image.save(file);
    }
}

void BpdeMainWindow::plot()
{
    std::string cmd = getSource("/home/ashamin/diploma/src/Bpde/gui/R/isoterms.R");
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
        line.replace(rx1, "<g");
        line.replace(rx2, "</g");
        out << line << "\n";
    }
    infile.close();
    outfile.close();
}
