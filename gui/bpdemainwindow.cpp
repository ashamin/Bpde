#include "bpdemainwindow.h"

#include <iostream>

BpdeMainWindow::BpdeMainWindow(RInside& R, QString sourceFile)
    : R(R),
    sourceFile(sourceFile)
{
    tempfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("tfile <- tempfile()")));
    svgfile = QString::fromStdString(Rcpp::as<std::string>(R.parseEval("sfile <- tempfile()")));

    QWidget *window = new QWidget;
    window->setWindowTitle("Bpde GUI");

    svg = new QSvgWidget;
    plot();

    QHBoxLayout *layout = new QHBoxLayout;
    layout->addWidget(svg);
    window->setLayout(layout);
    window->show();
}

void BpdeMainWindow::plot()
{
    //x <- seq(-10, 10, length= 30);y <- x;f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r };z <- outer(x, y, f);z[is.na(z)] <- 1;
    std::string cmd0 = "svg(width=6,height=6,pointsize=10,filename=tfile); ";
    std::string cmd = cmd0 + "library(\"fields\");image.plot(volcano);contour(volcano, add = TRUE);dev.off()";
    R.parseEvalQ(cmd);
    filterFile();           	// we need to simplify the svg file for display by Qt
    svg->load(svgfile);
}

void BpdeMainWindow::filterFile() {
    // cairoDevice creates richer SVG than Qt can display
    // but per Michaele Lawrence, a simple trick is to s/symbol/g/ which we do here
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
