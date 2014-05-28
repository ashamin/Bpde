#include "bpdemainwindow.h"
#include <QApplication>

#include <vector>

int main(int argc, char *argv[])
{
    RInside R(argc, argv);
//    std::string str = "print(x)";

//    std::vector<double> x;
//    x.push_back(1);
//    x.push_back(8);

//    R.assign(x, "x");

//    R.parseEval(str);

    QApplication a(argc, argv);

    QString file = "";
    QMainWindow * window = nullptr;
    try {
        file = QFileDialog::getOpenFileName(
                                NULL,
                                "Select file",
                                "/home",
                                "Bpde files (*.bde)");
        if (file != "")
            window = new BpdeMainWindow(R, file);
//            BpdeMainWindow window(R, file);
        else
            QMessageBox::information(NULL, "Hello World!", "Hi!");
    }catch(...) {
        QMessageBox::information(NULL, "Error", "Maybe u have wrong file format");
    }
    if (window)
        window->show();

    int ret = a.exec();
    delete window;
    return ret;
}
