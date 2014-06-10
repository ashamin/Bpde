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
                                "Выберите файл",
                                "/home",
                                "Bpde files (*.bde)");
        if (file != "")
            window = new BpdeMainWindow(R, file);
        else
            QMessageBox::information(NULL, "Ошибка!", "Во время работы произошел сбой.");
    }catch(...) {
        QMessageBox::information(NULL, "Ошибка", "Файл поврежден или имеет неверный формат");
    }
    if (window)
        window->show();

    int ret = a.exec();
    delete window;
    return ret;
}
