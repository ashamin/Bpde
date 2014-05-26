#include "bpdemainwindow.h"
#include <QApplication>

#include <QFileDialog>
#include <QMessageBox>

int main(int argc, char *argv[])
{
    RInside R(argc, argv);
    QApplication a(argc, argv);

    QString file = "";
    file = QFileDialog::getOpenFileName(
                            NULL,
                            "Select file",
                            "/home",
                            "Bpde files (*.bde)");
    if (file != "")
        BpdeMainWindow window(R, file);
    else
        QMessageBox::information(NULL, "Hello World!", "Hi!");

    return a.exec();
}
