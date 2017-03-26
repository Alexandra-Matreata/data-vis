#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void updateLegend();

private slots:
    void on_radioButton_clicked();

    void on_radioButton_2_clicked();

    void on_radioButton_3_clicked();

    void on_value_valueChanged(int value);

    void on_value_2_valueChanged(int value);

    void on_rho_clicked();

    void on_v_clicked();

    void on_f_clicked();

    void updateLabel(float lmax, float lmin);

    void on_draw_vecs_clicked();

    void on_no_draw_vecs_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
