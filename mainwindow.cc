#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <string>
using namespace std;

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::updateLegend(){

}

void MainWindow::on_radioButton_clicked()
{
    //black and white
    ui->openGLLegend->scalar_col= 0;
    ui->openGLLegend->update();

    ui->openGLWidget->scalar_col = 0;
    ui->openGLWidget->update();
}

void MainWindow::on_radioButton_2_clicked()
{
    //rainbow
    ui->openGLLegend->scalar_col= 1;
    ui->openGLLegend->update();

    ui->openGLWidget->scalar_col = 1;
    ui->openGLWidget->update();
}

void MainWindow::on_radioButton_3_clicked()
{
    //heatmap
    ui->openGLLegend->scalar_col= 2;
    ui->openGLLegend->update();

    ui->openGLWidget->scalar_col = 2;
    ui->openGLWidget->update();
}

void MainWindow::on_value_valueChanged(int value)
{
    ui->openGLLegend->increase_v = (float)(value/10);
    ui->openGLLegend->update();

    ui->openGLWidget->increase_v = (float)(value/10);
    ui->openGLWidget->update();
}

void MainWindow::on_value_2_valueChanged(int value)
{
    ui->openGLLegend->increase_s = (float)(value/10);
    ui->openGLLegend->update();

    ui->openGLWidget->increase_s = (float)(value/10);
    ui->openGLWidget->update();
}

void MainWindow::on_rho_clicked()
{
    ui->openGLWidget->dataset = 'r';
    ui->openGLWidget->update();
}

void MainWindow::on_v_clicked()
{
    ui->openGLWidget->dataset = 'v';
    ui->openGLWidget->update();
}

void MainWindow::on_f_clicked()
{
    ui->openGLWidget->dataset = 'f';
    ui->openGLWidget->update();
}


void MainWindow::updateLabel(float lmax, float lmin)
{
   // float lmax = ui->openGLWidget->d_simulation.get_rho_max();
    float unit = (lmax-lmin)/10;
    QString lbel = QString::number(lmax-unit);
    ui->lmax->setText(lbel);
     lbel = QString::number(lmax-unit*2);
    ui->lmax_2->setText(lbel);
     lbel = QString::number(lmax-unit*3);
    ui->lmax_3->setText(lbel);
     lbel = QString::number(lmax-unit*4);
    ui->lmax_4->setText(lbel);
     lbel = QString::number(lmax-unit*5);
    ui->lmax_5->setText(lbel);
     lbel = QString::number(lmax-unit*6);
    ui->lmax_6->setText(lbel);
     lbel = QString::number(lmax-unit*7);
    ui->lmax_7->setText(lbel);
     lbel = QString::number(lmax-unit*8);
    ui->lmax_8->setText(lbel);
     lbel = QString::number(lmax-unit*9);
    ui->lmax_9->setText(lbel);
     lbel = QString::number(lmin);
    ui->lmin->setText(lbel);
}


void MainWindow::on_draw_vecs_clicked()
{
    ui->openGLWidget->draw_vecs = 1;
    ui->openGLWidget->draw_smoke = 0;
    ui->openGLWidget->update();
}

void MainWindow::on_no_draw_vecs_clicked()
{
    ui->openGLWidget->draw_vecs = 0;
    ui->openGLWidget->draw_smoke = 1;
    ui->openGLWidget->update();
}

void MainWindow::on_scalar_rho_clicked()
{
    ui->openGLWidget->scalar_dataset = 'r';
    ui->openGLWidget->update();
}

void MainWindow::on_scalar_v_clicked()
{
    ui->openGLWidget->scalar_dataset = 'v';
    ui->openGLWidget->update();
}

void MainWindow::on_scalar_f_clicked()
{
    ui->openGLWidget->scalar_dataset = 'f';
    ui->openGLWidget->update();
}

void MainWindow::on_vector_v_clicked()
{
    ui->openGLWidget->vector_dataset = 'v';
    ui->openGLWidget->update();
}

void MainWindow::on_vector_f_clicked()
{
    ui->openGLWidget->vector_dataset = 'f';
    ui->openGLWidget->update();
}
