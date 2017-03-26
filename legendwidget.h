#ifndef LEGENDWIDGET_H
#define LEGENDWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMouseEvent>
#include <rfftw.h>             //the numerical simulation FFTW library
#include <fftw.h>             //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <cmath>               //for various math functions
#include <vector>
#include "simulation.h"


class LegendWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

    Simulation d_simulation;

    const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, banded
    const int COLOR_RAINBOW=1;
    const int COLOR_YELLOW=2;
    const int COLOR_BANDS=3;

    int DIM;
    int color_dir;            //use direction color-coding or not

    int winWidth, winHeight;      //size of the graphics window, in pixels

public:

    int scalar_col;           //method for scalar coloring
    float increase_s;
    float increase_v;

    LegendWidget(QWidget *parent = 0);

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    void rainbow(float value,float* R,float* G,float* B);
    void saturation(float value,float* R,float* G,float* B);
    void heatmap(float value, float* R, float* G, float* B);
    void rgb2hsv(float R, float G, float B, float* S, float* V, float* H);
    void hsv2rgb(float* R, float* G, float* B, float S, float V, float H);
    void set_colormap(float vy);

    void visualize(void);

public slots:
   void do_one_simulation_step();
private:

};

#endif // LEGENDWIDGET_H
