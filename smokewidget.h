#ifndef SMOKEWIDGET_H
#define SMOKEWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMouseEvent>
#include <vector>
#include "simulation.h"

QT_FORWARD_DECLARE_CLASS(QOpenGLShaderProgram)

class SmokeWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

    const int COLOR_BLACKWHITE=0;   //different types of color mapping: black-and-white, rainbow, banded
    const int COLOR_RAINBOW=1;
    const int COLOR_YELLOW=2;
    const int COLOR_BANDS=3;

    int DIM;
    int   winWidth, winHeight;      //size of the graphics window, in pixels
    int   color_dir;            //use direction color-coding or not
    float vec_scale;			//scaling of hedgehogs
   // int   scalar_col;           //method for scalar coloring
    int   frozen;               //toggles on/off the animation

    QPoint lastPos;
    bool select_points;
    std::vector<int> mouse_x;
    std::vector<int> mouse_y;

public:

    Simulation d_simulation;
    int scalar_col;
    float increase_s;
    float increase_v;
    char dataset;
    char scalar_dataset;
    char glyphs;
    char vector_dataset;
    int   draw_smoke;           //draw the smoke or not
    int   draw_vecs;            //draw the vector field or not

    SmokeWidget(QWidget *parent = 0);

    void initializeGL();
    void paintGL();
    void resizeGL(int width, int height);

    void rainbow(float value,float* R,float* G,float* B);
    void saturation(float value,float* R,float* G,float* B);
    void heatmap(float value, float* R, float* G, float* B);
    void rgb2hsv(float R, float G, float B, float* S, float* V, float* H);
    void hsv2rgb(float* R, float* G, float* B, float S, float V, float H);

    void set_colormap(float vy);
    void direction_to_color(float x, float y, int method);
    void length_to_color(fftw_real *x, fftw_real *y, int method);

    void visualize(void);
    void display(void);
    void reshape(int w, int h);

    void updateLabel();

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

signals:
   void upd(float lmax, float lmin);

public slots:
   void do_one_simulation_step();
private:


};

#endif // SMOKEWIDGET_H
