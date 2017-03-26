#include "legendwidget.h"
#include <GL/glut.h>
#include <QTimer>
#include <QObject>
#include <iostream>


LegendWidget::LegendWidget(QWidget *parent)
    : QOpenGLWidget(parent)
{
    DIM = 50;
    color_dir = 0;
    scalar_col = 0;
    increase_s = 0.8;
    increase_v = 0.8;

    QTimer *timer = new QTimer;
    timer->start(1);
    QObject::connect(timer,SIGNAL(timeout()),this,SLOT(do_one_simulation_step()));
}


void LegendWidget::initializeGL()
{
    initializeOpenGLFunctions();
}

void LegendWidget::paintGL()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_TABLE);
    glEnable(GL_SMOOTH);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable (GL_BLEND);
    glClearColor(0.0,0.0,0.0,0.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    visualize();
    glFlush();
}

void LegendWidget::resizeGL(int w, int h)
{
    glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
    winWidth = w; winHeight = h;

}

void LegendWidget::do_one_simulation_step()
{
    update();
}

void LegendWidget::visualize()
{
    int        i, j;
    fftw_real  wn = winWidth;   // Grid cell width
    fftw_real  hn = (fftw_real)winHeight / (fftw_real)(winHeight);  // Grid cell heigh

        float idx0, idx1, idx2, idx3;
        double px0, py0, px1, py1, px2, py2, px3, py3;
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBegin(GL_TRIANGLES);

            for (j = 0; j < winHeight; j++)
            {
                px0 = 0;
                py0 = (fftw_real)j * hn;
                idx0 = (float)j/(float)winHeight;

                px1 = 0;
                py1 = (fftw_real)(j + 1) * hn;
                idx1 = (float)j/(float)winHeight;


                px2 = wn;
                py2 = (fftw_real)(j + 1) * hn;
                idx2 = (float)j/(float)winHeight;


                px3 = wn;
                py3 = (fftw_real)j * hn;
                idx3 = (float)j/(float)winHeight;


                set_colormap(idx0);    glVertex3f(px0, py0, 0.0);
                set_colormap(idx1);    glVertex3f(px1, py1, 0.0);
                set_colormap(idx2);    glVertex3f(px2, py2, 0.0);


                set_colormap(idx0);    glVertex3f(px0, py0, 0.0);
                set_colormap(idx2);    glVertex3f(px2, py2, 0.0);
                set_colormap(idx3);    glVertex3f(px3, py3, 0.0);
            }
        glEnd();
}

void LegendWidget::rainbow(float value, float *R, float *G, float *B)
{
    const float dx=0.8;
    if (value<0)
        value=0;
    if (value>1)
        value=1;
    value = (6-2*dx)*value+dx;
    *R = d_simulation.max(0.0,(3-fabs(value-4)-fabs(value-5))/2);
    *G = d_simulation.max(0.0,(4-fabs(value-2)-fabs(value-4))/2);
    *B = d_simulation.max(0.0,(3-fabs(value-1)-fabs(value-2))/2);
}


void LegendWidget::saturation(float value, float *R, float *G, float *B)
{
   float color_r = 0.5;									//The base color whose saturation we change (green).
   float color_g = 0.5;									//Try different colors!
   float color_b = 0;

   if (value<0.5)										//value in [0,0.5]: modulate the luminance from black to the base-color.
   {
       *R = 2*value*color_r;
       *G = 2*value*color_g;
       *B = 2*value*color_b;
   }
   else													//value in [0.5,1]: modulate the saturation from base-color to white.
   {
       value = 2*(value-0.5);
       *R = (1-value)*color_r + 1*value;
       *G = (1-value)*color_g + 1*value;
       *B = (1-value)*color_b + 1*value;
   }
}


void LegendWidget::heatmap(float value, float* R, float* G, float* B)
{
        if (value<0) value = 0; if (value>1) value = 1;
        //all low start, red high midle, green high end
        *B = 0;//std::max(0.0, -((value - 0.2)*(value - 0.2)) + 1);
        *G = std::max(0.0, -((value-1.4)*(value-1.4)) + 1);
        *R = std::max(0.0, -((value - 0.8)*(value - 0.8)) + 1);
}

void LegendWidget::rgb2hsv(float R, float G, float B, float* S, float* V, float* H)
{
    float M = d_simulation.max(R, d_simulation.max(G, B));
    float m = d_simulation.min(R, d_simulation.min(G, B));
    float d = M - m;
    *V = M; //value = max( r ,g ,b)
    *S = (M>0.00001) ? d / M : 0; //saturation
    if(S == 0) H = 0; //achromatic case , hue=0 by convention
    else //chromatic case
    {
        if(R == M) *H = (G - B) / d;
    else if(G == M) *H = 2 + (B - R) / d;
    else *H = 4 + (R - G) / d;
    *H /= 6;
    if(H<0) *H += 1;
    }
}


void LegendWidget::hsv2rgb(float* R, float* G, float* B, float S, float V, float H)
{
    int hueCase = (int) (H*6);
    float frac = 6*H - hueCase;
    float lx = V*(1 - S);
    float ly = V*(1 - S * frac);
    float lz = V*(1 - S *(1 - frac));
    switch(hueCase)
    {
case 0:
case 6: *R = V; *G = lz; *B = lx; break; // 0<hue<1/6
case 1: *R = ly; *G = V; *B = lx; break; // 1/6<hue<2/6
case 2: *R = lx; *G = V; *B = lz; break; // 2/6<hue<3/6
case 3: *R = lx; *G = ly; *B = V; break; // 3/6<hue/4/6
case 4: *R = lz; *G = lx; *B = V; break; // 4/6<hue<5/6
case 5: *R = V; *G = lx; *B = ly; break; // 5/6<hue<1
    }
}



void LegendWidget::set_colormap(float vy)
{
    float R,G,B;
    float H, S, V;

    if (scalar_col==COLOR_BLACKWHITE)
        R = G = B = vy;
    else if (scalar_col==COLOR_RAINBOW)
        rainbow(vy,&R,&G,&B);
    else if (scalar_col==COLOR_YELLOW)
        heatmap(vy,&R,&G,&B);
    else if (scalar_col==COLOR_BANDS)
    {
        const int NLEVELS = 7;
        vy *= NLEVELS;
        vy = (int)(vy);
        vy/= NLEVELS;
        rainbow(vy,&R,&G,&B);
    }

    rgb2hsv(R, G, B, &S, &V, &H);
    hsv2rgb(&R, &G, &B, (float)(increase_s), (float)(increase_v), H);

    glColor3f(R,G,B);
}
