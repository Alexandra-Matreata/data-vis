#include "smokewidget.h"
#include <GL/glut.h>
#include <QTimer>
#include <QObject>
#include <iostream>
#include <math.h>

SmokeWidget::SmokeWidget(QWidget *parent)
: QOpenGLWidget(parent)
{
    DIM = 50;

    color_dir = 0;
    vec_scale = 1000;
    draw_smoke = 0;
    draw_vecs = 1;
    scalar_col = 1;
    frozen = 0;

    increase_s = 0.8;
    increase_v = 0.8;
    dataset = 'r';
    scalar_dataset = 'r';
    vector_dataset = 'v';
    glyphs = 'c';

    select_points = 1;

    QTimer *timer = new QTimer;
    timer->start(1);
    QObject::connect(timer,SIGNAL(timeout()),this,SLOT(do_one_simulation_step()));

}

void SmokeWidget::initializeGL()
{
    initializeOpenGLFunctions();

    d_simulation.init_simulation(DIM);

}

void SmokeWidget::paintGL()
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
    //glutSwapBuffers();

}

void SmokeWidget::resizeGL(int w, int h)
{
    glViewport(0.0f, 0.0f, (GLfloat)w, (GLfloat)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, (GLdouble)w, 0.0, (GLdouble)h);
    winWidth = w; winHeight = h;

}

void SmokeWidget::do_one_simulation_step()
{
    if (!d_simulation.get_frozen())
    {
        d_simulation.set_forces(DIM);
        d_simulation.solve(DIM);
        d_simulation.diffuse_matter(DIM);

    }
    update();
}



void SmokeWidget::visualize()
{
    int        i, j, idx;
    fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
    fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

    if (draw_smoke)
    {
        int idx0, idx1, idx2, idx3;
        double px0, py0, px1, py1, px2, py2, px3, py3;
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glBegin(GL_TRIANGLES);
        for (j = 0; j < DIM - 1; j++)            //draw smoke
        {
            for (i = 0; i < DIM - 1; i++)
            {
                px0 = wn + (fftw_real)i * wn;
                py0 = hn + (fftw_real)j * hn;
                idx0 = (j * DIM) + i;


                px1 = wn + (fftw_real)i * wn;
                py1 = hn + (fftw_real)(j + 1) * hn;
                idx1 = ((j + 1) * DIM) + i;


                px2 = wn + (fftw_real)(i + 1) * wn;
                py2 = hn + (fftw_real)(j + 1) * hn;
                idx2 = ((j + 1) * DIM) + (i + 1);


                px3 = wn + (fftw_real)(i + 1) * wn;
                py3 = hn + (fftw_real)j * hn;
                idx3 = (j * DIM) + (i + 1);

                switch (dataset) {
                case 'r':
                    set_colormap(d_simulation.get_rho()[idx0]);    glVertex3f(px0, py0, 0.0);
                    set_colormap(d_simulation.get_rho()[idx1]);    glVertex3f(px1, py1, 0.0);
                    set_colormap(d_simulation.get_rho()[idx2]);    glVertex3f(px2, py2, 0.0);


                    set_colormap(d_simulation.get_rho()[idx0]);    glVertex3f(px0, py0, 0.0);
                    set_colormap(d_simulation.get_rho()[idx2]);    glVertex3f(px2, py2, 0.0);
                    set_colormap(d_simulation.get_rho()[idx3]);    glVertex3f(px3, py3, 0.0);
                    break;
                case 'f':
                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx0]*100,2)+pow(d_simulation.get_fy()[idx0]*100,2)));    glVertex3f(px0, py0, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx1]*100,2)+pow(d_simulation.get_fy()[idx1]*100,2)));    glVertex3f(px1, py1, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx2]*100,2)+pow(d_simulation.get_fy()[idx2]*100,2)));    glVertex3f(px2, py2, 0.0);


                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx0]*100,2)+pow(d_simulation.get_fy()[idx0]*100,2)));    glVertex3f(px0, py0, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx2]*100,2)+pow(d_simulation.get_fy()[idx2]*100,2)));    glVertex3f(px2, py2, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_fx()[idx3]*100,2)+pow(d_simulation.get_fy()[idx3]*100,2)));    glVertex3f(px3, py3, 0.0);
                    break;
                case 'v':
                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx0]*100,2)+pow(d_simulation.get_vy()[idx0]*100,2)));    glVertex3f(px0, py0, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx1]*100,2)+pow(d_simulation.get_vy()[idx1]*100,2)));    glVertex3f(px1, py1, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx2]*100,2)+pow(d_simulation.get_vy()[idx2]*100,2)));    glVertex3f(px2, py2, 0.0);


                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx0]*100,2)+pow(d_simulation.get_vy()[idx0]*100,2)));    glVertex3f(px0, py0, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx2]*100,2)+pow(d_simulation.get_vy()[idx2]*100,2)));    glVertex3f(px2, py2, 0.0);
                    set_colormap(sqrt(pow(d_simulation.get_vx()[idx3]*100,2)+pow(d_simulation.get_vy()[idx3]*100,2)));    glVertex3f(px3, py3, 0.0);
                    break;
                default:
                    break;
                }
            }
        }
        glEnd();
         upd(d_simulation.round(d_simulation.get_rho_max()), d_simulation.round(d_simulation.get_rho_min()));
    }

    if (draw_vecs)
        {//LINES
            /*glBegin(GL_LINES);				//draw velocities
            for (i = 0; i < DIM; i++)
                for (j = 0; j < DIM; j++)
                {
                    idx = (j * DIM) + i;
                    direction_to_color(d_simulation.get_vx()[idx],d_simulation.get_vy()[idx],color_dir);
                    glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
                    glVertex2f((wn + (fftw_real)i * wn) + vec_scale * d_simulation.get_vx()[idx], (hn + (fftw_real)j * hn) + vec_scale * d_simulation.get_vy()[idx]);
                }
            glEnd();

            glColor3f(1.0, 1.0, 1.0);*/
           // glClear(GL_COLOR_BUFFER_BIT);

        //TRIANGLES
           /* glBegin(GL_TRIANGLES);
            for (i = 0; i < DIM; i++)
                for (j = 0; j < DIM; j++)
                {
                    idx = (j * DIM) + i;
            //glColor3f(1.0, 1.0, 1.0);
                    glVertex3f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn, 0.0);
                    //glColor3f(1.0, 1.0, 0.0);
                    glVertex3f((wn + (fftw_real)i * wn) + vec_scale * d_simulation.get_vx()[idx], (hn + (fftw_real)j * hn) + vec_scale * d_simulation.get_vy()[idx], 0.0);
                    //glColor3f(1.0, 0, 0);
                    glVertex3f(wn + 5+ (fftw_real)i * wn, hn+5 + (fftw_real)j * hn, 0.0);
                }

            glEnd();*/


        //CONES
            int cell_width = 10;
            int cell_height = 10;
            float radius;



            glClear(GL_COLOR_BUFFER_BIT);
            glColor3f(1.0,1.0,0.0);

           for (i = 0; i < DIM; i++)
               for (j = 0; j < DIM; j++)
                {
                   idx = (j * DIM) + i;
                   //compute radius according to scalar dataset:
                   switch (scalar_dataset) {
                       case 'r':
                           radius = d_simulation.get_rho()[idx]; // calculate radius
                       break;
                       case 'v':
                           radius = (sqrt(pow(d_simulation.get_vx()[idx],2) + pow(d_simulation.get_vy()[idx],2)))*200; // calculate radius
                       break;
                       case 'f':
                           radius = (sqrt(pow(d_simulation.get_fx()[idx],2) + pow(d_simulation.get_fy()[idx],2)))*200; // calculate radius
                       break;
                   }
                   float direction;
                   switch (vector_dataset) {
                       case 'v':
                            direction_to_color(d_simulation.get_vx()[idx],-d_simulation.get_vy()[idx],color_dir);
                            direction = atan2 (d_simulation.get_vx()[idx], -d_simulation.get_vy()[idx]) * (180 / M_PI);
                       break;
                       case 'f':
                            direction_to_color(d_simulation.get_fx()[idx],-d_simulation.get_fy()[idx],color_dir);
                            direction = atan2 (d_simulation.get_fx()[idx], -d_simulation.get_fy()[idx]) * (180 / M_PI);
                       break;
                   }

                  //float length = sqrt(pow(d_simulation.get_vx()[idx],2) + pow(d_simulation.get_vy()[idx],2));

                    glPushMatrix();
                    glTranslatef(wn + (fftw_real)i * wn,hn + (fftw_real)j * hn,0);
                    glRotatef(direction,0,0,1);
                    glBegin(GL_TRIANGLE_FAN);
                         glVertex2f(cell_width, cell_height); // draw cone point/tip
                         for (int angle = 1; angle <= 360; angle++) {
                              glColor4f(1,0,0,0.5/color_dir); // colors(R, G, B, alpha)
                              glVertex2f(sin(angle) * radius, cos(angle) * radius); // draw cone base (circle)
                         }
                    glEnd();
                   // float radius = d_simulation.get_rho()[idx]; // calculate radius
                    //int slices = 30; int stacks = 10;
                    //glutSolidCone(radius, cell_height, slices, stacks);
                    glPopMatrix();
                 }


    }
}


void SmokeWidget::display()
{

}


void SmokeWidget::reshape(int w, int h)
{

}



void SmokeWidget::direction_to_color(float x, float y, int method)
{
    float r,g,b,f;
    if (method)
    {
        f = atan2(y,x) / 3.1415927 + 1;
        r = f;
        if(r > 1)
            r = 2 - r;
        g = f + .66667;
        if(g > 2)
            g -= 2;
        if(g > 1)
            g = 2 - g;
        b = f + 2 * .66667;
        if(b > 2)
            b -= 2;
        if(b > 1)
            b = 2 - b;
    }
    else
    {
        r = g = b = 1;
    }
    glColor3f(r,g,b);
}


/*void SmokeWidget::length_to_color(fftw_real *x, fftw_real *y, int method)
{
    int        i, j;
    float max, min = 0;
    fftw_real  wn = (fftw_real)winWidth / (fftw_real)(DIM + 1);   // Grid cell width
    fftw_real  hn = (fftw_real)winHeight / (fftw_real)(DIM + 1);  // Grid cell heigh

    for (j = 0; j < DIM - 1; j++)
    {
        for (i = 0; i < DIM - 1; i++)
        {
            float tmp = sqrt(x[j]*x[j]+y[i]*y[i]);
            if (tmp > max)
                max = tmp;
        }
    }

    float range = max - min;

    for (j = 0; j < DIM - 1; j++)
    {
        for (i = 0; i < DIM - 1; i++)
        {
            float tmp = sqrt(x[j]*x[j]+y[i]*y[i]);
            set_colormap((tmp-min)/range);
            glVertex2f(wn + (fftw_real)i * wn, hn + (fftw_real)j * hn);
            glVertex2f((wn + (fftw_real)i * wn) + vec_scale * x[j], (hn + (fftw_real)j * hn) + vec_scale * y[i]);
        }
    }

}
*/



void SmokeWidget::rainbow(float value, float *R, float *G, float *B)
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


void SmokeWidget::heatmap(float value, float* R, float* G, float* B)
{
    if (value<0) value = 0; if (value>1) value = 1;
    //all low start, red high midle, green high end
    *B = 0;//std::max(0.0, -((value - 0.2)*(value - 0.2)) + 1);
    *G = std::max(0.0, -((value-1.4)*(value-1.4)) + 1);
    *R = std::max(0.0, -((value - 0.8)*(value - 0.8)) + 1);
}

void SmokeWidget::rgb2hsv(float R, float G, float B, float* S, float* V, float* H)
{
    float M = std::max(R, std::max(G, B));
    float m = std::min(R, std::min(G, B));
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


void SmokeWidget::hsv2rgb(float* R, float* G, float* B, float S, float V, float H)
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


void SmokeWidget::saturation(float value, float *R, float *G, float *B)
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


void SmokeWidget::set_colormap(float vy)
{
    float R,G,B;
    float H,S,V;

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

void SmokeWidget::updateLabel(){

}

void SmokeWidget::mousePressEvent(QMouseEvent *event)
{
    lastPos = event->pos();
    if(select_points){
        mouse_x.insert(mouse_x.end(), lastPos.x()*2); //
        mouse_y.insert(mouse_y.end(), winHeight - lastPos.y()*2);//
    }
}

void SmokeWidget::mouseMoveEvent(QMouseEvent *event)
{
    int mx = event->x();// - lastposition gets calculated in drag(), could save a step by using lastPos.x/y but leaving it like this is safer
    int my = event->y();
    //simulation.drag(mx,my, DIM, windowWidth, windowHeight);  // Works for Freerk when using external display
    d_simulation.drag(mx,my,DIM,winHeight,winWidth); // Works for Niek
}
