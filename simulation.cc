#include "simulation.h"
#include <iostream>

Simulation::Simulation()
{

}

void Simulation::init_simulation(int n)
{
    int i; size_t dim;

    dt = 0.4;				//simulation time step
    visc = 0.001;
    dim     = n * 2*(n/2+1)*sizeof(fftw_real);        //Allocate data structures
    vx       = (fftw_real*) malloc(dim);
    vy       = (fftw_real*) malloc(dim);
    vx0      = (fftw_real*) malloc(dim);
    vy0      = (fftw_real*) malloc(dim);
    dim     = n * n * sizeof(fftw_real);
    fx      = (fftw_real*) malloc(dim);
    fy      = (fftw_real*) malloc(dim);
    rho     = (fftw_real*) malloc(dim);
    rho0    = (fftw_real*) malloc(dim);
    plan_rc = rfftw2d_create_plan(n, n, FFTW_REAL_TO_COMPLEX, FFTW_IN_PLACE);
    plan_cr = rfftw2d_create_plan(n, n, FFTW_COMPLEX_TO_REAL, FFTW_IN_PLACE);
    rho_max = -100; rho_min = 100;

    for (i = 0; i < n * n; i++)                      //Initialize data structures to 0
    { vx[i] = vy[i] = vx0[i] = vy0[i] = fx[i] = fy[i] = rho[i] = rho0[i] = 0.0f; }
}

void Simulation::FFT(int direction, void *vx)
{
    if(direction==1)
        rfftwnd_one_real_to_complex(plan_rc,(fftw_real*)vx,(fftw_complex*)vx);
    else
        rfftwnd_one_complex_to_real(plan_cr,(fftw_complex*)vx,(fftw_real*)vx);
}


void Simulation::solve(int n)
{
    fftw_real x, y, x0, y0, f, r, U[2], V[2], s, t;
    int i, j, i0, j0, i1, j1;

    for (i=0;i<n*n;i++)
    {
        vx[i] += dt*vx0[i];
        vx0[i] = vx[i];
        vy[i] += dt*vy0[i];
        vy0[i] = vy[i];
    }

    for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
        for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
        {
            x0 = n*(x-dt*vx0[i+n*j])-0.5f;
            y0 = n*(y-dt*vy0[i+n*j])-0.5f;
            i0 = clamp(x0); s = x0-i0;
            i0 = (n+(i0%n))%n;
            i1 = (i0+1)%n;
            j0 = clamp(y0); t = y0-j0;
            j0 = (n+(j0%n))%n;
            j1 = (j0+1)%n;
            vx[i+n*j] = (1-s)*((1-t)*vx0[i0+n*j0]+t*vx0[i0+n*j1])+s*((1-t)*vx0[i1+n*j0]+t*vx0[i1+n*j1]);
            vy[i+n*j] = (1-s)*((1-t)*vy0[i0+n*j0]+t*vy0[i0+n*j1])+s*((1-t)*vy0[i1+n*j0]+t*vy0[i1+n*j1]);
        }

    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
        {
            vx0[i+(n+2)*j] = vx[i+n*j];
            vy0[i+(n+2)*j] = vy[i+n*j];
        }

    FFT(1,vx0);
    FFT(1,vy0);

    for (i=0;i<=n;i+=2)
    {
        x = 0.5f*i;
        for (j=0;j<n;j++)
        {
            y = j<=n/2 ? (fftw_real)j : (fftw_real)j-n;
            r = x*x+y*y;
            if ( r==0.0f ) continue;
            f = (fftw_real)exp(-r*dt*visc);
            U[0] = vx0[i  +(n+2)*j];
            V[0] = vy0[i  +(n+2)*j];
            U[1] = vx0[i+1+(n+2)*j];
            V[1] = vy0[i+1+(n+2)*j];

            vx0[i  +(n+2)*j] = f*((1-x*x/r)*U[0]     -x*y/r *V[0]);
            vx0[i+1+(n+2)*j] = f*((1-x*x/r)*U[1]     -x*y/r *V[1]);
            vy0[i+  (n+2)*j] = f*(  -y*x/r *U[0] + (1-y*y/r)*V[0]);
            vy0[i+1+(n+2)*j] = f*(  -y*x/r *U[1] + (1-y*y/r)*V[1]);
        }
    }

    FFT(-1,vx0);
    FFT(-1,vy0);

    f = 1.0/(n*n);
    for (i=0;i<n;i++)
        for (j=0;j<n;j++)
        {
            vx[i+n*j] = f*vx0[i+(n+2)*j];
            vy[i+n*j] = f*vy0[i+(n+2)*j];
        }
}



void Simulation::diffuse_matter(int n)
{
    fftw_real x, y, x0, y0, s, t;
        int i, j, i0, j0, i1, j1, a;

        for ( x=0.5f/n,i=0 ; i<n ; i++,x+=1.0f/n )
            for ( y=0.5f/n,j=0 ; j<n ; j++,y+=1.0f/n )
            {
                x0 = n*(x-dt*vx[i+n*j])-0.5f;
                y0 = n*(y-dt*vy[i+n*j])-0.5f;
                i0 = clamp(x0);
                s = x0-i0;
                i0 = (n+(i0%n))%n;
                i1 = (i0+1)%n;
                j0 = clamp(y0);
                t = y0-j0;
                j0 = (n+(j0%n))%n;
                j1 = (j0+1)%n;
                rho[i+n*j] = (1-s)*((1-t)*rho0[i0+n*j0]+t*rho0[i0+n*j1])+s*((1-t)*rho0[i1+n*j0]+t*rho0[i1+n*j1]);
            }

        for(a=0;a<i+n*j;a++){
            rho_max = rho[a]>rho_max?rho[a]:rho_max;
            rho_min = rho[a]<rho_min?rho[a]:rho_min;
        }

        //rho_max = max(rho);
}


void Simulation::set_forces(int dim)
{
    int i;
        for (i = 0; i < dim * dim; i++)
        {
            rho0[i]  = 0.995 * rho[i];
            fx[i] *= 0.85;
            fy[i] *= 0.85;
            vx0[i]    = fx[i];
            vy0[i]    = fy[i];
        }
}


void Simulation::drag(int mx, int my, int dim,int winHeight,int winWidth)
{
    int xi,yi,X,Y;
    double  dx, dy, len;
    static int lmx=0,lmy=0;				//remembers last mouse location

    // Compute the array index that corresponds to the cursor location
    xi = (int)clamp((double)(dim + 1) * ((double)mx / (double)winWidth));
    yi = (int)clamp((double)(dim + 1) * ((double)(winHeight - my) / (double)winHeight));

    X = xi;
    Y = yi;

    if (X > (dim - 1))
        X = dim - 1;
    if (Y > (dim - 1))
        Y = dim - 1;
    if (X < 0)
        X = 0;
    if (Y < 0)
        Y = 0;

    // Add force at the cursor location
    my = winHeight - my;
    dx = mx - lmx;
    dy = my - lmy;
    len = sqrt(dx * dx + dy * dy);
    if (len != 0.0)
    {
        dx *= 0.1 / len;
        dy *= 0.1 / len;
    }
    fx[Y * dim + X] += dx;
    fy[Y * dim + X] += dy;
    rho[Y * dim + X] = 10.0f;
    lmx = mx;
    lmy = my;
}
