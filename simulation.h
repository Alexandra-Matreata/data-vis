#ifndef SIMULATION_H
#define SIMULATION_H

#include <rfftw.h>             //the numerical simulation FFTW library
#include <fftw.h>             //the numerical simulation FFTW library
#include <stdio.h>              //for printing the help text
#include <cmath>               //for various math functions
#include <GL/glut.h>            //the GLUT graphics library
#include <algorithm>


class Simulation
{
    int n;
    int frozen;               //toggles on/off the animation
    float v_magnitude_min;
    float f_magnitude_min;
    float v_magnitude_max;
    float f_magnitude_max;

    double dt;				//simulation time step
    float visc;				//fluid viscosity
    fftw_real *vx, *vy;             //(vx,vy)   = velocity field at the current moment
    fftw_real *vm;             // velocity magnitude at the current moment
    fftw_real *vx0, *vy0;           //(vx0,vy0) = velocity field at the previous moment
    fftw_real *fx, *fy;	            //(fx,fy)   = user-controlled simulation forces, steered with the mouse
    fftw_real *rho, *rho0;			//smoke density at the current (rho) and previous (rho0) moment
    rfftwnd_plan plan_rc, plan_cr;  //simulation domain discretization
    fftw_real vx_min, vx_max, vy_min, vy_max;
    fftw_real fx_min, fx_max, fy_min, fy_max;
    fftw_real rho_min, rho_max;

    public:
        Simulation();

        void init_simulation(int n);
        void FFT(int direction,void* vx);
        void solve(int n);
        void diffuse_matter(int n);
        void set_forces(int dim);
        void drag(int mx, int my, int dim, int winHeight,int winWidth);


        int clamp(float x);
        float max(float x, float y);
        float min(float x, float y);
        float round(float f);

        float rho_minimal;
        float rho_maximal;

        // Getter and Setter
        bool get_frozen() const;
        double get_dt() const;				//simulation time step
        float get_visc() const;				//fluid viscosity
        fftw_real* get_rho() const;
        fftw_real* get_fy() const;
        fftw_real* get_fx() const;
        fftw_real* get_rho0() const;
        fftw_real* get_vx() const;
        fftw_real* get_vm() const;
        fftw_real* get_vy() const;
        fftw_real* get_vx0() const;
        fftw_real* get_vy0() const;
        fftw_real get_rho_min() const;
        fftw_real get_rho_max() const;
        float get_v_magnitude_min() const;
        float get_v_magnitude_max() const;
        float get_f_magnitude_min() const;
        float get_f_magnitude_max() const;
        //Mutator functions
        void set_frozen(bool);
        void set_visc(float);
        void set_dt(double);
};

inline int Simulation::clamp(float x)
{
    return ((x)>=0.0?((int)(x)):(-((int)(1-(x)))));
}

inline float Simulation::max(float x, float y)
{
    return x > y ? x : y;
}

inline float Simulation::min(float x, float y)
{
    return x > y ? y : x;
}

inline float Simulation::round(float f)
{
    return floor(f * 5 + 0.5) / 5;
    // return std::round(f * 5) / 5; // C++11
}

inline bool Simulation::get_frozen() const
{
    return frozen;
}

inline double Simulation::get_dt() const
{
    return dt;
}

inline float Simulation::get_visc() const
{
    return visc;
}

inline fftw_real* Simulation::get_fx() const
{
    return fx;
}

inline fftw_real* Simulation::get_fy() const
{
    return fy;
}

inline fftw_real* Simulation::get_rho() const
{
    return rho;
}


inline fftw_real* Simulation::get_rho0() const
{
    return rho0;
}

inline fftw_real* Simulation::get_vx() const
{
    return vx;
}

inline fftw_real* Simulation::get_vy() const
{
    return vy;
}

inline fftw_real* Simulation::get_vm() const
{
    return vm;
}

inline fftw_real Simulation::get_rho_min() const
{
    return rho_min;
}

inline fftw_real Simulation::get_rho_max() const
{
    return rho_max;
}

inline float Simulation::get_f_magnitude_min() const
{
    return f_magnitude_min*100;
}

inline float Simulation::get_v_magnitude_min() const
{
    return v_magnitude_min*100;
}

inline float Simulation::get_f_magnitude_max() const
{
    return f_magnitude_max*100;
}

inline float Simulation::get_v_magnitude_max() const
{
    return v_magnitude_max*100;
}

inline fftw_real* Simulation::get_vx0() const
{
    return vx0;
}

inline fftw_real* Simulation::get_vy0() const
{
    return vy0;
}

inline void Simulation::set_frozen(bool new_frozen)
{
    frozen = new_frozen;
}

inline void Simulation::set_dt(double new_dt)
{
    dt = new_dt;
}
inline void Simulation::set_visc(float new_visc)
{
    visc = new_visc;
}



#endif // SIMULATION_H
