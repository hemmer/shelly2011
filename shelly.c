
#include <math.h>
#include <stdio.h>
#include <fourier.h>

#include "shelly.h"

#define M_PI   3.14159265358979323846264338327950288
#define PI M_PI


// flow field
const double eta = 1.0;
const int isStressSymmetric = TRUE;

int numv, numK, numPi, numf;
double *v, *K, *Pi, *f;

//////////////////////////////////////////
// polymer ///////////////////////////////
int numS = 3;                           // symmetric so only 3 entries are required
double *S = 0, *SRe = 0, *SIm = 0;      // memory is allocated in setupDiffusiveSolver
fftw_plan S_FORWARD, S_BACKWARD;        // Fourier transform plans
double *S_LUT = 0;                      // the LUT (lookup table) for solving ∂t S = ∇²S in Fourier space
const double D = 1e-3;                  // diffusive constant

// numerical
double dt = 0.01;
/*int Nx = 16, Ny = 16, Nz = 1;*/
int Nx = 64, Ny = 64, Nz = 1;
/*int Nx = 128, Ny = 128, Nz = 1;*/



int main(int argc, char **argv)
{
    SUPPRESS_UNUSED_WARNING(argc);
    SUPPRESS_UNUSED_WARNING(argv);


    setupSolvers();

    initialiseSystem();


    solveFlow();
    printFlow();


    /*const double tmax = 100; double t; int iterations;*/
    /*for (t = 0.0, iterations = 0; t < tmax; t += dt, ++iterations)*/
    /*{*/
        /*solveDiffusion(numS, SRe, SIm, &S_FORWARD, &S_BACKWARD, S_LUT);*/

        /*if (iterations % 100 == 0)*/
        /*{*/
            /*fprintf (stderr, "%f\n", t);*/
            /*printSlice();*/
        /*}*/
    /*}*/






    freeSolvers();

    return 0;
}

void initialiseSystem()
{

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            S(i, j, XX) = S(i, j, YY) = 1.0;
            S(i, j, XY) = 0.0;

            const double xx = (double) i / (double) Nx;
            const double yy = (double) j / (double) Ny;

            Pi(i, j, XX) = Pi(i, j, XY) = Pi(i, j, YY) = 0.0;

            f(i, j, X) = -2.0 * sin(2.0 * PI * xx) * cos(2.0 * PI * yy);
            f(i, j, Y) = +2.0 * cos(2.0 * PI * xx) * sin(2.0 * PI * yy);
        }
}

void printSlice()
{
    for (int j = 0; j < Ny; ++j)
        printf ("%f %f\n", j / (double) Ny, S(Nx / 2, j, XX));
    printf ("\n");
}


void setupSolvers()
{
    // Set up general system / numerical parameters
    setupSystem(Nx, Ny, Nz, dt);

    // Set up solver for Navier-Stokes (using Oseen). This allocates memory for v, K etc
    // and returns the size of these tensors.
    setupFlowSolver(eta, isStressSymmetric, &v, &numv, &K, &numK, &Pi, &numPi, &f, &numf);

    // Setup our polymer order parameter S,
    setupDiffusiveSolver(numS, &S, &SRe, &SIm, D, &S_FORWARD, &S_BACKWARD, &S_LUT);
}

void freeSolvers()
{

    freeFlowSolver(&v, &K, &Pi, &f);
    freeDiffusiveSolver(&S, &SRe, &SIm, &S_FORWARD, &S_BACKWARD, &S_LUT);

    // This *must* be performed last, as it performs the 
    // final FFTW cleanup (may get memory leaks otherwise)
    freeSystem();
}
