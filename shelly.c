
#include <math.h>
#include <stdio.h>
#include <fourier.h>
#include <gsl/gsl_rng.h>

#include "shelly.h"

#define M_PI   3.14159265358979323846264338327950288
#define PI M_PI


//////////////////////////////////////////
// flow field ////////////////////////////
const double mu = 1.0;                  // viscosity
const double rho = 1.0;                 // density (not used in practice)
const int isStressSymmetric = TRUE;     // we work with a symmetric stress PiS
double F;                               // magnitude of forcing

int numv, numK, numPi, numf;            // size of vectors for velocity, force etc
double *v, *K, *Pi, *f;                 // data for velocity, force etc


//////////////////////////////////////////
// polymer ///////////////////////////////
int numS = 3;                           // symmetric so only 3 entries are required
double *S = 0, *SRe = 0, *SIm = 0;      // memory is allocated in setupDiffusiveSolver
fftw_plan S_FORWARD, S_BACKWARD;        // Fourier transform plans
double *S_LUT = 0;                      // the LUT (lookup table) for solving ∂t S = ∇²S in Fourier space
const double D = 1e-3;                  // diffusive constant
const double tau = 1.0, G = 0.5;        // polymer relaxation time and modulus
double Wi, beta;                        // dimensionless paramters


//////////////////////////////////////////
// numerical parameters //////////////////
int iterations;                         // how many loops have we made
double dt, t, tmax;                     // timestep, time, max time
double t_kick_2D = 10;                  // when is perturbation added

gsl_rng * rng;     // random number generator

int perturbationAdded = FALSE;
int Nx, Ny, Nz;
const double Lx = 2.0 * PI, Ly = 2.0 * PI, Lz = 1.0;



int main(int argc, char **argv)
{
    processArgs(argc, argv);

    F = Wi / Lx;
    beta = G / (Lx * F * rho);

    printboth ("#S# %d x %d\n", Nx, Ny);
    printboth ("#L# %g %g\n", Lx, Ly);
    printboth ("#T# %g %g\n", dt, tmax);
    printboth ("#D# Wi %g\tβ %g\t Wi.β %g\n", Wi, beta, Wi * beta);
    fprintf (stderr, "FORCE: %f\n", F);

    setupSolvers();

    initialiseSystem();



    if (TRUE)
    {
        const int printFreq = (int) (0.1 / dt);
        for (t = 0.0, iterations = 0; t < tmax; t += dt, ++iterations)
        {
            updateSystemFull();

            if (iterations % printFreq == 0)
            {
                fprintf (stderr, "%f\n", t);

                print1D();
                print2D();

                if (isnan(S(Nx / 4, Ny / 4, XX)))
                {
                    fprintf (stderr, "Nan encountered!\n");
                    break;
                }
            }

            if (t > t_kick_2D && !perturbationAdded)
                addPerturbation2D();
        }
    }

    freeSolvers();

    return 0;
}


void updateSystemFull()
{
    solveConvection(numS, S);

    solveOldroydBLocal();

    /*solveDiffusion(numS, SRe, SIm, &S_FORWARD, &S_BACKWARD, S_LUT);*/

    solveOldroydBStress();

    solveFlowField();
}

void initialiseSystem()
{

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            S(i, j, XX) = S(i, j, YY) = 1.0;
            S(i, j, XY) = 0.0;

            const double xx = Lx * (double) i / (double) Nx;
            const double yy = Ly * (double) j / (double) Ny;

            Pi(i, j, XX) = Pi(i, j, XY) = Pi(i, j, YY) = 0.0;


            f(i, j, X) = +2.0 * sin(xx) * cos(yy);
            f(i, j, Y) = -2.0 * cos(xx) * sin(yy);
        }

    solveFlowField();

}

void addPerturbation2D()
{
    // initialise a Mersenne Twister RNG
    rng = gsl_rng_alloc ( gsl_rng_mt19937 );
    // change this to try a different 2D perturbation
    unsigned long int seed = 13; gsl_rng_set(rng, seed);
    const double delta = 1e-3;

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
            for (int d = 0; d < numS; ++d)
                S(i, j, d) += delta * (gsl_rng_uniform(rng) - 0.5);

    // avoid adding this again
    perturbationAdded = TRUE;
}


void solveOldroydBLocal()
{
    double dS_dt[numS];                 // for each point
    const double tinv = 1.0 / Wi;       // calculate once at start

    for (int i = 0; i < Nx; ++i)
    {
        for (int j = 0; j < Ny; ++j)
        {

            // store pointer to local tensor array mainly for
            // aesthetic / readability reasons
            double *Sij = &S(i, j, 0);
            double *Kij = &K(i, j, 0);

            dS_dt[XX]  = Sij[XX] * (2.0*Kij[XX] - tinv);
            dS_dt[XX] += Sij[XY] * Kij[XY]*2.0;
            dS_dt[XX] += tinv;

            dS_dt[XY]  = 0.5 * Sij[XX] * Kij[YX]*2.0;
            dS_dt[XY] += Sij[XY] * (Kij[XX] + Kij[YY] - tinv);
            dS_dt[XY] += 0.5 * Sij[YY] * Kij[XY]*2.0;

            dS_dt[YY]  = Sij[YY] * (2.0*Kij[YY] - tinv);
            dS_dt[YY] += Sij[XY] * Kij[YX]*2.0;
            dS_dt[YY] += tinv;


            // advance in time
            Sij[XX] += dS_dt[XX] * dt;
            Sij[XY] += dS_dt[XY] * dt;
            Sij[YY] += dS_dt[YY] * dt;

        }       // end of spacial loop (y)
    }       // end of spacial loop (x)
}

void solveOldroydBStress()
{

    for (int i = 0; i < Nx; ++i)
        for (int j = 0; j < Ny; ++j)
        {
            // NOTE: that stress IS renormalised here (for a != 1)
            Pi(i, j, XX) = beta * (S(i, j, XX) - 1.0);
            Pi(i, j, YY) = beta * (S(i, j, YY) - 1.0);
            Pi(i, j, XY) = beta * S(i, j, XY);             // = Pi(i, j, YX) as C symmetric

        }   // end of spacial loop
}



void print1D()
{
    for (int i = 0; i < Nx; ++i)
        printf ("%f %f\n", Lx * i / (double) Nx, v(i, 0, X));

    /*printf ("%f %f\n", j / (double) Ny, S(Nx / 2, j, XX));*/
    printf ("\n");
}

void print2D()
{
    if (Nz != 1) return;                // this is designed for 2D plots only

    const int xresolution = max(Nx / (256 * Lx), 1);
    const int yresolution = max(Ny / (256 * Ly), 1);

    // first determine the largest vector magnitude
    double vmagMax = 0.0;
    for (int i = 0; i < Nx; i += xresolution)
        for (int j = 0; j < Ny; j += yresolution)
        {
            // find magnitudes of velocity vectors
            const double vmag = sqrt(v(i, j, X)*v(i, j, X) + v(i, j, Y)*v(i, j, Y));
            if (vmag > vmagMax) vmagMax = vmag;
        }


    FILE *outputFile;
    char outputName[128];

    // files have different suffix for each time
    // NOTE: here the files are output in the same directory as the executable
    // so make sure you copy the executable somewhere sensible

    sprintf(outputName, "st-%010.3f", t);
    outputFile = fopen(outputName, "w+");

    // print header (at top of each file)
    fputs ("#h# #x \t y \t Qxx \t Qxy \t Cxx \t Cxy \t gdot \t ωz \t vmag \t nxny2\n", outputFile);

    for (int i = 0; i < Nx; i += xresolution)
    {
        for (int j = 0; j < Ny; j += yresolution)
        {
            fprintf (outputFile, "%11.9g \t %11.9g \t ", i * Lx / Nx, j * Ly / Ny);


            fprintf (outputFile, "%11.9g \t %11.9g %11.9g\t ", S(i, j, XX), S(i, j, XY), S(i, j, YY));
            fprintf (outputFile, "%11.9g \t ", S(i, j, XX) + S(i, j, YY));


            // print velocity magnitudes and stream function
            const double vmag = sqrt(v(i, j, X)*v(i, j, X) + v(i, j, Y)*v(i, j, Y));
            fprintf (outputFile, "%11.9g \t ", vmag);

            // print velocities
            const double du = v(i, j, X) / vmagMax;
            const double dv = v(i, j, Y) / vmagMax;
            fprintf (outputFile, "%f\t%f\t", du, dv);

            // vorticity (Shelly defines vorticity with a factor of -1
            const double vorticity = -(K(i, j, XY) - K(i, j, YX));
            fprintf (outputFile, "%11.9g \t", vorticity);

            const double fmag = sqrt( f(i, j, X)*f(i, j, X) + f(i, j, Y)*f(i, j, Y));
            fprintf (outputFile, "%f\t", fmag);
            fprintf (outputFile, "%f\t%f\t", f(i, j, X) / fmag, f(i, j, Y) / fmag);

            fprintf (outputFile, "\n");
        }

        fprintf (outputFile, "\n");
    }

    fclose(outputFile);

}



void setupSolvers()
{
    // Set up general system / numerical parameters
    setupSystem(Nx, Ny, Nz, Lx, Ly, Lz, dt);

    // Set up solver for Navier-Stokes (using Oseen). This allocates memory for v, K etc
    // and returns the size of these tensors.
    const double eta = 1.0;
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


void processArgs(int argc, char **argv)
{
    const int argsNeeded = 4;
    const int argsGiven = argc - 1;
    if (argsNeeded != argsGiven)
    {
        fprintf (stderr, "Need %d args, given %d:"
                "\n*N *dt *tmax *Wi\n", argsNeeded, argsGiven);
        exit(-1);
    }
    else
    {
        int N;
        sscanf(*(++argv), "%d", &N);
        Nx = Ny = N; Nz = 1;
        sscanf(*(++argv), "%lf", &dt);
        sscanf(*(++argv), "%lf", &tmax);
        sscanf(*(++argv), "%lf", &Wi);
    }
}
