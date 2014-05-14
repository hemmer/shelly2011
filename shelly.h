
#ifndef SHELLY_H_
#define SHELLY_H_


#define S(i, j, d) S[(d) + numS*((j) + Ny * (i))]

void updateSystemFull();

void solveOldroydBLocal();
void solveOldroydBStress();

void initialiseSystem();
void addPerturbation2D();

void print1D();
void print2D();

void setupSolvers();
void freeSolvers();
void processArgs(int argc, char **argv);


#endif
