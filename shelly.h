
#ifndef SHELLY_H_
#define SHELLY_H_

void initialiseSystem();
void setupSolvers();
void freeSolvers();
void printSlice();

#define S(i, j, d) S[(d) + numS*((j) + Ny * (i))]

#endif
