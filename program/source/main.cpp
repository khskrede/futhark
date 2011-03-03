#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "headers/cgsolver.h"
#include "headers/diagmatrix.h"

using namespace std;

/* Functions */

void CalculateWaveField(float *f, int nx, int ny, int nz, float dt);
void CalculateBoundaryConditions(float *b, int nx, int ny, int nz, float hx, float hy, float hz, float dt);
void DumpState(float *x, int nx, int ny, int nz, int i);

int
main(int argc, char *argv[]) {

   // Initialize variables
   float Lx, Ly, Lz, Lt;
   int nx, ny, nz, nt;
   float dx, dy, dz, dt;
   float r, diagA, diagB, hx, hy, hz;
   float initial_temp, outside_temp, alpha;

   // initialize problem
   outside_temp = 100;
   initial_temp = 10;
   Lx = 1.0; Ly = 1.0; Lz = 1.0; Lt = 1000;
   nx = 3; ny = 3; nz = 3; nt = 100;

   dx = Lx / (float) (nx+1);
   dy = Ly / (float) (ny+1);
   dz = Lz / (float) (nz+1);
   dt = Lt / (float) (nt+1);

   // Calculate diagonals of the A and B matrices
   alpha = 1;
   r = alpha * 2 * ( dt/dx/dx + dt/dy/dy + dt/dz/dz );
   diagA = 1 + r;
   diagB = 1 - r;
   hx = alpha * dt/dx/dx;
   hy = alpha * dt/dy/dy;
   hz = alpha * dt/dz/dz;

   // Initialize vectors
   int big_n = nx*ny*nz;
   float *b0 = new float[big_n];
   float *f = new float[big_n];
   float *b = new float[big_n];
   float *x = new float[big_n];

   // Set initial temperature
   for ( int i = 0; i < big_n; i++ ) {
      x[i] = initial_temp;
   }

   // calculate heat field
   CalculateWaveField( f, nx, ny, nz, dt );

   // Calculate boundary conditions
   CalculateBoundaryConditions( b0, nx, ny, nz, hx, hy, hz, outside_temp );
   
   // Initialize A matrix
   int* offsetsA = new int[7];
   offsetsA[0] = -nx*ny;
   offsetsA[1] = -nx;
   offsetsA[2] = -1;
   offsetsA[3] = 0;
   offsetsA[4] = 1;
   offsetsA[5] = nx;
   offsetsA[6] = nx*ny;
   float* valuesA = new float[7];
   valuesA[0] = hz;
   valuesA[1] = hy;
   valuesA[2] = hx;
   valuesA[3] = diagA;
   valuesA[4] = hx;
   valuesA[5] = hy;
   valuesA[6] = hz;
   DIAG_MATRIX *A = new DIAG_MATRIX(valuesA, offsetsA, 7, nx, ny, nz);

   // Initialize B matrix
   int* offsetsB = new int[7];
   offsetsB[0] = -nx*ny;
   offsetsB[1] = -nx;
   offsetsB[2] = -1;
   offsetsB[3] = 0;
   offsetsB[4] = 1;
   offsetsB[5] = nx;
   offsetsB[6] = nx*ny;
   float* valuesB = new float[7];
   valuesB[0] = -hz;
   valuesB[1] = -hy;
   valuesB[2] = -hx;
   valuesB[3] = diagB;
   valuesB[4] = -hx;
   valuesB[5] = -hy;
   valuesB[6] = -hz;
   DIAG_MATRIX *B = new DIAG_MATRIX(valuesB, offsetsB, 7, nx, ny, nz);

   // Print A matrix
   (*A).Print();

   // Print cross-section of initial state to file
   DumpState(x, nx, ny, nz, 0);
   
   // Initialize conjugate gradient solver
   CG_SOLVER solver(b, x, nx, ny, nz, A);

   for ( int i = 1; i < nt; i++ ) {
   
      // b <= Bx + bi + b(i+1) + f
      (*B).MultiplyVector(b, x, big_n);
      for( int j = 0; j < big_n; j++ ) {
         b[j] += 2*b0[j]; // + f[j];
      }

      // solve Ax = b for x
      solver.Solve();

      // Print cross-section to file
      DumpState(x, nx, ny, nz, i);
   }

   delete[] offsetsA;
   delete[] valuesA;
   delete A;
   delete[] offsetsB;
   delete[] valuesB;
   delete B;
   
   delete[] b0;
   delete[] f;
   delete[] b;
   delete[] x;

   return 0;
}


// Print cross-section
void
DumpState(float *x, int nx, int ny, int nz, int i) {

   cout << "snapshot at iteration " << i << "\n";

   char str[10];
   sprintf (str, "%d", i);
   ofstream outfile;
   outfile.open( (string("data/") + string(str) + string(".dat")).c_str() , ios::out);

   int k=nz/2;
      for (int j = 0; j<ny; j++ ) {
         for (int i = 0; i<nx; i++ ) {
            outfile << j << " " << i << " " << x[k*nz*ny + j*nx + i] << "\n";
         }
         outfile << "\n";
   }

   outfile.close();
}


void
CalculateWaveField(float *f, int nx, int ny, int nz, float dt) {

   float effect = nx*ny*nz*dt*30.1;
   float sum = 0;
   float r, r_max;

   // max radius
   r_max = sqrt( (nx / 2) * (nx / 2)
             + (ny / 2) * (ny / 2)
             + (nz / 2) * (nz / 2) );

   // Calculate microwave field
   for (int k = 0; k<nz; k++ ) {
      for (int j = 0; j<ny; j++ ) {
         for (int i = 0; i<nx; i++ ) {

            r = 16 * sqrt( (nx / 2 - i) * (nx / 2 - i)
                      + (ny / 2 - j) * (ny / 2 - j)
                      + (nz / 2 - k) * (nz / 2 - k) ) / r_max;

            f[k*nz*ny + j*nx + i] = - 1.33789
                                    + 4.54712 * r
                                    - 1.25814 * r * r
                                    + 0.127253 * r * r* r
                                    - 0.00469462 * r * r * r * r
                                    + 0.0000385689 * r * r * r * r * r;
            sum += f[k*nz*ny + j*nx + i];
         }
      }
   }

   // Normalize values in f
   for (int k = 0; k<nz; k++ ) {
      for (int j = 0; j<ny; j++ ) {
         for (int i = 0; i<nx; i++ ) {
            f[k*nz*ny + j*nx + i] = f[k*nz*ny + j*nx + i] / sum * effect;
         }
      }
   }

}

void
CalculateBoundaryConditions(float *b, int nx, int ny, int nz, float hx, float hy, float hz, float temp) {

   for ( int k = 0; k < nz; k++ ) {
      for ( int j = 0; j < ny; j++ ) {
         for ( int i = 0; i < nx; i++ ) {
            b[ k*nx*ny + j*nx + i ] = 0;
            if ( i == 0 || i == nx-1 ) {
               b[ k*nx*ny + j*nx + i ] += hx * temp;
            }
            if ( j == 0 || j == ny-1 ) {
               b[ k*nx*ny + j*nx + i ] += hy * temp;
            }
            if ( k == 0 || k == nz-1 ) {
               b[ k*nx*ny + j*nx + i ] += hz * temp;
            }
         }
      }
   }

}
