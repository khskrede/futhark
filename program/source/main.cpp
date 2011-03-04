#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "headers/cgsolver.h"
#include "headers/diagmatrix.h"

using namespace std;

/* Functions */

void CalculateWaveField(float *f, int nx, int ny, int nz, float dx, float dy, float dz, float dt, float effect);
void CalculateBoundaryConditions(float *b, int nx, int ny, int nz, float hx, float hy, float hz, float dt);
void DumpState(float *x, int nx, int ny, int nz, int i);

int
main(int argc, char *argv[]) {

   // Set number of processors
   omp_set_num_threads(4);

   // Initialize variables
   float Lx, Ly, Lz, Lt;
   int nx, ny, nz, nt;
   float dx, dy, dz, dt;
   float r, diagA, diagB, hx, hy, hz;
   float initial_temp, outside_temp, alpha, microwave_effect, cfl;

   // initialize problem
   outside_temp = 20;
   initial_temp = 20;
   microwave_effect = 300;
   
   Lx = 100.0; Ly = 100.0; Lz = 100.0; Lt = 1000.0;
   nx = 20; ny = 20; nz = 20; nt = 10000;

   dx = Lx / (float) (nx+1);
   dy = Ly / (float) (ny+1);
   dz = Lz / (float) (nz+1);
   dt = Lt / (float) (nt+1);

   cfl = 1.0/12.0;

   // Calculate diagonals of the A and B matrices
   alpha = 0.5;

   hx = alpha * dt / 2 / dx / dx;
   hy = alpha * dt / 2 / dy / dy;
   hz = alpha * dt / 2 / dz / dz;

   r = 2 * ( hx + hy + hz );

   diagA = 1 + r;
   diagB = 1 - r;

   // Check Courant-Friedrichs-Lewy (CFL) condition
   if ( dt/dx/dx > cfl || dt/dy/dy > cfl || dt/dz/dz > cfl  ) {
      std::cout << "Error: Courant-Friedrichs-Lewy (CFL) condition broken\n";
      std::cout << "The calculation was stopped because accuracy can be inaccurate and oscillations can occur. \n\n";
      return 1;
   }

   std::cout << "hx: " << hx << " hy: " << hy << " hz: " << hz << "\n";

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
   CalculateWaveField( f, nx, ny, nz, dx, dy, dz, dt, microwave_effect );

   // Calculate boundary conditions
   CalculateBoundaryConditions( b0, nx, ny, nz, hx, hy, hz, outside_temp );
   
   // Initialize A matrix
   int offsetsA[7] = { -nx*ny, -nx, -1, 0, 1, nx, nx*ny };
   float valuesA[7] = { -hz, -hy, -hx, diagA, -hx, -hy, -hz };
   DIAG_MATRIX *A = new DIAG_MATRIX(valuesA, offsetsA, 7, nx, ny, nz);

   // Initialize B matrix
   int offsetsB[7] = { -nx*ny, -nx, -1, 0, 1, nx, nx*ny };
   float valuesB[7] = { hz, hy, hx, diagB, hx, hy, hz };
   DIAG_MATRIX *B = new DIAG_MATRIX(valuesB, offsetsB, 7, nx, ny, nz);

   // Print cross-section of initial state to file
   DumpState(x, nx, ny, nz, 0);
   
   // Initialize conjugate gradient solver
   CG_SOLVER solver(b, x, nx, ny, nz, A);

   // Main loop
   for ( int i = 1; i < nt; i++ ) {
   
      // b <= Bx + bi + b(i+1) + f
      (*B).MultiplyVector(b, x, big_n);
      for( int j = 0; j < big_n; j++ ) {
         b[j] += 2*b0[j] + f[j];
      }

      // solve Ax = b for x
      solver.Solve();

      // Print cross-section to file
      if( i % 10 == 0)
         DumpState(x, nx, ny, nz, i/10);
   }

   delete A;
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
CalculateWaveField(float *f, int nx, int ny, int nz, float dx, float dy, float dz, float dt, float effect) {

   float added_heat = 5 * effect * dt;
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

            f[k*nx*ny + j*nx + i] = - 1.33789
                                    + 4.54712 * r
                                    - 1.25814 * r * r
                                    + 0.127253 * r * r* r
                                    - 0.00469462 * r * r * r * r
                                    + 0.0000385689 * r * r * r * r * r;
            sum += f[k*nx*ny + j*nx + i];
         }
      }
   }

   // Normalize values in f
   for (int k = 0; k<nz; k++ ) {
      for (int j = 0; j<ny; j++ ) {
         for (int i = 0; i<nx; i++ ) {
            f[k*nx*ny + j*nx + i] = f[k*nx*ny + j*nx + i] / sum * added_heat;
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
               b[ k*nx*ny + j*nx + i ] += hz * temp;
            }
            if ( j == 0 || j == ny-1 ) {
               b[ k*nx*ny + j*nx + i ] += hz * temp;
            }
            if ( k == 0 || k == nz-1 ) {
               b[ k*nx*ny + j*nx + i ] += hx * temp;
            }
         }
      }
   }

}
