#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "headers/cgsolver.h"

using namespace std;

/* Functions */

void CalculateWaveField(float *f, int nx, int ny, int nz, float dt);

void DumpState(float *x, int nx, int ny, int nz, int i);

int
main(int argc, char *argv[]) {

   // Initialize variables
   float Lx, Ly, Lz, Lt;
   int nx, ny, nz, nt;
   float dx, dy, dz, dt;
   float r, diag, hx, hy, hz;
   float initial_temp, outside_temp, alpha;

   // initialize problem
   outside_temp = 10;
   initial_temp = 20;
   Lx = 10; Ly = 10; Lz = 10; Lt = 10;
   nx = 50; ny = 50; nz = 50; nt = 20;

   dx = Lx / (float)nx;
   dy = Ly / (float)ny;
   dz = Lz / (float)nz;
   dt = Lt / (float)nt;

   // Calculate diagonals of the A matrix
   alpha = 1;
   r = alpha * dt / 2 * ( 1/dx/dx + 1/dy/dy + 1/dz/dz );
   diag = 1 + 2*r;
   hx = -r;
   hy = -r;
   hz = -r;

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
   for ( int k = 0; k < nz; k++ ) {
      for ( int j = 0; j < ny; j++ ) {
         for ( int i = 0; i < nx; i++ ) {
            b0[ k*nx*ny + j*nx + i ] = 0;
            if ( i == 0 || i == nx-1 ) {
               b0[ k*nx*ny + j*nx + i ] += hx;
            }
            if ( j == 0 || j == ny-1 ) {
               b0[ k*nx*ny + j*nx + i ] += hy;
            }
            if ( k == 0 || k == nz-1 ) {
               b0[ k*nx*ny + j*nx + i ] += hz;
            }
//            cout << b0[k*nx*ny + j*nx + i] << " ";
         }
//         cout << "\n";
      }
//      cout << "\n";
   }

   // Initialize conjugate gradient solver
   CG_SOLVER solver(b, x, nx, ny, nz);
   int* offsets = new int[4];
   offsets[0] = 0;
   offsets[1] = 1;
   offsets[2] = nx;
   offsets[3] = nx*ny;
   float* values = new float[4];
   values[0] = diag;
   values[1] = hx;
   values[2] = hy;
   values[3] = hz;
   solver.Diagonalize(values, offsets, 4);

   //solver.PrintA();

   // Print cross-section of initial state to file
   DumpState(x, nx, ny, nz, 0);

   // print h values
   std::cout << "diag: " << diag << " hx: " << hx << " hy: " << hy << " hz: " << hz << "\n";

   for ( int i = 1; i < nt; i++ ) {
   
      // b <= Ax + b0 + f
      solver.MultiplyDiagAVector(b, x, big_n);
      for( int j = 0; j < big_n; j++ ) {
         b[j] = -b[j] + b0[j] + f[j];
      }

      // solve Ax = b for x
      solver.Solve();

      // Print cross-section to file
      if( i % 1 == 0)
         DumpState(x, nx, ny, nz, i);
   }

   delete[] offsets;
   delete[] values;
   
   delete[] b0;
   delete[] f;
   delete[] b;
   delete[] x;

   return 0;
}


void
DumpState(float *x, int nx, int ny, int nz, int i) {

   cout << "snapshot at iteration " << i << "\n";

   char str[10];
   sprintf (str, "%d", i);
   ofstream outfile;
   outfile.open( (string("data/") + string(str) + string(".dat")).c_str() , ios::out);

   // Print cross-section
   int k=nz/2;
//   for (int k = 0; k<nz; k++ ) {
      for (int j = 0; j<ny; j++ ) {
         for (int i = 0; i<nx; i++ ) {
            outfile << j << " " << i << " " << x[k*nz*ny + j*nx + i] << "\n";
         }
         outfile << "\n";
//      }
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
