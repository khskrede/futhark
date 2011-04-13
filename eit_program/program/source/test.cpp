
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "headers/phys_sys.h"
#include "headers/phys_consts.h"
#include "headers/cgsolver.h"

using namespace std;

/* Functions */
void DumpState(float *x, int nx, int ny, int nz, int i);
void MultiplyMatrixVector(float *p, float* v);

int
main(int argc, char *argv[]) {

   // -------------------------------------------
   // Set number of threads
   // -------------------------------------------

   omp_set_num_threads(4);

   // -------------------------------------------
   // initialize problem
   // -------------------------------------------

   // Initial temperatures
   float outside_temp = 40;
   float initial_temp = 20;

   // Microwave effect
   float microwave_effect = 300;

   // dimensions of bacon (m)
   float Lx = 1;
   float Ly = 1;
   float Lz = 1;

   // Fat / Meat partition, Fat at y > DLy
   float DLy = 0.5;

   // Set length of cooking (s)
   float Lt = 60.0;

   // Set size of mesh
   int nx = 3;
   int ny = 3;
   int nz = 3;
   int nt = 1000;

   // -------------------------------------------
   // Calculate delta values
   // -------------------------------------------

   float dx = Lx / (float) (nx+1);
   float dy = Ly / (float) (ny+1);
   float dz = Lz / (float) (nz+1);
   float dt = Lt / (float) (nt+1);

   // Fat starts at Dny =
   int dny = DLy / Ly * ny;

   // -------------------------------------------
   // Initialize system
   // -------------------------------------------

   int n = nx*ny*nz;

   float *temperatures = new float[n];
   float *alphas = new float[n];
   float *betas = new float[n];
   float *boundarys = new float[n];
   float *diagonals = new float[n];
   float *microfield = new float[n];

   phys_sys::Init(nx, ny, nz, nt, dny,
                  dx, dy, dz, dt,
                  temperatures, alphas, betas,
                  boundarys, diagonals, microfield);

   cout  << endl << "------------------------" << endl;
   cout << "Starting tests:" << endl;
   cout << "------------------------" << endl << endl;

   // ------------------------------
   // Test InitializeTemperature
   // ------------------------------

   phys_sys::InitializeTemperature( initial_temp );
   bool test = true;
   for ( int i = 0; i < n; i++ ) {
      if ( temperatures[i] != initial_temp ) {
         test = false;
         break;
      }
   }
   if ( !test ) {
      cout << "Failed: InitializeTemperature" << endl << endl;
   }
   else {
      cout << "Passed: InitializeTemperature" << endl << endl;
   }


   // ------------------------------
   // Test CalculateWaveField
   // ------------------------------

   // Only tests if all values in the field is greater or equal to 0
   // and that they sum to approximately 1

   phys_sys::CalculateWaveField( microwave_effect );
   test = true;
   float sum=0;
   for ( int i = 0; i < n; i++ ) {
      if ( microfield[i] < 0 ) {
         test = false;
         cout << microfield[i] << endl;
         break;
      }
      sum += microfield[i];
   }
   if ( !test || sum != 1.0 ) {
      cout << "Failed: CalculateWaveField, sum was: " << sum << endl << endl;
   }
   else {
      cout << "Passed: CalculateWaveField, sum was: " << sum << endl << endl;
   }

   // ------------------------------
   // Test CalculateBoundaryConditions
   // ------------------------------

   cout << "Printing Boundary conditions: " << endl << endl;

   phys_sys::CalculateBoundaryConditions( outside_temp );

   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {
            int i = z*nx*ny + y*nx + x;

            cout << boundarys[i] << " ";
         }
         cout << endl;
      }
      cout << endl;
   }
   cout << "Are these correct?" << endl << endl;

   // ------------------------------
   // Test UpdateAlphaValues
   // ------------------------------

   cout << "Printing Alpha values: " << endl << endl;

   // Change temperatures
   for ( int i = 0; i < n; i++ ) {
      temperatures[i] = 2*i;
   }

   phys_sys::UpdateAlphaBetaValues();

   test = true;

   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {
            int i = z*nx*ny + y*nx + x;

            // 1 = meat, 2=fat

            int substance = 1;

            if ( y > dny ) {
               substance = 2;

               if ( i*2 <= 40 ) {
                  if ( alphas[i] != phys_consts::alpha_values[substance][0] ) {
                     test = false;
                     cout << "solid fat (" << phys_consts::alpha_values[substance][0] << ") ";
                  }
               }
               else if ( i*2 <= 50 ) {
                  if ( alphas[i] != phys_consts::alpha_values[substance][1] ) {
                     test = false;
                     cout << "melting fat (" << phys_consts::alpha_values[substance][1] << ") ";
                  }
               }
               else if ( i*2 <= 2000 ) {
                  if ( alphas[i] != phys_consts::alpha_values[substance][2] ) {
                     test = false;
                     cout << "liquid fat (" << phys_consts::alpha_values[substance][2] << ") ";
                  }
               }
            }
            else {
               substance = 1;

               if ( alphas[i] != phys_consts::alpha_values[substance][0] ) {
                  test = false;
                  cout << "meat (" << phys_consts::alpha_values[substance][0] << ") ";
               }
            }

            cout << alphas[i] << " ";
         }
         cout << endl;
      }
      cout << endl;
   }

   if ( !test ) {
      cout << "Failed: UpdateAlphaBetaValues automatic test " << endl;
   }
   else {
      cout << "Passed: UpdateAlphaBetaValues automatic test " << endl;
   }
   cout << "Do the values look ok?: " << endl << endl;


   // ------------------------------
   // Test CalculateDiagonal
   // ------------------------------

   phys_sys::InitializeTemperature( initial_temp );
   phys_sys::UpdateAlphaBetaValues( );
   phys_sys::CalculateDiagonal( );

   cout << "Printing Diagonal values: " << endl << endl;

   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {
            int i = z*nx*ny + y*nx + x;

            cout << diagonals[i] << " ";
         }
         cout << endl;
      }
      cout << endl;
   }
   cout << "Are these correct?" << endl << endl;

   // ------------------------------
   // Test conjugate gradient solver
   // ------------------------------

   cout << "Testing conjugate gradient solver: " << endl << endl;

   float* b = new float[n];
   for ( int i = 0; i < n; i++ ) {
      b[i] = 30;
   }
   cg_solver solver(b, temperatures, n, MultiplyMatrixVector);
   solver.Solve();

   cout << "Printing temperatures:" << endl << endl;
   test = true;
   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {
            int i = z*nx*ny + y*nx + x;
            cout << temperatures[i] << " ";
            if ( temperatures[i] != b[i] )
               test = false;
         }
         cout << endl;
      }
      cout << endl;
   }
   if ( test ) {
      cout << "Passed: Conjugate gradient solver automatic test" << endl << endl;
   }
   else {
      cout << "Failed: Conjugate gradient solver automatic test" << endl << endl;
   }
   cout << "Does the results look ok?" << endl << endl;

   // ------------------------------
   // Delete allocated space
   // ------------------------------

   delete[] boundarys;
   delete[] temperatures;
   delete[] alphas;
   delete[] betas;
   delete[] microfield;

   cout << endl;

   return 0;
}

void
MultiplyMatrixVector(float* p, float* v) {

   int nx=3, ny=3, nz=3;
   int n = nx*ny*nz;

   for ( int i = 0; i < n; i++ ) {
      p[i] = v[i];
   }

}
