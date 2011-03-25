
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
   float outside_temp = 200;
   float initial_temp = 20;

   // Microwave effect (W)
   float microwave_effect = 750;

   // dimensions of bacon (m)
   float Lx = 0.15;
   float Ly = 0.15;
   float Lz = 0.15;

   // Fat / Meat partition, Fat at y > DLy
   float DLy = 0.5;

   // Set length of cooking (s)
   float Lt = 120.0;

   // Snapshots per second (25 = realtime video)
   float snapshots_ps = 1;

   // Set size of mesh
   int nx = 7;
   int ny = 7;
   int nz = 7;
   int nt = 1200000;

   // -------------------------------------------
   // Calculate delta values
   // -------------------------------------------

   float dx = Lx / (float) (nx+2);
   float dy = Ly / (float) (ny+2);
   float dz = Lz / (float) (nz+2);
   float dt = Lt / (float) (nt+2);

   // Fat starts at Dny =
   int dny = DLy / Ly * ny;

   // Iterations between snapshots
   int snapshots_pi = snapshots_ps / dt;

   // -------------------------------------------
   // Check courant-Friedrichs-Lewy conditions
   // -------------------------------------------

   float cfl = 1.0/6.0;

   // Check Courant-Friedrichs-Lewy (CFL) condition

   if ( dt/dx/dx > cfl || dt/dy/dy > cfl || dt/dz/dz > cfl ) {
      std::cout << "Error: Courant-Friedrichs-Lewy (CFL) condition broken\n";
      std::cout << "The calculation was stopped because inaccurate and oscillating solutions may occur. \n\n";

      std::cout << "dt/dx^2: " 
      << dt/dx/dx << " and dt/dy^2: " 
      << dt/dy/dy << " and dt/dz^2: " 
      << dt/dz/dz << " should be less than: " 
      << cfl << "\n";

      return 1;
   }

   // -------------------------------------------
   // Initialize system
   // -------------------------------------------

   int n = nx*ny*nz;

   float* temperatures = new float[n];
   float* alphas = new float[n];
   float* betas = new float[n];
   float* boundarys = new float[n];
   float* diagonals = new float[n];
   float* microfield = new float[n];
   float* b = new float[n];

   phys_sys::Init(nx, ny, nz, nt, dny,
                  dx, dy, dz, dt,
                  temperatures, alphas, betas,
                  boundarys, diagonals, microfield);

   // Set initial temperature
   phys_sys::InitializeTemperature( initial_temp );

   // Calculate boundary conditions
   phys_sys::CalculateBoundaryConditions( outside_temp );

   // calculate microwave field
   phys_sys::CalculateWaveField( microwave_effect );

   // Initialize conjugate gradient solver
   cg_solver solver(b, temperatures, n, phys_sys::MultiplyMatrixVector);

   // Print cross-section of initial state to file
   DumpState(temperatures, nx, ny, nz, 0);

   // Main loop
   for ( int i = 1; i < nt; i++ ) {
      // Update alpha values
      phys_sys::UpdateAlphaBetaValues( );

      // Calculate diagonal
      phys_sys::CalculateDiagonal( );

      // b <= Bx
      phys_sys::SetRightSigns();
      phys_sys::MultiplyMatrixVector(b, temperatures);

      // b += 2*boundarys
      for ( int j = 0; j < n; j++ ) {
         b[j] += 2*boundarys[j]; // + 1000 * microfield[j] * betas[j] * microwave_effect * dt ;
      }

      // solve Ax = b for x
      phys_sys::SetLeftSigns();
      solver.Solve();

      // Print cross-section to file
      if( i%snapshots_pi == 0)
         DumpState(temperatures, nx, ny, nz, i/snapshots_pi);
   }

   // Delete allocated space
   delete[] boundarys;
   delete[] temperatures;
   delete[] alphas;
   delete[] betas;
   delete[] microfield;
   delete[] b;

   return 0;
}

// Print cross-section
void
DumpState(float *t, int nx, int ny, int nz, int i) {

   std::cout << "snapshot at iteration " << i << "\n";

   char str[10];
   sprintf (str, "%d", i);
   std::ofstream outfile;
   outfile.open( (std::string("data/") + std::string(str) + std::string(".dat")).c_str() , std::ios::out);

   int z=nz/2;
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {
            outfile << x << " " << y << " " << t[z*nz*ny + y*nx + x] << "\n";
         }
         outfile << "\n";
   }

   outfile.close();
}
