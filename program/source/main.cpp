
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

#include "headers/phys_sys.h"
#include "headers/phys_consts.h"
#include "headers/cgsolver.h"

using namespace std;

/* Functions declarations */
void DumpHeat(float *x, int nx, int ny, int nz, float dx, float dy, float dz, int i);
void DumpFlow(float *x, int nx, int ny, int nz, float dx, float dy, float dz, int i);

int
main(int argc, char *argv[]) {

   time_t start_time = time(0);

   // -------------------------------------------
   // Set number of threads
   // -------------------------------------------

   omp_set_num_threads(4);

   // -------------------------------------------
   // initialize problem
   // -------------------------------------------

   // Initial temperatures
   float outside_temp = 20;
   float initial_temp = 20;

   // Microwave effect (W)
   float microwave_effect = 750;

   // dimensions of bacon (m)
   float Lx = 0.15;
   float Ly = 0.05;
   float Lz = 0.01;

   // Fat / Meat partition, Fat at y > DLy
   float DLy = 0.03;

   // Set length of cooking (s)
   float Lt = 60.0;

   // Snapshots per second (25 = realtime video)
   float snapshots_ps = 25;

   // Set size of mesh
   int nx = 15;
   int ny = 15;
   int nz = 15;
   int nt = 600000;

   // -------------------------------------------
   // Calculate delta values
   // -------------------------------------------

   float dx = Lx / (float) (nx);
   float dy = Ly / (float) (ny);
   float dz = Lz / (float) (nz);
   float dt = Lt / (float) (nt);

   // Print delta values
   cout << "dt " << dt << " dx " << dx << " dy " << dy << " dz " << dz << "\n";

   // Fat starts at Dny =
   int dny = DLy / Ly * ny;

   // Iterations between snapshots
   int snapshots_pi = 1 / snapshots_ps / dt;

   // -------------------------------------------
   // Check courant-Friedrichs-Lewy conditions
   // -------------------------------------------

   float cfl = 0.5;
   // Check Courant-Friedrichs-Lewy (CFL) condition

   float greatest_alpha = 0.04;
   float s = greatest_alpha * ( dt/dx/dx + dt/dy/dy + dt/dz/dz );

   if ( s > cfl ) {
      std::cout << "\nError: Courant-Friedrichs-Lewy (CFL) condition broken\n"
                << "The calculation was stopped because " 
                << "inaccurate and oscillating solutions may occur. \n\n"
                << "a * ( dt/dx^2 + dt/dy^2 + dt/dz^2 ) " << s
                << " should be less than: " << cfl << "\n"
                << "Do you want to continue? [y/n]\n";
      
      std::string answer;
      std::cin >> answer;

      if ( answer.c_str()[0] != 'y' ) {
         return 1;
      }
   }

   // -------------------------------------------
   // Initialize system
   // -------------------------------------------

   int n = nx*ny*nz;

   // Allocate space
   float* temperatures = new float[n];
   float* alphas = new float[n];
   float* betas = new float[n];
   float* boundarys = new float[n];
   float* diagonals = new float[n];
   float* microfield = new float[n];
   float* b = new float[n];
   float* flow = new float[n];
   float* temp_flow = new float[n];
   float* epsilon = new float[n];
   float* prev_eps = new float[n];

   // Initialize physical system
   phys_sys::Init(nx, ny, nz, nt, dny,
                  dx, dy, dz, dt,
                  temperatures, alphas, betas,
                  boundarys, diagonals, microfield,
                  flow, temp_flow, epsilon, prev_eps);

   // Set initial temperature
   phys_sys::InitializeTemperature( initial_temp );

   // Calculate heat boundary conditions
   phys_sys::CalculateBoundaryConditions( outside_temp );

   // calculate microwave field
   phys_sys::CalculateWaveField( microwave_effect );

   // Initialize conjugate gradient solver
   cg_solver solver(b, temperatures, n, phys_sys::MultiplyMatrixVector);

   // Print cross-sections of heat and flow state to file
   DumpHeat(temperatures, nx, ny, nz, dx, dy, dz, 0);
   DumpFlow(flow, nx, ny, nz, dx, dy, dz, 0);

   // Main loop
   for ( int i = 1; i < nt; i++ ) {

      // Update values depending on the current temperature
      phys_sys::UpdateAlphaBetaValues( );

      // Calculate diagonal of the matrix for the heat equation
      phys_sys::CalculateDiagonal( );

      // solve Bx for right-hand side of heat equation
      phys_sys::SetRightSigns();
      phys_sys::MultiplyMatrixVector(b, temperatures);

      // add boundary conditions and microwave heat to right-
      // hand side of heat equation

      for ( int j = 0; j < n; j++ ) {
         b[j] += 2 * boundarys[j]
                  + microfield[j] * betas[j] 
                  * microwave_effect * dt / dx / dy / dz;
      }

      // solve heat equation using the conjugate gradient method
      phys_sys::SetLeftSigns();
      solver.Solve();

      // calculate flow with updated heat values
      phys_sys::CalculateFlow();
      
      // Print cross-sections of temperature and flow to files
      if( i%snapshots_pi == 0) {
         DumpHeat(temperatures, nx, ny, nz, dx, dy, dz, i/snapshots_pi);
         DumpFlow(flow, nx, ny, nz, dx, dy, dz, i/snapshots_pi);
      }
   }

   time_t end_time = time(0);

   time_t timed = end_time - start_time;

   cout << timed << "\n";

   // Delete allocated space
   delete[] boundarys;
   delete[] temperatures;
   delete[] alphas;
   delete[] betas;
   delete[] microfield;
   delete[] b;
   delete[] flow;
   delete[] temp_flow;
   delete[] epsilon;
   delete[] prev_eps;
   return 0;
}

// Print cross-section of heat
void
DumpHeat(float *t, int nx, int ny, int nz, float dx, float dy, float dz, int i) {

   std::cout << "snapshot at iteration " << i << "\n";

   stringstream ss;
   ss << i;
   string str;
   ss >> str;


   int str_length = str.length();
   for (int j = 0; j < 7 - str_length; j++)
      str = "0" + str;

   //sprintf (str, "%d", i);
   std::ofstream outfile;
   outfile.open( (string("data/") + string(str) + string(".dat")).c_str(), 
                 std::ios::out);

   int z=nz/2;
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {
            outfile << x*dx << " " << y*dy << " " << t[z*nz*ny + y*nx + x] << "\n";
         }
         outfile << "\n";
   }

   outfile.close();
}

// Print cross-section of flow field
void
DumpFlow(float *t, int nx, int ny, int nz, float dx, float dy, float dz, int i) {
   
   stringstream ss;
   ss << i;
   string str;
   ss >> str;


   int str_length = str.length();
   for (int j = 0; j < 7 - str_length; j++)
      str = "0" + str;
   
   std::ofstream outfile;
   outfile.open( (std::string("data/flow") + std::string(str) + std::string(".flo")).c_str() , std::ios::out);

   int z=nz/2;
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {
            outfile << x*dx << " " << y*dy << " " << t[z*nz*ny + y*nx + x] << "\n";
         }
         outfile << "\n";
   }

   outfile.close();
}
