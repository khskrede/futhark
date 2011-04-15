
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>

#include <math.h>
#include <time.h>
#include <omp.h>

#include "headers/phys_sys.h"
#include "headers/phys_consts.h"

void
phys_sys::Init(int t_nx, int t_ny, int t_nz, int t_nt, int t_dny,
               float t_dx, float t_dy, float t_dz, float t_dt,
               float* t_temperatures, float* t_alphas, float* t_betas,
               float* t_boundarys, float* t_diagonals, float* t_microfield,
               float* t_flow, float* t_temp_flow, float* t_epsilon, float* t_f1, float* t_f2, float* t_prev_eps ) {

   nx = t_nx;
   ny = t_ny;
   nz = t_nz;
   nt = t_nt;
   dny = t_dny;

   n = nx*ny*nz;

   dx = t_dx;
   dy = t_dy;
   dz = t_dz;
   dt = t_dt;

   temperatures = t_temperatures;
   alphas = t_alphas;
   betas = t_betas;
   boundarys = t_boundarys;
   diagonals = t_diagonals;
   microfield = t_microfield;
   flow = t_flow;
   temp_flow = t_temp_flow;
   epsilon = t_epsilon;
   prev_eps = t_prev_eps;
   f1 = t_f1;
   f2 = t_f2;

   offsets[0] = -nx*ny;
   offsets[1] = -nx;
   offsets[2] = -1;
   offsets[3] = 0;
   offsets[4] = 1;
   offsets[5] = nx;
   offsets[6] = nx*ny;

   deltas[0] = dt / 2 / dz / dz;
   deltas[1] = dt / 2 / dy / dy;
   deltas[2] = dt / 2 / dx / dx;
   deltas[3] = 0;
   deltas[4] = dt / 2 / dx / dx;
   deltas[5] = dt / 2 / dy / dy;
   deltas[6] = dt / 2 / dz / dz;

   phys_consts::CalculateValues();
   srand(time(0));
}

void
phys_sys::InitializeTemperature( float initial_temp ) {

   #pragma omp parallel for
   for ( int i = 0; i < n; i++ ) {
      temperatures[i] = initial_temp;
   }
}

void
phys_sys::CalculateWaveField( float effect ) {

   float sum = 0;
   float r, r_max;

   // max radius
   r_max = sqrt( (nx / 2) * (nx / 2)
             + (ny / 2) * (ny / 2)
             + (nz / 2) * (nz / 2) );

   // Calculate microwave field
   for (int z = 0; z<nz; z++ ) {
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {

            int i = z*nx*ny + y*nx + x;

            r = 16 * sqrt( (nx / 2 - x) * (nx / 2 - x)
                      + (ny / 2 - y) * (ny / 2 - y)
                      + (nz / 2 - z) * (nz / 2 - z) ) / r_max;

            microfield[i] = 0.5
                                       + 2.5508 * r
                                       - 0.588013 * r * r
                                       + 0.032445 * r * r* r
                                       + 0.00124411 * r * r * r * r
                                       - 0.0000973516 * r * r * r * r * r;
            sum += microfield[i];
         }
      }
   }

   // Normalize values in f for effect * dt
   for (int z = 0; z<nz; z++ ) {
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {
            microfield[z*nx*ny + y*nx + x] /= sum;
         }
      }
   }
}

void
phys_sys::InitializeFlow( float initial_flow ) {
   for ( int i = 0; i<n; i++ ) {
      flow[i] = initial_flow;
      f1[i] = 0;
      f2[i] = 0;
   }
}

inline float
phys_sys::RungeKuttaFlow( int x, int y, int z, float f  ) {

   // get iterators along z axis
   int j1 = (z-1)*nx*ny + y*nx + x;
   int j2 = z*nx*ny + y*nx + x;
   int j3 = (z+1)*nx*ny + y*nx + x;

   // calculate my values (mirror boundary conditions)
   float a = 0.0002414;
   float b = 247.8; // C
   float c = 140;
   float my1 = a*exp(b/(temperatures[j2]+273.15-c));
   float my2 = a*exp(b/(temperatures[j2]+273.15-c));
   float my3 = a*exp(b/(temperatures[j2]+273.15-c));
   if ( z != 0 )
      my1 = a*exp(b/(temperatures[j1]+273.15-c));
   if ( z != nz-1 )
      my3 = a*exp(b/(temperatures[j3]+273.15-c));

   // get flow (mirror boundary conditions)
   float flow1 = flow[j2];
   float flow2 = flow[j2];
   float flow3 = flow[j2];
   if ( z != 0 ) 
      flow1 = flow[j1];
   if ( z != nz-1 ) 
      flow3 = flow[j3];
   
   // gravity constant
   float g = 9.81;

   // get epsilon values (mirror boundary conditions)
   float epsilon1=epsilon[j2];
   float epsilon2=epsilon[j2];
   float epsilon3=epsilon[j2];
   if ( z != 0 )
      epsilon1 = epsilon[j1];
   if ( z != nz-1 )
      epsilon3 = epsilon[j3];

   // get previous epsilon values (mirror boundary conditions)
   float prev_eps1=prev_eps[j2];
   float prev_eps2=prev_eps[j2];
   float prev_eps3=prev_eps[j2];
   if ( z != 0 )
      prev_eps1 = prev_eps[j1];
   if ( z != nz-1 )
      prev_eps3 = prev_eps[j3];


   // calculate k values
   float d = 0.000001;
   float k = 0;
   k += flow1 * 1600*(1-epsilon1)*(1-epsilon1)
              /( (epsilon1+d) * (epsilon1+d) * (epsilon1+d));
   k += flow2 * 1600*(1-epsilon2)*(1-epsilon2)
              /( (epsilon2+d) * (epsilon2+d) * (epsilon2+d));
   k += flow3 * 1600*(1-epsilon3)*(1-epsilon3)
              /( (epsilon3+d) * (epsilon3+d) * (epsilon3+d));

   k /= 3;

   // calculate eps values
   float eps = 0;
   eps += epsilon1 - prev_eps1;
   eps += epsilon2 - prev_eps2;
   eps += epsilon3 - prev_eps3;
   eps /= 3;
   eps /= dt;
   // 50/50 chance of +/-
   eps *= (rand()%2) - 1;

   return - ( my1 * my1 * flow1 - 2 * my2 * my2 * flow2 + my3 * my3 * flow3 ) * dt / dz / dz
          - 1 / 2 / dz * (flow3*flow3 - flow1*flow1) + eps ;
}


void
phys_sys::CalculateFlow( ) {

   // Reset f1 and f2
   #pragma omp parallel for
   for ( int i = 0; i < n; i++ ) {
      f1[i] = 0;
      f2[i] = 0;
   }

   // Calculate initial flow
   #pragma omp parallel for
   for ( int x = 0; x < nz; x++ ) {
      for ( int y = dny; y<ny; y++ ) {
	    for ( int z = 0; z < nz; z++ ) {
               int j = z*nx*ny + y*nx + x;
               f1[j] = RungeKuttaFlow(x,y,z,0);
               f2[j] = RungeKuttaFlow(x,y,z,0);
               temp_flow[j] = flow[j] + (f1[j] + f2[j]) / 2 ;
            }
         }
   }


   // Uphold v*e = v*e
   #pragma omp parallel for
   for ( int x = 0; x < nz; x++ ) {
      for ( int y = dny; y<ny; y++ ) {
         for ( int z = 0; z < nz; z++ ) {

            int j2 = z*nx*ny + y*nx + x;
            int j3 = (z+1)*nx*ny + y*nx + x;

            if ( z != nz - 1 ) {

               if ( epsilon[j2] != 0 ) {
                  temp_flow[j2] = temp_flow[j3] * epsilon[j3] / epsilon[j2] ;
	       }
               else if ( epsilon[j3] != 0 ) {
                  temp_flow[j3] = temp_flow[j2] * epsilon[j2] / epsilon[j3] ;
               }
               else {
                  temp_flow[j2] = 0;
               }
            }

         }
      }
   }


   #pragma omp parallel for
   for ( int i = 0; i < n; i++ ) {
      flow[i] = temp_flow[i];
   }


}

void
phys_sys::CalculateBoundaryConditions( float outside_temp ) {

   // The boundary conditions only consider values outside the mesh
   // these values are all in: Air(0) and phase: Gass(4)

   float a = phys_consts::alpha_values[0][4];

   #pragma omp parallel for
   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {

            int i = z*nx*ny + y*nx + x;
            boundarys[i] = 0;

            if ( x == 0 || x == nx-1 )
               boundarys[i] += deltas[2] * a * outside_temp;

            if ( y == 0 || y == ny-1 )
               boundarys[i] += deltas[1] * a * outside_temp;

            if ( z == 0 || z == nz-1 )
               boundarys[i] += deltas[0] * a * outside_temp;

         }
      }
   }
}

void
phys_sys::UpdateAlphaBetaValues( ) {
   float t1 = 40;
   float t2 = 50;

   float* temp = prev_eps;
   prev_eps = epsilon;
   epsilon = temp;


   // For all points in mesh
   #pragma omp parallel for
   for ( int i = 0; i < n; i++ ) {
      alphas[i]=0;
      betas[i]=0;

      int y = ( i / nx ) % ny;

      int substance = 1;

      // Where in the mesh / in what substance are we?
      if ( y>dny ) {
         // We are in fat
         substance=2;
      }
      else {
         // We are in meat
         substance=1;
      }

      // What phase is the substance in?
      for ( int j = 0; j<phys_consts::phases; j++ ) {
         if ( temperatures[ i ] > phys_consts::substance_temps[substance][j] ) {
            alphas[i] = phys_consts::alpha_values[substance][j];
            betas[i] = phys_consts::beta_values[substance][j];
         }
         else {
            break;
         }
      }

      // What is the epsilon value?

      if ( temperatures[i] >= t2 )
         epsilon[i] = 1;
      else if ( temperatures[i] < t1 )
         epsilon[i] = 0;
      else if ( temperatures[i] < t2 )
         epsilon[i] = ( temperatures[i] - t1) / (t2 - t1);

   }
}

void
phys_sys::MultiplyMatrixVector( float* product, float* vector ) {

   #pragma omp parallel for
   for ( int i = 0; i < n; i++ ) {
      product[i] = 0;

      int x0 = i % nx;
      int y0 = ( i / nx ) % ny;
      int z0 = i / nx / ny;

      for ( int t = 0; t < 7; t++ ) {
         int j = i + offsets[t];

         int x = x0 + stencil[t][0];
         int y = y0 + stencil[t][1];
         int z = z0 + stencil[t][2];

         // Check to see if we are inside the mesh
         if ( !isOutside(x,y,z) ) {

            // If we are not on the diagonal
            if ( t != 3 ) {
               product[i] += signs[t]
                           * deltas[t]
                           * alphas[j]
                           * vector[j];
            }
            // If we are on the diagonal
            else {
               product[i] += (1
                           + signs[t]
                           * diagonals[j])
                           * vector[j];
            }
         }
      }
   }
}

void
phys_sys::CalculateDiagonal( ) {

   // The diagonal contains outside elements in addition to
   // the elements in the mesh

   // UpdateAlphaValues() has to be called before this
   // function is called

   // We do not fully calculate the diagonal values, only the part after 1 +/-
   // in the crank nicolson discretization

   for ( int i = 0; i < n; i++ ) {
      diagonals[i] = 0;

      int x0 = i % nx;
      int y0 = ( i / nx ) % ny;
      int z0 = i / nx / ny;

      for ( int t = 0; t < 7; t++ ) {

         int j = i + offsets[t];

         int x = x0 + stencil[t][0];
         int y = y0 + stencil[t][1];
         int z = z0 + stencil[t][2];

         float a;

         if (t != 3) {
            if (isOutside(x,y,z)) {
               // Outside is always Air in Gass form
               a = phys_consts::alpha_values[0][4];
            }
            // If not outside, snatch values from alphas array
            else {
               a = alphas[j];
            }
            diagonals[i] += deltas[t] * a;
         }
      }
   }
}

bool
phys_sys::isOutside(int x, int y, int z) {

   if ( z<0 || z>nz-1 ||
        y<0 || y>ny-1 ||
        x<0 || x>nx-1 ) {
      return true;
   }
   else {
      return false;
   }
}

void
phys_sys::SetLeftSigns( ) {
   signs = (float*)&phys_consts::left_signs;
}

void
phys_sys::SetRightSigns( ) {
   signs = (float*)&phys_consts::right_signs;
}

int phys_sys::offsets[7] = { 0, 0, 0, 0, 0, 0, 0 };
float phys_sys::deltas[7] = { 0, 0, 0, 0, 0, 0, 0 };

int phys_sys::nx=0;
int phys_sys::ny=0;
int phys_sys::nz=0;
int phys_sys::nt=0;
int phys_sys::dny=0;
int phys_sys::n=0;

float phys_sys::dx=0;
float phys_sys::dy=0;
float phys_sys::dz=0;
float phys_sys::dt=0;

float* phys_sys::temperatures=0;
float* phys_sys::alphas=0;
float* phys_sys::betas=0;
float* phys_sys::boundarys=0;
float* phys_sys::diagonals=0;
float* phys_sys::microfield=0;
float* phys_sys::flow=0;
float* phys_sys::temp_flow=0;
float* phys_sys::epsilon=0;
float* phys_sys::prev_eps=0;
float* phys_sys::f1=0;
float* phys_sys::f2=0;

float* phys_sys::signs=0;

int phys_sys::stencil[7][3] = { {0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,0},
                              {1,0,0}, {0,1,0}, {0,0,1} };

