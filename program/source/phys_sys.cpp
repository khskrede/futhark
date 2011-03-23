
#include <math.h>
#include <omp.h>

#include "headers/phys_sys.h"
#include "headers/phys_consts.h"

void
phys_sys::Init(int t_nx, int t_ny, int t_nz, int t_nt, int t_dny,
               float t_dx, float t_dy, float t_dz, float t_dt,
               float* t_temperatures, float* t_alphas, float* t_betas,
               float* t_boundarys, float* t_diagonals, float* t_microfield ) {

   nx=t_nx;
   ny=t_ny;
   nz=t_nz;
   nt=t_nt;
   dny=t_dny;

   n=nx*ny*nz;

   dx=t_dx;
   dy=t_dy;
   dz=t_dz;
   dt=t_dt;

   temperatures=t_temperatures;
   alphas=t_alphas;
   betas=t_betas;
   boundarys=t_boundarys;
   diagonals=t_diagonals;
   microfield=t_microfield;

   offsets[0] = -nx*ny;
   offsets[1] = -nx;
   offsets[2] = -1;
   offsets[3] = 0;
   offsets[4] = 1;
   offsets[5] = nx;
   offsets[6] = nx*ny;

   deltas[0] = dt/2/dz/dz;
   deltas[1] = dt/2/dy/dy;
   deltas[2] = dt/2/dx/dx;
   deltas[3] = 0;
   deltas[4] = dt/2/dx/dx;
   deltas[5] = dt/2/dy/dy;
   deltas[6] = dt/2/dz/dz;

   phys_consts::CalculateValues();
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

            r = 16 * sqrt( (nx / 2 - x) * (nx / 2 - x)
                      + (ny / 2 - y) * (ny / 2 - y)
                      + (nz / 2 - z) * (nz / 2 - z) ) / r_max;

            microfield[z*nx*ny + y*nx + x] = 0.5
                                       + 2.5508 * r
                                       - 0.588013 * r * r
                                       + 0.032445 * r * r* r
                                       + 0.00124411 * r * r * r * r
                                       - 0.0000973516 * r * r * r * r * r;
            sum += microfield[z*nx*ny + y*nx + x];
         }
      }
   }

   // Normalize values in f for effect * dt
   for (int z = 0; z<nz; z++ ) {
      for (int y = 0; y<ny; y++ ) {
         for (int x = 0; x<nx; x++ ) {
            microfield[z*nx*ny + y*nx + x] = microfield[z*nx*ny + y*nx + x] / sum;
         }
      }
   }
}

void
phys_sys::CalculateBoundaryConditions( float outside_temp ) {

   // The boundary conditions only consider values outside the mesh
   // these values are all in: Air(0) and phase: Gass(4)

   float alpha = phys_consts::alpha_values[0][4];

   for ( int z = 0; z < nz; z++ ) {
      for ( int y = 0; y < ny; y++ ) {
         for ( int x = 0; x < nx; x++ ) {

            int i = z*nx*ny + y*nx + x;
            boundarys[ i ] = 0;

            if ( x == 0 || x == nx-1 )
               boundarys[ i ] += deltas[2] * alpha * outside_temp;

            if ( y == 0 || y == ny-1 )
               boundarys[ i ] += deltas[1] * alpha * outside_temp;

            if ( z == 0 || z == nz-1 )
               boundarys[ i ] += deltas[0] * alpha * outside_temp;

         }
      }
   }
}

void
phys_sys::UpdateAlphaBetaValues( ) {

   // Set substance to meat
   int substance = 1;

   // For all points in mesh
   for ( int i = 0; i < n; i++ ) {
      int y = ( i / nx ) % ny;

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
   }
}

void
phys_sys::MultiplyMatrixVector( float* product, float* vector ) {

   for ( int j = 0; j < n; j++ ) {
      product[j] = 0;

      int x0 = j % nx;
      int y0 = ( j / nx ) % ny;
      int z0 = j / nx / ny;

      for ( int t = 0; t < 7; t++ ) {
         int i = j + offsets[t];

         if ( i < n && i>-1 ) {

            int x = x0+stencil[t][0];
            int y = y0+stencil[t][1];
            int z = z0+stencil[t][2];

            // Check to see if we are inside the mesh
            if ( !isOutside(x,y,z) ) {

               // If we are not on the diagonal
               if (t != 3 ) {
                  product[j] += signs[t]
                              * deltas[t]
                              * alphas[i]
                              * vector[i];
               }
               // If we are on the diagonal
               else {
                  product[j] += (1
                              + signs[t]
                              * diagonals[i])
                              * vector[i];
               }
            }
         }
      }
   }
}

void
phys_sys::InitializeTemperature( float initial_temp ) {
   for ( int i = 0; i < n; i++ ) {
      temperatures[i] = initial_temp;
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

inline bool
phys_sys::isOutside(int x, int y, int z) {

   if ( z>=0 && z<nz && y>=0 && y<ny && x>=0 && x <nx ) {
      return false;
   }
   else {
      return true;
   }
}

void
phys_sys::CalculateDiagonal( ) {

   // The diagonal contains outside elements in addition to
   // the elements in the mesh

   // UpdateAlphaValues() and SetSigns() has to be called before this
   // function is called

   // We do not fully calculate the diagonal values, only the part after 1 +/-

   for ( int i = 0; i < n; i++ ) {

      int x0 = i % nx;
      int y0 = ( i / nx ) % ny;
      int z0 = i / nx / ny;

      for ( int t = 0; t < 7; t++ ) {

         int j = i+t;

         int x = x0+stencil[t][0];
         int y = y0+stencil[t][1];
         int z = z0+stencil[t][2];

         float alpha;

         if (isOutside(x,y,z)) {
            // Outside is always Air in Gass form
            alpha = phys_consts::alpha_values[0][4];
         }
         // If not outside, snatch values from alphas array
         else {
            alpha = alphas[j];
         }
         diagonals[i] += deltas[t] * alpha;
      }
   }

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

float* phys_sys::signs=0;

int phys_sys::stencil[7][3] = { {0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,0},
                              {1,0,0}, {0,1,0}, {0,0,1} };
