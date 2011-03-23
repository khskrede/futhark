
#include <iostream>
#include <omp.h>
#include "headers/diagmatrix.h"

diag_matrix::diag_matrix(int substances,
                         int phases,
                         float* alphas,
                         float* deltas,
                         float* signs,
                         float* temps,
                         int* dims,
                         int dny,
                         int *offsets,
                         float *a) {

   m_substances = substances;
   m_phases = phases;

   m_alphas = new float*[substances];
   m_temps = new float*[substances];
   for (int i = 0; i<substances; i++) {
      m_alphas[i] = new float[phases];
      m_temps[i] = new float[phases];
      for ( int j = 0; j < phases; j++ ) {
         m_alphas[i][j] = alphas[i*substances+j];
         m_temps[i][j] = alphas[i*substances+j];
      }
   }

   m_deltas = deltas;
   m_signs = signs;
   m_dny = dny;

   m_nx = dims[0];
   m_ny = dims[2];
   m_nz = dims[3];

   m_n = m_nx*m_ny*m_nz;

   m_values = a;
}

diag_matrix::~diag_matrix() {

   for (int i = 0; i<m_substances; i++) {
      delete m_temps[i];
      delete m_alphas[i];
   }
   delete[] m_alphas;
   delete[] m_temps;
   delete[] m_values;

}

void
diag_matrix::update( float* temps ) {

   int substance = 0;

   // Update matrix with alpha values
   for ( int i = 0; i < m_n; i++ ) {

      int y = ( i / m_nx ) % m_ny;

      if ( y>m_dny ) {
         substance=1;
      }
      else {
         substance=0;
      }

      float temp = temps[ i ];

      for ( int j = 0; j<m_phases; j++ ) {
         if ( temp < m_temps[substance][j] ) {
            m_values[i] = m_alphas[substance][j];
         }
         else {
            break;
         }
      }
   }

   // Update matrix with diagonal values
   for( int j = 0; j < m_n; j++ ) {
      for ( int t = 0; t < 7; t++ ) {
         int i = j + m_offsets[t];

         m_diagonals[j] += m_deltas[m_delta_off[t]]
                         * m_values[i];
      }
      m_diagonals[j] = 1 + m_signs[3] * m_diagonals[j];
   }
}

// Function used to multiply the diagonal matrix with a vector
void
diag_matrix::multiply(float *product, float *vector, int n) {

   #pragma omp parallel for
   for ( int j = 0; j < m_n; j++ ) {
      product[j] = 0;

      int i;
      int x0, y0, z0, x, y, z;
      bool test;

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int t = 0; t < 7; t++ ) {
         i = j + m_offsets[t];

         if ( i < m_n && i>-1 ) {

            test = false;

            x = x0+m_disps[t][0];
            y = y0+m_disps[t][1];
            z = z0+m_disps[t][2];

            // Check to see if we are inside the mesh
            if ( (x!=-1 && x!=m_nx) &&
                 (y!=-1 && y!=m_ny) &&
                 (z!=-1 && z!=m_nz) ) {
               test = true;
            }

            if (test) {
               // If we are not the diagonal
               if (m_delta_off[t] > -1) {
                  product[j] += m_signs[t]
                              * m_deltas[m_delta_off[t]]
                              * m_values[i]
                              * vector[i];
               }
               else {
                  product[j] += m_signs[t]
                              * m_diagonals[i]
                              * m_values[i]
                              * vector[i];
               }
            }
         }
      }
   }
}

int
diag_matrix::m_disps[7][3] = { {0,0,-1}, {0,-1,0}, {-1,0,0}, {0,0,0},
                       {1,0,0}, {0,1,0}, {0,0,1} };

int
diag_matrix::m_delta_off[7] = { 2, 1, 0, -1, 0, 1, 2 };


/*
// Function used to multiply the diagonal matrix with a vector
void
diag_matrix::Print() {

   int disps[7][3] = { {-1,0,0}, {0,-1,0}, {0,0,-1}, {0,0,0},
                       {0,0,1}, {0,1,0}, {1,0,0} };

   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) {

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int i = 0; i < m_n; i++ ) {
         test = false;

         if ( true ) {

            test = false;
            for ( int p = 0; p < 7; p++ ) {
               x = x0+disps[p][0];
               y = y0+disps[p][1];
               z = z0+disps[p][2];
               if ( (i == x + y*m_nx + z*m_nx*m_ny) &&
                    !((x==-1 || x==m_nx) ||
                      (y==-1 || y==m_ny) ||
                      (z==-1 || z==m_nz)) ) {
                  test = true;
                  break;
               }
            }

            if (test) {
               std::cout << "1 " ;
            }
            else {
               std::cout << "0 ";
            }

         }
      }
      std::cout << "\n";
   }
} */


