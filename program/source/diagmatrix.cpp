
#include <iostream>
#include "headers/diagmatrix.h"

DIAG_MATRIX::DIAG_MATRIX(float *values, int *offsets, int n, int nx, int ny, int nz) {
   m_size = n;
   m_values = new float[n];
   m_offsets = new int[n];
   for ( int i = 0; i < n; i++ ) {
      m_values[i] = values[i];
      m_offsets[i] = offsets[i];
   }
   m_n = nx*ny*nz;
   m_nx=nx;
   m_ny=ny;
   m_nz=nz;
}

DIAG_MATRIX::~DIAG_MATRIX() {

   delete[] m_values;
   delete[] m_offsets;

}


// Function used to multiply the diagonal matrix with a vector
void
DIAG_MATRIX::MultiplyVector(float *product, float *vector, int n) {

   int disps[7][3] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
                       {-1,0,0}, {0,-1,0}, {0,0,-1} };

   int i, k;
   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) {
      product[j] = 0;

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int t = 0; t < m_size; t++ ) {
         k = m_offsets[t];
         i = j + k;
         test = false;
         
         if ( i < m_n && i >= 0 ) {
         
            test = false;
            for ( int p = 0; p < 7; p++ ) {
               x = x0+disps[p][0]; 
               y = y0+disps[p][1];
               z = z0+disps[p][2];
               if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
                   !( (x==-1 || x==m_nx) || 
                      (y==-1 || y==m_ny) ||
                      (z==-1 || z==m_nz) ) ) {
                  test = true;
                  break;
               }
            }         

            if (test && i>=0 && i<m_n) 
               product[j] += m_values[ t ] * vector[ i ];

         }
      }
   }
}

// Function used to multiply the diagonal matrix with a vector
void
DIAG_MATRIX::Print() {

   int disps[7][3] = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
                       {-1,0,0}, {0,-1,0}, {0,0,-1} };

   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) {

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int i = 0; i < m_n; i++ ) {
         test = false;
         
         if ( i < m_n && i >= 0 ) {
         
            test = false;
            for ( int p = 0; p < 7; p++ ) {
               x = x0+disps[p][0]; 
               y = y0+disps[p][1];
               z = z0+disps[p][2];
               if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
                   !( (x==-1 || x==m_nx) || 
                      (y==-1 || y==m_ny) ||
                      (z==-1 || z==m_nz) ) ) {
                  test = true;
                  break;
               }
            }         

            if (test && i>=0 && i<m_n) {
               std::cout << "1 ";
            }
            else {
               std::cout << "0 ";
            }

         }
      }
      std::cout << "\n";
   }
}
