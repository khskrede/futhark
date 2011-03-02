#include <iostream>

#include "headers/cgsolver.h"

CG_SOLVER::CG_SOLVER(float *b, float *x, int nx, int ny, int nz) {

   m_n = nx*ny*nz;
   m_nx = nx;
   m_ny = ny;
   m_nz = nz;

   m_b = b;
   m_x = x;

   m_d = new float[m_n];
   m_Ad = new float[m_n];
   m_Ax = new float[m_n];
   m_r = new float[m_n];

   m_beta = 0;
   m_alpha = 0;


}

CG_SOLVER::~CG_SOLVER() {

   delete[] m_d;
   delete[] m_r;
   delete[] m_Ad;
   delete[] m_Ax;

}

void
CG_SOLVER::Solve() {

   // Acceptable error
   float error = 0.0;

   // variable declarations
   float theta_old = 0, theta_new = 0, theta_0 = 0, temp;
   int iterations = 0;

   // r <= b - Ax
   // d <= r
   MultiplyDiagAVector(m_Ax, m_x, m_n);
   for ( int i = 0; i < m_n; i++ ) {
      m_r[i] = m_b[i] - m_Ax[i];
      m_d[i] = m_r[i];
   }

   // theta_new <= r*r
   for ( int i = 0; i < m_n; i++ ) {
      theta_new += m_r[i] * m_r[i];
   }

   // theta_0 <= theta_new
   theta_0 = theta_new;
   
   while ( theta_new > error * error * theta_0 ) {

      // alpha <= theta_new / d*Ad
      MultiplyDiagAVector(m_Ad, m_d, m_n);
      temp = 0;
      for ( int i = 0; i < m_n; i++) {
         temp += m_d[i] * m_Ad[i];
      }
      m_alpha = theta_new / temp;

      // x <= x + alpha * d
      for ( int i = 0; i < m_n; i++ ) {
         m_x[i] += m_alpha * m_d[i];
      }

      // if 50 iterations have passed
      if ( iterations % 10 == 0 ) {

         // r <= b - Ax
         MultiplyDiagAVector(m_Ax, m_x, m_n);
         for ( int i = 0; i < m_n; i++ ) {
            m_r[i] = m_b[i] - m_Ax[i];
         }

      }

      // else
      else {

         // r <= r - alpha * Ad
         for ( int i = 0; i < m_n; i++ ) {
            m_r[i] -= m_alpha * m_Ad[i] ;
         }

      }

      // theta_old <= theta_new
      theta_old = theta_new;

      // theta_new <= r*r
      theta_new = 0;
      for ( int i = 0; i < m_n; i++ ) {
         theta_new += m_r[i] * m_r[i];
      }

      // m_beta <= theta_new / theta_old
      m_beta = theta_new / theta_old;

      // d <= r + beta * d
      for ( int i = 0; i < m_n; i++ ) {
         m_d[i] = m_r[i] + m_beta * m_d[i];
      }

      // iterations <= iterations + 1
      iterations++;
   }
}

// Set the required information for the symmetric-diagonal matrix A
void
CG_SOLVER::DiagonalizeA(float *values, int *offsets, int n) {
   m_A_size = n;
   m_A_values = new float[n];
   m_A_offsets = new int[n];
   for ( int i = 0; i < n; i++ ) {
      m_A_values[i] = values[i];
      m_A_offsets[i] = offsets[i];
   }
}

// Set the required information for the symmetric-diagonal matrix B
void
CG_SOLVER::DiagonalizeB(float *values, int *offsets, int n) {
   m_B_size = n;
   m_B_values = new float[n];
   m_B_offsets = new int[n];
   for ( int i = 0; i < n; i++ ) {
      m_B_values[i] = values[i];
      m_B_offsets[i] = offsets[i];
   }
}

// Function used to multiply the diagonal-symmetric matrix A with a vector
void
CG_SOLVER::MultiplyDiagAVector(float *product, float *vector, int n) {
   int i, k;
   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) {
      product[j] = 0;

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int t = 0; t < m_A_size; t++ ) {
         k = m_A_offsets[t];
         i = j - k;

         test = false;
         x = x0; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0-1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0+1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0-1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0+1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0-1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0+1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;


         if (test) {
               product[j] += m_A_values[ k ] * vector[ i ];
         }
         
         i = j + k;
         test = false;
         x = x0; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0-1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0+1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0-1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0+1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0-1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0+1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
                  
         
         if (test) {
               product[j] += m_A_values[ k ] * vector[ i ];
         }
         
         
      }
   }
}

// Function used to multiply the diagonal-symmetric matrix B with a vector
void
CG_SOLVER::MultiplyDiagBVector(float *product, float *vector, int n) {
   int i, k;
   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) {
      product[j] = 0;

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int t = 0; t < m_B_size; t++ ) {
         k = m_B_offsets[t];
         i = j - k;

         test = false;
         x = x0; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0-1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0+1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0-1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0+1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0-1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0+1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;


         if (test) {
               product[j] += m_B_values[ k ] * vector[ i ];
         }
         
         i = j + k;

         test = false;
         x = x0; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0-1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0+1; 
         y = y0;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0-1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0+1;
         z = z0;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0-1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         x = x0; 
         y = y0;
         z = z0+1;
         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
                  
         
         if (test) {
               product[j] += m_B_values[ k ] * vector[ i ];
         }
         
         
      }
   }
}


void
CG_SOLVER::PrintA() {

   int x0, y0, z0, x, y, z;
   bool test;

   for ( int j = 0; j < m_n; j++ ) { 

      x0 = j % m_nx;
      y0 = ( j / m_nx ) % m_ny;
      z0 = j / m_nx / m_ny;

      for ( int i = 0; i < m_n; i++ ) {

         test = false;

         x = x0; 
         y = y0;
         z = z0;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
            
         x = x0-1; 
         y = y0;
         z = z0;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
         
         x = x0+1; 
         y = y0;
         z = z0;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;

         x = x0; 
         y = y0-1;
         z = z0;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
                  
         x = x0; 
         y = y0+1;
         z = z0;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;

         x = x0; 
         y = y0;
         z = z0-1;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;

         x = x0; 
         y = y0;
         z = z0+1;

         if ( ( i == x + y*m_nx + z*m_nx*m_ny ) &&
               !( (x==-1 || x==m_nx) || 
                  (y==-1 || y==m_ny) ||
                  (z==-1 || z==m_nz) ) )
            test = true;
                  
         
         if (test) {
            std::cout << "1 ";
         }
         else {
            std::cout << "0 ";
         }

      }
      std::cout << "\n";
   }
}
