
#include "headers/cgsolver.h"

CG_SOLVER::CG_SOLVER(float *b, float *x, int nx, int ny, int nz, DIAG_MATRIX *A) {

   m_n = nx*ny*nz;
   m_nx = nx;
   m_ny = ny;
   m_nz = nz;

   m_A = A;

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
   float error = 0.00001;

   // variable declarations
   float theta_old = 0, theta_new = 0, theta_0 = 0, temp;
   int iterations = 0;

   // r <= b - Ax
   // d <= r
   (*m_A).MultiplyVector(m_Ax, m_x, m_n);
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
      (*m_A).MultiplyVector(m_Ad, m_d, m_n);
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
         (*m_A).MultiplyVector(m_Ax, m_x, m_n);
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

