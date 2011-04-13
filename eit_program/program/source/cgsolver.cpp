

#include "headers/cgsolver.h"

cg_solver::cg_solver(float *b, float *x, int n, 
                     void (*p_MVMultiply)(float*, float*)) {

   m_n = n;

   m_b = b;
   m_x = x;

   m_d = new float[m_n];
   m_Ad = new float[m_n];
   m_Ax = new float[m_n];
   m_r = new float[m_n];

   m_beta = 0;
   m_alpha = 0;

   MVMultiply = p_MVMultiply;
}

cg_solver::~cg_solver() {

   delete[] m_d;
   delete[] m_r;
   delete[] m_Ad;
   delete[] m_Ax;

}

void
cg_solver::Solve() {

   // Acceptable error
   float error = 0.00000001;

   // variable declarations
   float theta_old = 0, theta_new = 0, theta_0 = 0, temp;
   int iterations = 0;

   // r <= b - Ax
   // d <= r
   MVMultiply(m_Ax, m_x);
   #pragma omp parallel for
   for ( int i = 0; i < m_n; i++ ) {
      m_r[i] = m_b[i] - m_Ax[i];
      m_d[i] = m_r[i];
   }

   // theta_new <= r*r
   #pragma omp parallel for reduction(+:theta_new)
   for ( int i = 0; i < m_n; i++ ) {
      theta_new += m_r[i] * m_r[i];
   }

   // theta_0 <= theta_new
   theta_0 = theta_new;

   while ( theta_new > error * error * theta_0 ) {

      // alpha <= theta_new / d*Ad
      MVMultiply(m_Ad, m_d);
      temp = 0;
      #pragma omp parallel for reduction(+:temp)
      for ( int i = 0; i < m_n; i++) {
         temp += m_d[i] * m_Ad[i];
      }
      m_alpha = theta_new / temp;

      // x <= x + alpha * d
      #pragma omp parallel for
      for ( int i = 0; i < m_n; i++ ) {
         m_x[i] += m_alpha * m_d[i];
      }

      // if 50 iterations have passed
      if ( iterations % 50 == 0 ) {
         // r <= b - Ax
         MVMultiply(m_Ax, m_x);
         #pragma omp parallel for
         for ( int i = 0; i < m_n; i++ ) {
            m_r[i] = m_b[i] - m_Ax[i];
         }
      }

      // else
      else {
         // r <= r - alpha * Ad
         #pragma omp parallel for
         for ( int i = 0; i < m_n; i++ ) {
            m_r[i] -= m_alpha * m_Ad[i] ;
         }
      }

      // theta_old <= theta_new
      theta_old = theta_new;

      // theta_new <= r*r
      theta_new = 0;
      #pragma omp parallel for reduction(+:theta_new)
      for ( int i = 0; i < m_n; i++ ) {
         theta_new += m_r[i] * m_r[i];
      }

      // m_beta <= theta_new / theta_old
      m_beta = theta_new / theta_old;

      // d <= r + beta * d
      #pragma omp parallel for
      for ( int i = 0; i < m_n; i++ ) {
         m_d[i] = m_r[i] + m_beta * m_d[i];
      }

      // iterations <= iterations + 1
      iterations++;
   }
}

