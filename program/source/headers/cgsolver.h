
#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "diagmatrix.h"

class CG_SOLVER {

   private:

   DIAG_MATRIX* m_A;

   int m_n; // size of problem
   int m_nx;
   int m_ny;
   int m_nz;

   float *m_b; // known vector of size n.
   float *m_x; // unknown vector of size n.

   float *m_d;
   float *m_r;

   float *m_Ad;
   float *m_Ax;

   float m_beta;
   float m_alpha;

   public:

   CG_SOLVER(float *b, float *x, int nx, int ny, int nz, DIAG_MATRIX *A);
   ~CG_SOLVER();

   void Solve();
};

#endif
