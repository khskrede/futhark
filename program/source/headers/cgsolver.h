
#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

class CG_SOLVER {

   private:


   int m_n; // size of problem
   int m_nx;
   int m_ny;
   int m_nz;

   float *m_b; // known vector of size n.
   float *m_x; // unknown vector of size n.

   // values used for the matrix A
   float *m_A_values;
   int *m_A_offsets;
   int m_A_size;
   
   // values used for the matrix A
   float *m_B_values;
   int *m_B_offsets;
   int m_B_size;

   float *m_d;
   float *m_r;

   float *m_Ad;
   float *m_Ax;

   float m_beta;
   float m_alpha;

   public:

   CG_SOLVER(float* b, float* x, int nx, int ny, int nz);
   ~CG_SOLVER();

   void Solve();
   void DiagonalizeA(float *values, int *offsets, int n);
   void DiagonalizeB(float *values, int *offsets, int n);
   void MultiplyDiagAVector(float *product, float *vector, int n);
   void MultiplyDiagBVector(float *product, float *vector, int n);
   void PrintA();
};

#endif
