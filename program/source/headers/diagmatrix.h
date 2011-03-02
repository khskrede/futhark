
#ifndef DIAG_MATRIX_H
#define DIAG_MATRIX_H

class DIAG_MATRIX {

   private:


   int m_n; // size of problem
   int m_nx;
   int m_ny;
   int m_nz;

   // values used for the matrix A
   float *m_values;
   int *m_offsets;
   int m_size;

   public:

   DIAG_MATRIX(float *values, int *offsets, int n, int nx, int ny, int nz);
   ~DIAG_MATRIX();

   void MultiplyVector(float *product, float *vector, int n);
   void Print();
};

#endif
