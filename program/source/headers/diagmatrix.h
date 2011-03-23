
#ifndef DIAGONAL_MATRIX_H
#define DIAGONAL_MATRIX_H

class diag_matrix {

   private:

   int m_n;
   int m_nx;
   int m_ny;
   int m_nz;
   int m_dny;

   int m_substances;
   int m_phases;

   // values used for the matrix A
   float **m_alphas;
   float **m_temps;
   float *m_deltas;
   float *m_signs;

   float *m_values;
   float *m_diagonals;

   int *m_offsets;

   // Positions of x, y, and z in stencil relative to current x, y and z position
   static int m_disps[7][3];

   // Position of delta value in m_deltas array for x, y and z direction
   static int m_delta_off[7];

   public:

   diag_matrix(int substances,
               int phases,
               float* alphas,
               float* deltas,
               float* signs,
               float* temps,
               int* dims,
               int dny,
               int *offsets,
               float *a);

   ~diag_matrix();

   void update( float *temps );

   void multiply(float *product, float *vector, int n);
   void Print();
};


#endif
