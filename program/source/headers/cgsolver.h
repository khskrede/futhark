
#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

class cg_solver {

   private:

   void (*MVMultiply)(float* p, float* v);

   int m_n; // size of problem
   int m_nx;
   int m_ny;
   int m_nz;

   float *m_b;
   float *m_x;

   float *m_d;
   float *m_r;

   float *m_Ad;
   float *m_Ax;

   float m_beta;
   float m_alpha;

   public:

   cg_solver(float *b, float *x, int n, void (*MVMultiply)(float*, float*));
   ~cg_solver();

   void Solve();
};

#endif
