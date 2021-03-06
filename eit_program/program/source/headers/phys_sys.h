
#ifndef PHYSICAL_SYSTEM_H
#define PHYSICAL_SYSTEM_H

class phys_sys {

   private:

   static int nx, ny, nz, nt, dny, n;
   static float dx, dy, dz, dt;

   static float* temperatures;
   static float* alphas;
   static float* betas;
   static float* boundarys;
   static float* diagonals;
   static float* microfield;
   static float* flow;
   static float* temp_flow;
   static float* epsilon;
   static float* prev_eps;
   static float* f1;
   static float* f2;

   static int stencil[7][3];
   static int offsets[7];
   static float deltas[7];

   static float* signs;

   // Helper function
   static bool
   isOutside(int x, int y, int z);

   public:

   static void
   Init( int nx, int ny, int nz, int nt, int dny,
         float dx, float dy, float dz, float dt,
         float* temperatures, float* alphas, float* betas,
         float* boundarys, float* diagonals, float* microfield,
         float* flow, float* temp_flow, float* epsilon, float* f1, float* f2 , float* prev_eps);

   static void
   CalculateWaveField( float effect );

   static void
   InitializeFlow( float flow );

   inline static float
   RungeKuttaFlow( int x, int y, int z, float f );

   static void
   CalculateFlow();

   static void
   InitializeTemperature( float initial_temp );

   static void
   CalculateBoundaryConditions( float outside_temp );

   static void
   CalculateDiagonal( );

   static void
   UpdateAlphaBetaValues( );

   static void
   SetLeftSigns( );

   static void
   SetRightSigns( );

   static void
   MultiplyMatrixVector( float *product,
                         float *vector );

};

#endif
