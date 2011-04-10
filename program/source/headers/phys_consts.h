
#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

class phys_consts {

   private:

   // ---------------
   // Air
   // ---------------

   static const float F_air = 1.1839; // kg/m^3
   static const float Cp_air = 1.012*1000; // J/kg*C
   static const float theta_air = 0.025; // W/m*C

   // ---------------
   // Meat
   // ---------------

   static const float F_meat = 1.042 * 1000; // kg/m^3
   static const float Cp_meat = 1.51*1000; // J/kg*C
   static const float theta_meat = 0.6; // W/m*C

   // ---------------
   // Fat
   // ---------------

   static const float F_fat = 1.000 * 1000; // kg/m^3
   static const float Cp_fat = 1.00*1000; // J/kg*C
   static const float theta_fat = 0.9; // W/m*C
   static const float L_fat = 11.9*1000; // J/kg*C

   public:

   static void CalculateValues( );

   // ----------------
   // 3 substances:
   // ----------------

   // 0: Air
   // 1: Meat
   // 2: Fat

   static const int substances = 3;

   // ----------------
   // 5 phases:
   // ----------------

   // 0: Solid
   // 1: Solid to Liquid
   // 2: Liquid
   // 3: Liquid to gass
   // 4: gass

   static const int phases = 5;

   // Alpha values: Heat equation coefficients
   static float alpha_values[substances][phases];

   // Beta values: Microwave effect coefficients
   static float beta_values[substances][phases];

   // theta_values: Mass transport coefficients
   static float theta_values[substances][phases];

   // Temperatures at state transitions
   static float substance_temps[substances][phases];

   // Signs used in heat equation:
   static const float left_signs[7];
   static const float right_signs[7];
};

#endif
