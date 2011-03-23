
#include "iostream"
#include "headers/phys_consts.h"

void
phys_consts::CalculateValues( ) {

   std::cout << "Calculating Physical constants \n";

   // ---------------
   // 0: Air
   // ---------------

   // 0:0: solid (irrelevant)
   alpha_values[0][0] = 0;
   beta_values[0][0] = 0;
   substance_temps[0][0] = -2000;

   // 0:1: solid-liquid (irrelevant)
   alpha_values[0][1] = 0;
   beta_values[0][1] = 0;
   substance_temps[0][1] = -2000;

   // 0:2: liquid (irrelevant)
   alpha_values[0][2] = 0;
   beta_values[0][2] = 0;
   substance_temps[0][2] = -2000;

   // 0:3: liquid-gass (irrelevant)
   alpha_values[0][3] = 0;
   beta_values[0][3] = 0;
   substance_temps[0][3] = -2000;

   // 0:4: gass
   alpha_values[0][4] = theta_air / F_air / Cp_air;
   beta_values[0][4] = 1 / F_air / Cp_air;
   substance_temps[0][4] = -196;

   // ---------------
   // 1: Meat
   // ---------------

   // 1:0: solid
   alpha_values[1][0] = theta_meat / F_meat / Cp_meat;
   beta_values[1][0] = 1 / F_meat / Cp_meat;
   substance_temps[1][0] = -200;

   // 1:1: solid-liquid (irrelevant)
   alpha_values[1][1] = 0;
   beta_values[1][1] = 0;
   substance_temps[1][1] = 2000;

   // 1:2: liquid (irrelevant)
   alpha_values[1][2] = 0;
   beta_values[1][2] = 0;
   substance_temps[1][2] = 2000;

   // 1:3: liquid-gass (irrelevant)
   alpha_values[1][3] = 0;
   beta_values[1][3] = 0;
   substance_temps[1][3] = 2000;

   // 1:4: gass (irrelevant)
   alpha_values[1][4] = 0;
   beta_values[1][4] = 0;
   substance_temps[1][4] = 2000;

   // ---------------
   // 2: Fat
   // ---------------

   // 2:0: solid
   alpha_values[2][0] = theta_fat / F_fat / Cp_fat;
   beta_values[2][0] = 1 / F_fat / Cp_fat;
   substance_temps[1][0] = -2000;

   // 2:1: solid-liquid
   alpha_values[2][1] = theta_fat / F_fat / (Cp_fat + L_fat);
   beta_values[2][1] = 1 / F_fat / (Cp_fat + L_fat);
   substance_temps[2][1] = 40;

   // 2:2: liquid
   alpha_values[2][2] = theta_fat / F_fat / Cp_fat;
   beta_values[2][2] = 1 / F_fat / Cp_fat;
   substance_temps[2][2] = 50;

   // 2:3: liquid-gass (irrelevant)
   alpha_values[2][3] = 0;
   beta_values[2][3] = 0;
   substance_temps[2][3] = 2000;

   // 2:4: gass (irrelevant)
   alpha_values[2][4] = 0;
   beta_values[2][4] = 0;
   substance_temps[2][4] = 2000;

}

const float phys_consts::left_signs[7] = { -1, -1, -1, 1, -1, -1, -1 };
const float phys_consts::right_signs[7] = { 1, 1, 1, -1, 1, 1, 1 };

float phys_consts::alpha_values[phys_consts::substances][phys_consts::phases] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
float phys_consts::beta_values[phys_consts::substances][phys_consts::phases] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
float phys_consts::substance_temps[phys_consts::substances][phys_consts::phases] = {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}};
