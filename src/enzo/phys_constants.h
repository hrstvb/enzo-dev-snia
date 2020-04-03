#ifndef __phys_constants_h_
#define __phys_constants_h_
/***********************************************************************
/
/ DEFINITION OF PHYSICAL CONSTANTS
/
/ written by: Elizabeth Tasker (renamed by Daniel Reynolds)
/ date:       May, 2005
/
/ Note: CGS units
/
*********************************************************************/

/* Physics constants */

/************************************************/

/* Boltzmann's constant [cm2gs-2K-1] or [ergK-1] */
#define kboltz                          1.3806504e-16

/* Avogadro's number [1/mole] */
#define N_Avog (6.022141e+23)

/* Atomic mass unit, amu [g * mole] */
#define amu (1/N_A)

/* Gas constant, R = kB * NA */
#define R_gas (kboltz * N_Avog)

/* Mass of hydrogen [g] */
#define mh                              1.67262171e-24

/* Mass of an electron [g] */
#define me                              9.10938215e-28

/* Pi */
#define pi                              3.14159265358979323846

/* ergs per eV */
#define erg_eV                          1.602176e-12


/************************************************/

/* Astronomical constant */

/************************************************/

/* Speed of light [cms-1] */

#define clight                          2.99792458e10

/* Gravitational constant [cm3g-1s-2]*/

#define GravConst                       6.67428e-8

/* Solar mass [g] */

#define SolarMass                       1.9891e33

/* Megaparsec [cm] */

#define Mpc                             3.0857e24
#define kpc                             3.0857e21
#define pc                              3.0857e18



/************************************************/

/* Nuclear physics constants */

/************************************************/

#define eV_per_amu  ( 931.494061e6 )		// 1 amu = 931.494061(21)x10^6 eV		//[BH]
#define keV_per_amu ( eV_per_amu/1000 )	// 1 amu = 931.494061(21)x10^3 keV		//[BH]
#define kg_per_amu  ( 1.660538921e-27 )	// 1 amu = 1.660538921(73)x10^-27 kg		//[BH]
#define g_per_amu   ( kg_per_amu * 1000 )	// 1 amu = 1.660538921(73)x10^-24 g		//[BH]
#define eV_per_erg  ( 6.2415e11)  		// 1 erg = 6.2415x10^11 eV			//[BH]
#define erg_per_amu ( eV_per_amu / eV_per_erg )						//[BH]

/*												//[BH]
 # 12C												//[BH]
 */												//[BH]
#define A_12C     ( 12 )									//[BH]
#define M_12C_amu ( 12.0 )	// Atomic Mass: 12.0000000 +- 0.0000000 amu			//[BH]
#define M_12C_cgs ( M_12C_amu * g_per_amu )							//[BH]
				// Binding Energy: 92161.753 +- 0.014 keV			//[BH]
				// Excess Mass: 0.000 +- 0.000 keV				//[BH]
				// Beta Decay Energy: B- -17338.083 +- 1.000 keV		//[BH]

/*												//[BH]
 * O16												//[BH]
 */												//[BH]
#define A_16O	  ( 16 )									//[BH]
#define M_16O_amu ( 15.9949146 )	// Atomic Mass: 15.9949146 +- 0.0000000 amu		//[BH]
#define M_16O_cgs ( M_16O_amu * g_per_amu )							//[BH]
					// Excess Mass: -4736.998 +- 0.001 keV			//[BH]
					// Binding Energy: 127619.336 +- 0.019 keV		//[BH]

/*												//[BH]
 * 56Ni												//[BH]
 */												//[BH]
#define A_56Ni     ( 56 )									//[BH]
#define M_56Ni_amu ( 55.9421363 )	// Atomic Mass: 55.9421363 +- 0.0000119 amu		//[BH]
#define M_56Ni_cgs ( M_56Ni_amu * g_per_amu )							//[BH]
					// Excess Mass: -53899.645 +- 11.130 keV		//[BH]
					// Binding Energy: 483987.827 +- 11.131 keV		//[BH]
					// Beta Decay Energy: B- -15299.000 +- 140.000 keV #	//[BH]




#endif
