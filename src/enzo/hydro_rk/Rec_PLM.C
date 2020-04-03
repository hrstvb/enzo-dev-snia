/***********************************************************************
/
/  PLM RECONSTRUCTION
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <exception>
#include "../DebugMacros.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ReconstructionRoutines.h"
#include "fortran.def"


float plm_lr(int leftright, float* field, size_t i, float* rhoField)
{
	size_t il = (leftright==-1) ? (i - 1) : (i + 2) ;
	size_t im = il + leftright;
	size_t ir = im + leftright;
	float vm1 = field[il];
	float v = field[im];
	float vp1 = field[ir];

	if(rhoField)
	{
		vm1 /= rhoField[il];
		v /= rhoField[im];
		vp1 /= rhoField[ir];
	}

	float dv_l = (v - vm1) * Theta_Limiter;
	float dv_m = 0.5*(vp1-vm1);
	float dv_r = (vp1 - v) * Theta_Limiter;

	float dv = minmod(dv_l, dv_r, dv_m);
	return v - leftright * 0.5 * dv;
}

inline float plm_l(float vm1, float v, float vp1)
{

  float dv_l, dv_r, dv_m, dv;

  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);

  dv = minmod(dv_l, dv_r, dv_m);

  return v + 0.5*dv;
}

inline float plm_r(float vm1, float v, float vp1)
{

  float dv_l, dv_r, dv_m, dv;

  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);

  dv = minmod(dv_l, dv_r, dv_m);

  return v - 0.5*dv;
}

inline float plm_point(float vm1, float v, float vp1)
{

  float dv_l, dv_r, dv_m, dv;

  dv_l = (v-vm1) * Theta_Limiter;
  dv_r = (vp1-v) * Theta_Limiter;
  dv_m = 0.5*(vp1-vm1);

  dv = minmod(dv_l, dv_r, dv_m);

  return v + 0.5*dv;

}

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq, int debug)
{
  int iprim;

  const int offset = NumberOfGhostZones - 1;
    for (int field = 0; field < Neq + (bool)(UseBurning); field++) {
      iprim = offset;
      for (int i = 0; i < ActiveSize+1; i++, iprim++) {
        priml[field][i] = plm_point(prim[field][iprim-1], prim[field][iprim  ], prim[field][iprim+1]);
        primr[field][i] = plm_point(prim[field][iprim+2], prim[field][iprim+1], prim[field][iprim]);
      }
  }

  for (int i = 0; i < ActiveSize+1; i++) {
    priml[0][i] = max(priml[0][i], SmallRho);
    //priml[1][i] = max(priml[1][i], SmallEint);
    primr[0][i] = max(primr[0][i], SmallRho);
    //primr[1][i] = max(primr[1][i], SmallEint);
  }

  return SUCCESS;
}

int plm_species(float **prim, int is, float **species, float *flux0, int ActiveSize)
  /* is : starting index of species field in prim */
{

  int iprim;
  const int offset = NumberOfGhostZones - 1;
  static float sum[MAX_ANY_SINGLE_DIRECTION];

  for (int field = 0; field < NSpecies; field++) {
    iprim = offset;
    for (int n = 0; n < ActiveSize+1; n++, iprim++) {
      if (flux0[n] >= 0) {
	species[field][n] = plm_l(prim[field+is][iprim-1], prim[field+is][iprim], prim[field+is][iprim+1]);
      } else {
	species[field][n] = plm_r(prim[field+is][iprim], prim[field+is][iprim+1], prim[field+is][iprim+2]);
      }
    }
  } // ENDFOR field

  /* renormalize species field */

  if (NoMultiSpeciesButColors != TRUE) {
    for (int n = 0; n < ActiveSize+1; n++)
      sum[n] = 0.0f;
    for (int field = 0; field < NSpecies; field++) {
      for (int n = 0; n < ActiveSize+1; n++)
        sum[n] += species[field][n];
    }

    for (int field = 0; field < NSpecies; field++) {
      for (int n = 0; n < ActiveSize+1; n++)
        species[field][n] /= sum[n];
    }
  }

  return SUCCESS;

}

int plm_color(float **prim, int is, float **color, float *flux0, int ActiveSize)
  /* is : starting index of *species* field in prim */
{

  int iprim;
  const int offset = NumberOfGhostZones - 1;
  for (int field = is+NSpecies; field < is+NSpecies+NColor; field++) {
    iprim = offset;
    for (int n = 0; n < ActiveSize+1; n++, iprim++) {
      if (flux0[n] >= 0) {
	color[field-is-NSpecies][n] = plm_l(prim[field][iprim-1], prim[field][iprim], prim[field][iprim+1]);
      } else {
	color[field-is-NSpecies][n] = plm_r(prim[field][iprim], prim[field][iprim+1], prim[field][iprim+2]);
    }
    }
  }

  return SUCCESS;

}

int plm(float **prim, float **priml, float **primr, int ActiveSize, int Neq)
{
	return plm(prim, priml, primr, ActiveSize, Neq, 0);
}
