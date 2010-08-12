/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Multi-Frequency, Multi-species, Implicit Problem Class
/  Dump routine
/  
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified1:  
/
/  PURPOSE: Outputs entire problem state to stderr and files.  Useful 
/           upon solve failure.
/
************************************************************************/
#ifdef TRANSFER
#include "MFProb.h"


int MFProb::Dump(EnzoVector *ucur)
{

  if (debug) 
    fprintf(stderr,"\n\nDumping MFProb:\n");

  if (debug) {
    fprintf(stderr,"  told = %g\n",told);
    fprintf(stderr,"  tnew = %g\n",tnew);
    fprintf(stderr,"  dt = %g\n",dt);
    fprintf(stderr,"  theta = %g\n",theta);
    fprintf(stderr,"  LimImp = %"ISYM"\n",LimImp);
    fprintf(stderr,"  LimType = %"ISYM"\n",LimType);
    fprintf(stderr,"  Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"  Model = %"ISYM"\n",Model);
    fprintf(stderr,"  AnalyticChem = %"ISYM"\n",AnalyticChem);
    fprintf(stderr,"  AprxJacobian = %"ISYM"\n",approx_jac);
    fprintf(stderr,"  ESpectrum = %"ISYM"\n",ESpectrum);
    fprintf(stderr,"  aUnits = %g\n",aUnits);
    fprintf(stderr,"  fsUnits = %g\n",fsUnits);
    fprintf(stderr,"  rUnits = %g\n",rUnits);
    fprintf(stderr,"  rUnits0 = %g\n",rUnits0);
    fprintf(stderr,"  eUnits = %g\n",eUnits);
    fprintf(stderr,"  nUnits = %g\n",nUnits);
    fprintf(stderr,"  nUnits0 = %g\n",nUnits0);
    fprintf(stderr,"  rScale = %g\n",rScale);
    fprintf(stderr,"  eScale = %g\n",eScale);
    fprintf(stderr,"  nScale = %g\n",nScale);
    fprintf(stderr,"  DenUnits = %g\n",DenUnits);
    fprintf(stderr,"  DenUnits0 = %g\n",DenUnits0);
    fprintf(stderr,"  LenUnits = %g\n",LenUnits);
    fprintf(stderr,"  LenUnits0 = %g\n",LenUnits0);
    fprintf(stderr,"  TimeUnits = %g\n",TimeUnits);
    fprintf(stderr,"  VelUnits = %g\n",VelUnits);
  }

  if (debug) {
    fprintf(stderr,"Dumping MF module parameters to file MFdump.params\n");
    FILE *fptr = fopen("MFdump.params", "w");
    this->WriteParameters(fptr);
    fclose(fptr);
  }

  if (debug) {
    fprintf(stderr,"EnzoVector indices:\n");
    fprintf(stderr,"   Radiation frequency 1: %"ISYM"\n",iE1);
    fprintf(stderr,"   Radiation frequency 2: %"ISYM"\n",iE2);
    fprintf(stderr,"   Radiation frequency 3: %"ISYM"\n",iE3);
    fprintf(stderr,"   Gas Energy Correction: %"ISYM"\n",iec);
    if (Nchem > 0)
      fprintf(stderr,"              HI Density: %"ISYM"\n",iHI);
    if (Nchem > 1) {
      fprintf(stderr,"             HeI Density: %"ISYM"\n",iHeI);
      fprintf(stderr,"            HeII Density: %"ISYM"\n",iHeII);
    }
  }

  float rmstmp, inftmp;
  for (int ns=0; ns<Nchem+4; ns++) {
    rmstmp = U0->rmsnorm_component(ns);
    inftmp = U0->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"\n  U0(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);

    rmstmp = ucur->rmsnorm_component(ns);
    inftmp = ucur->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"  u(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);

    rmstmp = L_e->rmsnorm_component(ns);
    inftmp = L_e->infnorm_component(ns);
    if (debug)  fprintf(stderr,"  L_e(%"ISYM"): rms = %g, max = %g\n",
			ns,rmstmp,inftmp);

    rmstmp = L_HI->rmsnorm_component(ns);
    inftmp = L_HI->infnorm_component(ns);
    if (debug)  fprintf(stderr,"  L_HI(%"ISYM"): rms = %g, max = %g\n",
			ns,rmstmp,inftmp);

    if (Nchem == 3) {
      rmstmp = L_HeI->rmsnorm_component(ns);
      inftmp = L_HeI->infnorm_component(ns);
      if (debug)  fprintf(stderr,"  L_HeI(%"ISYM"): rms = %g, max = %g\n",
			  ns,rmstmp,inftmp);
      
      rmstmp = L_HeII->rmsnorm_component(ns);
      inftmp = L_HeII->infnorm_component(ns);
      if (debug)  fprintf(stderr,"  L_HeII(%"ISYM"): rms = %g, max = %g\n",
			  ns,rmstmp,inftmp);
    }

  }
  

  char *ofile = new char[12];
  char *tmp_str = new char[3];
  for (int ns=0; ns<Nchem+4; ns++) {

    sprintf(tmp_str,"%"ISYM,ns);

    strcpy(ofile,"u0_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing U0(%"ISYM") to %s\n",ns,ofile);
    U0->writeall(ofile,ns);

    strcpy(ofile,"u_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing u(%"ISYM") to %s\n",ns,ofile);
    ucur->writeall(ofile,ns);

    strcpy(ofile,"L_e_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing L_e(%"ISYM") to %s\n",ns,ofile);
    L_e->writeall(ofile,ns);

    strcpy(ofile,"L_HI_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing L_HI(%"ISYM") to %s\n",ns,ofile);
    L_HI->writeall(ofile,ns);

    if (Nchem == 3) {
      strcpy(ofile,"L_HeI_");
      strcat(strcat(ofile,tmp_str),".vec");
      if (debug) 
	fprintf(stderr,"  writing L_HI(%"ISYM") to %s\n",ns,ofile);
      L_HeI->writeall(ofile,ns);
      
      strcpy(ofile,"L_HeII_");
      strcat(strcat(ofile,tmp_str),".vec");
      if (debug) 
	fprintf(stderr,"  writing L_HeII(%"ISYM") to %s\n",ns,ofile);
      L_HeII->writeall(ofile,ns);
    }

  }

  // output opacities (create temporary vector to do output)
  int empty = 1;
  EnzoVector extras = EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
				 GhDims[0][0], GhDims[0][1], GhDims[1][0], 
				 GhDims[1][1], GhDims[2][0], GhDims[2][1], 
				 1, NBors[0][0], NBors[0][1], NBors[1][0], 
				 NBors[1][1], NBors[2][0], NBors[2][1], 
				 empty);
  extras.SetData(0, piHI);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  piHI: rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"piHI");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing piHI to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, GHI);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  GHI: rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"GHI");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing GHI to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E1_E1);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E1(E1): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E1_E1");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E1(E1) to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E2_E2);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E2(E2): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E2_E2");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E2(E2) to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E3_E3);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E3(E3): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E3_E3");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E3(E3) to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E1_HI);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E1(HI): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E1_HI");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E1(HI) to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E2_HI);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E2(HI): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E2_HI");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E2(HI) to %s\n",ofile);
  extras.writeall(ofile,0);

  extras.SetData(0, L_E3_HI);
  rmstmp = extras.rmsnorm_component(0);
  inftmp = extras.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  L_E3(HI): rms = %g, max = %g\n",rmstmp,inftmp);
  strcpy(ofile,"L_E3_HI");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing L_E3(HI) to %s\n",ofile);
  extras.writeall(ofile,0);

  if (Nchem == 3) {
    extras.SetData(0, piHeI);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  piHeI: rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"piHeI");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing piHeI to %s\n",ofile);
    extras.writeall(ofile,0);
    
    extras.SetData(0, piHeII);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  piHeII: rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"piHeII");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing piHeII to %s\n",ofile);
    extras.writeall(ofile,0);
    
    extras.SetData(0, GHeI);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  GHeI: rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"GHeI");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing GHeI to %s\n",ofile);
    extras.writeall(ofile,0);

    extras.SetData(0, GHeII);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  GHeII: rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"GHeII");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing GHeII to %s\n",ofile);
    extras.writeall(ofile,0);

    extras.SetData(0, L_E2_HeI);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  L_E2(HeI): rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"L_E2_HeI");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing L_E2(HeI) to %s\n",ofile);
    extras.writeall(ofile,0);
    
    extras.SetData(0, L_E3_HeI);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  L_E3(HeI): rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"L_E3_HeI");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing L_E3(HeI) to %s\n",ofile);
    extras.writeall(ofile,0);

    extras.SetData(0, L_E3_HeII);
    rmstmp = extras.rmsnorm_component(0);
    inftmp = extras.infnorm_component(0);
    if (debug) 
      fprintf(stderr,"\n  L_E3(HeII): rms = %g, max = %g\n",rmstmp,inftmp);
    strcpy(ofile,"L_E3_HeII");
    strcat(ofile,".vec");
    if (debug) 
      fprintf(stderr,"  writing L_E3(HeII) to %s\n",ofile);
    extras.writeall(ofile,0);

  }



  // output current residual (store in tmp3)
  this->nlresid(tmp3,ucur);
  for (int ns=0; ns<Nchem+4; ns++) {
    rmstmp = tmp3->rmsnorm_component(ns);
    inftmp = tmp3->infnorm_component(ns);
    if (debug)
      fprintf(stderr,"  f(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp, inftmp);
    sprintf(tmp_str,"%"ISYM,ns);
    strcpy(ofile,"fu_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing f(%"ISYM") to %s\n",ns,ofile);
    tmp3->writeall(ofile,ns);
  }
  delete[] ofile;
  delete[] tmp_str;

  return SUCCESS;
}
#endif
