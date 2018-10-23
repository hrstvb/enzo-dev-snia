

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
inline void EOS(float &p, float &rho, float &e, float &h, float &cs, float &dpdrho,
	       float &dpde, int eostype, int mode)
  /* 
     eostype: 
       0: ideal gas
       1: polytropic EOS
       2: another polytropic EOS
       3: isothermal 
       4: pseudocooling for Wengen 2 test
       5: Wengen 2 the original
       6: minimum pressure similar (but not equal to)
          http://adsabs.harvard.edu/abs/2004ApJ...606...32R
          equation 4

       7: polytropic: p = EOSPolytropicFactor * rho^EOSGamma            		//[BH]
       8: polytropic: p = EOSPolytropicFactor * rho^(4/3)				//[BH]
       9: polytropic: p = EOSPolytropicFactor * rho^(5/3)				//[BH]
       10:polytropic: p = EOSPolytropicFactor * rho^( 1 + 1/EOSPolytropicIndex )	//[BH]

     mode:  
       1: given p and rho, calculate others.
       2: given rho and e, calculate others.
  */
{

  float poverrho;

  switch( eostype ) {					//[BH]
  case  7:						//[BH]
  case  8:						//[BH]
  case  9:						//[BH]
  case 10:						//[BH]
    float gamma;					//[BH]

    switch( eostype ) {					//[BH]
    case  7: gamma = EOSGamma; break;			//[BH]
    case  8: gamma = 4.0/3.0;  break;			//[BH]
    case  9: gamma = 5.0/3.0;  break;			//[BH]
    case 10: gamma = 1.0 + 1.0/EOSPolytropicIndex;	//[BH]
    }							//[BH]

    float x, lenu, denu, tu, velu, tempu;		//[BH]
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);	//[BH]
    x /= EOSCriticalDensity/denu;			//[BH]
    p = EOSPolytropicFactor * pow( x, gamma );		//[BH]
    //TODO: Finish the rest of return values:		//[BH]
    e=h=cs=dpdrho=dpde=0;				//[BH]
    return;						//[BH]
  }							//[BH]


  if (eostype == 0) {
    
    if (mode == 1) {
      poverrho = p / rho;
      e = poverrho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e;
      poverrho = p / rho;
    }

    dpdrho = poverrho;
    dpde = (Gamma - 1) * rho;
    h = e + poverrho;
    cs = sqrt(Gamma*poverrho);

  }

  if (eostype == 1) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double c_s = EOSSoundSpeed;
    double rho_cr = EOSCriticalDensity;
    //    c_s /= velu;
    rho_cr /= denu;

    cs = c_s*sqrt(1.0 + EOSGamma*pow(rho/rho_cr, EOSGamma-1.0));
    p = rho*c_s*c_s*(1.0 + pow(rho/rho_cr, EOSGamma-1.0));

    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 2) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double c_s = EOSSoundSpeed;
    double rho_cr = EOSCriticalDensity;
    //    c_s /= velu;
    rho_cr /= denu;
    
    if (rho <= rho_cr) {
      cs = c_s*pow(rho/rho_cr, -0.125);
      p = rho*cs*cs;
    } else {
      cs = c_s*pow(rho/rho_cr, 0.05);
      p = rho*cs*cs;
    }
	
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;

  }

  if (eostype == 3) { // straight isothermal
    cs = EOSSoundSpeed;
    p = rho*cs*cs;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 4) { // Wengen 2 test wants pseudocooling
    cs = EOSSoundSpeed;
    // cooling only to 100 should reduce the resolution requirements
    // for the initial tests
    //    cs  = sqrt(1.e-3 + 1./(1.+pow(rho, 1.5)));
    cs *= sqrt(EOSCriticalDensity + 1./(1.+ rho*sqrt(rho)));
    p = rho * cs*cs;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 5) { // Wengen 2 test wants pseudocooling
    // this is the discontinuous one originally suggested
    cs = EOSSoundSpeed;
    // divided by 1000 is the suggested wengen EOS
    // doing to only 100 should reduce the resolution requirements
    // for the initial tests			
    cs = (rho > 1) ?  cs* sqrt(max(1./(rho*sqrt(rho)), 1.e-3)) : cs ;
    p = rho*cs*cs ;
    e = p / ((Gamma-1.0)*rho);
    dpdrho = 1;
    dpde = 1;
    h = e + p/rho;
  }

  if (eostype == 6) {
    float lenu, denu, tu, velu, tempu;
    GetUnits(&denu, &lenu, &tempu, &tu, &velu, 1);
    double rho_cr = EOSCriticalDensity;
    //    c_s /= velu;
    rho_cr /= denu;
    double Pmin = 1.e4*rho_cr/1.67e-24*1.381e-16/(velu*velu);
    Pmin *= (rho/rho_cr)*(rho/rho_cr); // effective gamma = 2 at high density
    if (mode == 1) {
      poverrho = p / rho;
      e = poverrho / (Gamma - 1);      
    } else if (mode == 2) {
      p = (Gamma - 1) * rho * e + Pmin;
      poverrho = p / rho;
    }
    dpdrho = poverrho;
    
    dpde = (Gamma - 1) * rho;
    
    h = e + poverrho;
    // this is somewhat inconsistent as we are using the gamma = 5/3 everywhere
    cs = sqrt(Gamma*poverrho);

  }



}
