function [etaf, eta1, eta2, eta3] = emissivities(NGamDot,dV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [etaf, eta1, eta2, eta3] = emissivities(NGamDot,dV)
%
% Script to compute the emissivity source values for a given 
% photon emission rate, NGamDot, and a given cell size, dV.
%
% Daniel R. Reynolds
% 1/17/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
hp = 6.6260693e-27;            % Planck's constant [ergs*s]
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
nu0_HI   = 13.6*ev2erg/hp;     % ionization threshold of HI (hz)
nu0_HeI  = 24.6*ev2erg/hp;     % ionization threshold of HeI (hz)
nu0_HeII = 54.4*ev2erg/hp;     % ionization threshold of HeII (hz)
tol = 1e-10;

% compute the numerical integrals of our output spectrum:
%    \int_0^\infty chi(nu)/h/nu dnu
Llimit = 0;
Ulimit = nu0_HeII;
[I1a,fcnt] = quadl('chinufun',Llimit,Ulimit,tol);
Llimit = 0;
Ulimit = 1;
[I1b,fcnt] = quadl('chinufun2',Llimit,Ulimit,tol);
I1 = I1a + I1b;

%    \int_{nu0_HI}^\infty chi(nu) dnu
Llimit = nu0_HI;
Ulimit = nu0_HeII;
[I2a,fcnt] = quadl('chifun',Llimit,Ulimit,tol);
Llimit = 0;
Ulimit = 1;
[I2b,fcnt] = quadl('chifun2',Llimit,Ulimit,tol);
I2 = I2a + I2b;

% evaluate the spectrum at the ionization thresholds
chi_HI   = chifun(nu0_HI);
chi_HeI  = chifun(nu0_HeI);
chi_HeII = chifun(nu0_HeII);

% form emissivity sources
etaf = NGamDot*I2/(dV*I1);
eta1 = NGamDot*chi_HI/(dV*I1);
eta2 = NGamDot*chi_HeI/(dV*I1);
eta3 = NGamDot*chi_HeII/(dV*I1);

% display components
disp(sprintf('Emissivities:'));
disp(sprintf('   I1 = %12e (%12e + %12e)',I1,I1a,I1b));
disp(sprintf('   I2 = %12e (%12e + %12e)',I2,I2a,I2b));
disp('  ')
disp(sprintf('   etaf = %12e',etaf));
disp(sprintf('   eta1 = %12e',eta1));
disp(sprintf('   eta2 = %12e',eta2));
disp(sprintf('   eta3 = %12e',eta3));
disp('  ')

end
