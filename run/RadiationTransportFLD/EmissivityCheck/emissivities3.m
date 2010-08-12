function [etaf, eta1, eta2, eta3] = emissivities3(NGamDot,dV,dx,tol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [etaf, eta1, eta2, eta3] = emissivities3(NGamDot,dV,dx,tol)
%
% Script to compute the emissivity source values for a given 
% photon emission rate, NGamDot, and a given cell size, dV.
%
% This routine uses a fixed mesh width dx to integrate the 
% functions with a composite midpoint rule (immune to the 
% singularity at nu=0), continuing along until the successive 
% integral updates are below 1e-4*tol/integral.
%
% Daniel R. Reynolds
% 1/17/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
h_nu1 = 13.6d0*ev2erg;         % ionization energy of HI [ergs]
h_nu2 = 24.6d0*ev2erg;         % ionization energy of HeI [ergs]
h_nu3 = 54.4d0*ev2erg;         % ionization energy of HeII [ergs]
hp = 6.6260693e-27;            % Planck's constant (ergs*s)
nu0_HI   = h_nu1/hp;           % ionization threshold of HI (hz)
nu0_HeI  = h_nu2/hp;           % ionization threshold of HeI (hz)
nu0_HeII = h_nu3/hp;           % ionization threshold of HeII (hz)

% get the radiation cross section at each frequency
chi1 = chifun(nu0_HI);
chi2 = chifun(nu0_HeI);
chi3 = chifun(nu0_HeII);

% compute the scaling factor for the output emissivity
chinuint = 0;
nu = dx/2;
finished = 0;
while (finished < 1) 
   chinuval = chifun(nu)/nu/hp;
   addval = dx*chinuval;
   chinuint = chinuint + addval;
   if (addval < tol*1e-4*max(chinuint))
      finished = 1;
   end
   nu = nu + dx;
end
disp(sprintf('chinuint: nu_final/nu0_HeII = %12e',nu/nu0_HeII))

% compute the scaling factor for the free-streaming emissivity
chiint = 0;
nu = nu0_HI + dx/2;
finished = 0;
while (finished < 1)
   chival = chifun(nu);
   addval = dx*chival;
   chiint = chiint + addval;
   if (addval < tol*1e-4*max(chiint))
      finished = 1;
   end
   nu = nu + dx;
end
disp(sprintf('chiint: nu_final/nu0_HeII = %12e',nu/nu0_HeII))


% compute the effective NGammaDot for the free-streaming emissivity
NGamDot_FS = NGamDot*chiint/chinuint/h_nu1;

% compute eta factors for given ionization source
etaf = h_nu1*NGamDot_FS/dV;
eta1 = NGamDot*chi1/chinuint/dV;
eta2 = NGamDot*chi2/chinuint/dV;
eta3 = NGamDot*chi3/chinuint/dV;

        
% display components
disp(sprintf('Emissivities:'));
disp(sprintf('   chinuint = %12e',chinuint));
disp(sprintf('   chiint   = %12e',chiint));
disp('  ')
disp(sprintf('   etaf = %12e',etaf));
disp(sprintf('   eta1 = %12e',eta1));
disp(sprintf('   eta2 = %12e',eta2));
disp(sprintf('   eta3 = %12e',eta3));
disp('  ')

end
