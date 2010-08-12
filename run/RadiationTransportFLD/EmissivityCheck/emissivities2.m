function [etaf, eta1, eta2, eta3] = emissivities2(NGamDot,dV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: [etaf, eta1, eta2, eta3] = emissivities2(NGamDot,dV)
%
% Script to compute the emissivity source values for a given 
% photon emission rate, NGamDot, and a given cell size, dV.
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
nnus = 99;                     % number of quadrature nodes

% get the radiation cross section at each frequency
chi1 = chifun(nu0_HI);
chi2 = chifun(nu0_HeI);
chi3 = chifun(nu0_HeII);

% compute the scaling factor for the output emissivity
%    set integration nodes
Llimit = 0.1;
Ulimit = 1-sqrt(eps);
for l=1:nnus
   nus(l,1) = 1 + (l-1)*(nu0_HeII-1)/(nnus-1);
%   nus(l,1) = 1 + (l-1)*(50*nu0_HeII-nu0_HI)/(nnus-1);
   etas(l,1) = Llimit + (l-1)*(Ulimit-Llimit)/(nnus-1);
   nusB(l,1) = nu0_HeII/etas(l,1);
end
chis = chifun(nus);
chisB = chifun(nusB);
%    set the quadrature weights (composite Simpson)
nbins = (nnus-1)/2;
wts = zeros(nnus,1);
wtsB = zeros(nnus,1);
for l=1:nbins
   i = 2*l-1;
   j = 2*l;
   k = 2*l+1;
   wts(i) = wts(i) +   (nus(k)-nus(i));
   wts(j) = wts(j) + 4*(nus(k)-nus(i));
   wts(k) = wts(k) +   (nus(k)-nus(i));
   wtsB(i) = wtsB(i) +   (etas(k)-etas(i))*nu0_HeII/etas(i)^2;
   wtsB(j) = wtsB(j) + 4*(etas(k)-etas(i))*nu0_HeII/etas(j)^2;
   wtsB(k) = wtsB(k) +   (etas(k)-etas(i))*nu0_HeII/etas(k)^2;
end
wts = wts/6;  wtsB = wtsB/6;
chinuinta = sum(wts.*chis./nus/hp);
chinuintb = sum(wtsB.*chisB./nusB/hp);
chinuint = chinuinta + chinuintb;
%chinuint = chinuinta;

%    set integration nodes for free-streaming emissivity
for l=1:nnus
   nus(l,1) = nu0_HI + (l-1)*(nu0_HeII-nu0_HI)/(nnus-1);
end
chis = chifun(nus);
%    set the quadrature weights (composite Simpson)
wts = zeros(nnus,1);
for l=1:nbins
   i = 2*l-1;
   j = 2*l;
   k = 2*l+1;
   wts(i) = wts(i) +   (nus(k)-nus(i));
   wts(j) = wts(j) + 4*(nus(k)-nus(i));
   wts(k) = wts(k) +   (nus(k)-nus(i));
end
wts = wts/6;
chiinta = sum(wts.*chis);
chiintb = sum(wtsB.*chisB);
chiint = chiinta + chiintb;


% compute the effective NGammaDot for the free-streaming emissivity
NGamDot_FS = NGamDot*chiint/chinuint/h_nu1;

% compute eta factors for given ionization source
etaf = h_nu1*NGamDot_FS/dV;
eta1 = NGamDot*chi1/chinuint/dV;
eta2 = NGamDot*chi2/chinuint/dV;
eta3 = NGamDot*chi3/chinuint/dV;

        
% display components
disp(sprintf('Emissivities:'));
disp(sprintf('   chinuint = %12e (%12e + %12e)',chinuint,chinuinta,chinuintb));
disp(sprintf('   chiint   = %12e (%12e + %12e)',chiint,chiinta,chiintb));
disp('  ')
disp(sprintf('   etaf = %12e',etaf));
disp(sprintf('   eta1 = %12e',eta1));
disp(sprintf('   eta2 = %12e',eta2));
disp(sprintf('   eta3 = %12e',eta3));
disp('  ')

end
