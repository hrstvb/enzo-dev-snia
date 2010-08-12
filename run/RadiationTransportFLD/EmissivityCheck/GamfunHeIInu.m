function f = GamfunHeIInu(eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: f = GamfunHeIInu(eta)
%
% Inputs: eta = nu0_HeII/nu - frequency to evaluate function
% Output: f  - integrand for computing the photo-ionization 
%              rate
%
% Daniel R. Reynolds
% 10.15.2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some parameters
hp = 6.6260693e-27;            % Planck's constant (ergs*s)
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
nu0_HI   = 13.6*ev2erg/hp;     % ionization threshold of HI (hz)
nu0_HeI  = 24.6*ev2erg/hp;     % ionization threshold of HeI (hz)
nu0_HeII = 54.4*ev2erg/hp;     % ionization threshold of HeII (hz)

% get nu
nu = nu0_HeII./eta;

% evaluate E
E = Erad(nu);

% evaluate cross-section
sig = zeros(size(nu));
for i=1:length(nu)
   if (nu(i) >= nu0_HeII)
      sig(i) = sigHeII(nu(i));
   end
end

% combine integrand
f = (nu0_HeII./eta./eta).*c.*sig.*E./hp./nu;


% end of function
