function f = GfunHInu(eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: f = GfunHInu(eta)
%
% Inputs: eta = nu0_HI/nu - frequency to evaluate function
% Output: f  - integrand for computing the photo-heating
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
nu = nu0_HI./eta;

% evaluate E
E = Erad(nu);

% evaluate cross-section
sig = zeros(size(nu));
for i=1:length(nu)
   if (nu(i) >= nu0_HI)
      sig(i) = sigHI(nu(i));
   end
end

% combine integrand
f = (nu0_HI./eta./eta).*c.*sig.*E.*(1 - nu0_HI./nu);


% end of function
