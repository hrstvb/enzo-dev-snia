function sig = sigHeII(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: sig = sigHeII(nu)
%
% Inputs:  nu - radiation frequency [hz] (array-valued)
% Outputs: sig - cross-section value (same size as nu)
%
% Daniel R. Reynolds
% 10/14/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
h = 6.6260693e-27;          % Planck's constant [ergs*s]
ev2erg = 1.60217653e-12;    % conversion constant from eV to ergs
nu0 = 54.4*ev2erg/h;        % ionization threshold of Hydrogen (hz)


% evaluate the cross-section
sig = zeros(size(nu));
for i=1:length(nu)
   if (nu(i) <= nu0)
      sig(i) = 1.575e-18;
   else
      eps = sqrt(nu(i)/nu0 - 1);
      sig(i) = 1.575d-18 * (nu0/nu(i))^4 * exp(4-4*atan(eps)/eps) / (1-exp(-2*pi/eps));
   end
end

                
