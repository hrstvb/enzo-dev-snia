function sig = sigHeI(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: sig = sigHeI(nu)
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
nu0 = 24.6*ev2erg/h;        % ionization threshold of Hydrogen (hz)


% evaluate the cross-section
sig = zeros(size(nu));
for i=1:length(nu)
   if (nu(i) <= nu0)
      sig(i) = 7.42e-18;
   else
      sig(i) = 7.42e-18 * (1.66*(nu0/nu(i))^(2.05) - 0.66*(nu0/nu(i))^(3.05));
   end
end

                
