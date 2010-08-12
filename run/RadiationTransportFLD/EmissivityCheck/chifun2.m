function chi = chifun2(eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: chi = chifun2(eta)
%
% Inputs:  eta - inverse radiation frequency, nu0/nu
%                (can be array-valued)
% Outputs: nu0*chi/eta^2 - scaled spectrum value 
%                (same size as eta)
%
% Daniel R. Reynolds
% 10/9/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the constant
nu0 = 54.4d0*1.60217653e-12/6.6260693e-27;

% evaluate the spectrum
chi = nu0*blackbody(nu0./eta)./eta./eta;



end
