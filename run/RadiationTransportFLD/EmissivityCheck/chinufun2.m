function fnu = chinufun2(eta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: fnu = chinufun(eta)
%
% Inputs:  eta - inverse radiation frequency, nu0/nu
%                (can be array-valued)
% Outputs: chi/eta - scaled spectrum value 
%                (same size as eta)
%
% Daniel R. Reynolds
% 10/9/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
h = 6.6260693e-27;          % Planck's constant [ergs*s]
nu0 = 54.4d0*1.60217653e-12/6.6260693e-27;

% evaluate the components
chival = blackbody(nu0./eta);

% evaluate the function
fnu = chival./(eta*h);

end
