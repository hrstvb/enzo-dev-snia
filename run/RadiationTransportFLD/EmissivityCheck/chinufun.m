function fnu = chinufun(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: fnu = chinufun(nu)
%
% Inputs:  nu - radiation frequency (can be array-valued)
% Outputs: fnu - chi/hnu (same size as nu)
%
% Daniel R. Reynolds
% 10/9/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
h = 6.6260693e-27;          % Planck's constant [ergs*s]

% evaluate the components
chival = blackbody(nu);

% evaluate the function
fnu = chival./(nu*h);

end
