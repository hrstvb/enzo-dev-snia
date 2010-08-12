function chi = chifun(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: chi = chifun(nu)
%
% Inputs:  nu - radiation frequency [hz] (can be array-valued)
% Outputs: chi - spectrum value (same size as nu)
%
% Daniel R. Reynolds
% 10/9/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evaluate the spectrum
chi = blackbody(nu);


end