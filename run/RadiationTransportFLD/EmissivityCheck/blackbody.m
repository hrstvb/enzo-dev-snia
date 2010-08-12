function chi = blackbody(nu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: chi = blackbody(nu)
%
% Inputs:  nu - radiation frequency [hz] (can be array-valued)
% Outputs: chi - spectrum value (same size as nu)
%
% Daniel R. Reynolds
% 10/9/2009
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some constants
h = 6.6260693e-27;          % Planck's constant [ergs*s]
kb = 1.3806504e-16;         % Boltzmann's constant [ergs/K]
c = 2.99792458e10;          % speed of light [cm/s]
ev2erg = 1.60217653e-12;    % conversion constant from eV to ergs
nu0 = 13.6*ev2erg/h;        % ionization threshold of Hydrogen (hz)
T = 1e5;                    % blackbody source temperature [K]

% evaluate the spectrum
chi = 8*pi*h*(nu/c).^3./(exp(h*nu/kb/T)-1);
