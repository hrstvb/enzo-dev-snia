function Eg = Egrey(Ef,E1,E2,E3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage:  E = Egrey(Ef,E1,E2,E3)
%
% Inputs: Ef - free-streaming radiation values (3D, CGS)
%         E1 - radiation frequency 1 values (3D, CGS)
%         E2 - radiation frequency 2 values (3D, CGS)
%         E3 - radiation frequency 3 values (3D, CGS)
% Output: Eg - integrated radiation energy density 
%
% Daniel R. Reynolds
% 7.8.2010
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize parameters & constants
nnus = 1001;                   % # of frequencies in each bin (must be odd)
hp = 6.6260693e-27;            % Planck's constant (ergs*s)
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
nu0_HI   = 13.6*ev2erg/hp;     % ionization threshold of HI (hz)
nu0_HeI  = 24.6*ev2erg/hp;     % ionization threshold of HeI (hz)
nu0_HeII = 54.4*ev2erg/hp;     % ionization threshold of HeII (hz)
Llimit = 0.1;                  % lower limit of int (shift away from 0)
Ulimit = 1 - eps^(0.5);        % upper limit of int (shift away from 1)

% set frequencies for evaluation
nusA = zeros(nnus,1);
nusB = zeros(nnus,1);
nusC = zeros(nnus,1);
etasC = zeros(nnus,1);
for l=1:nnus
   % nusA = linspace(nu0_HI,nu0_HeI,nnus)
   nusA(l) = nu0_HI + (l-1)*(0.9999999999*nu0_HeI-nu0_HI)/(nnus-1);
   % nusB = linspace(nu0_HeI,nu0_HeII,nnus)
   nusB(l) = nu0_HeI + (l-1)*(0.9999999999*nu0_HeII-nu0_HeI)/(nnus-1);
   % etasC = linspace(Llimit,Llimit,nnus)
   etasC(l) = Llimit + (l-1)*(Ulimit-Llimit)/(nnus-1);
   % nusC = nu0_HeII./etasC
   nusC(l) = nu0_HeII/etasC(l);
end

% evaluate frequency-dependent functions at nu values
%    shortcuts for cross-sections at ionization thresholds
s1n1 = sigHI(nu0_HI);
s1n2 = sigHI(nu0_HeI)/s1n1;
s1n3 = sigHI(nu0_HeII)/s1n1;
s2n2 = sigHeI(nu0_HeI);
s2n3 = sigHeI(nu0_HeII)/s2n2;
s3n3 = sigHeII(nu0_HeII);

%    species cross sections across intervals
fHIA   = sigHI(nusA);
fHIB   = sigHI(nusB);
fHIC   = sigHI(nusC);
fHeIB  = sigHeI(nusB);
fHeIC  = sigHeI(nusC);
fHeIIC = sigHeII(nusC);

%    radiation spectrum across intervals
chiA = chifun(nusA);
chiB = chifun(nusB);
chiC = chifun(nusC);

%    exponents for radiation terms
S3C = fHeIIC/s3n3;
S2C = fHeIC/s2n2 - s2n3*S3C;
S2B = fHeIB/s2n2;
S1C = fHIC/s1n1  - s1n3*S3C - s1n2*S2C;
S1B = fHIB/s1n1  - s1n2*S2B;
S1A = fHIA/s1n1;

% compute the quadrature weights across the intervals
wtsA = zeros(nnus,1);
wtsB = zeros(nnus,1);
wtsC = zeros(nnus,1);
nbins = (nnus-1)/2;
for l=1:nbins
   i = 2*l-1;
   j = 2*l;
   k = 2*l+1;
   wtsA(i) = wtsA(i) + (nusA(k)-nusA(i))/6;
   wtsA(j) = wtsA(j) + (nusA(k)-nusA(i))*2/3;
   wtsA(k) = wtsA(k) + (nusA(k)-nusA(i))/6;
   wtsB(i) = wtsB(i) + (nusB(k)-nusB(i))/6;
   wtsB(j) = wtsB(j) + (nusB(k)-nusB(i))*2/3;
   wtsB(k) = wtsB(k) + (nusB(k)-nusB(i))/6;
   wtsC(i) = wtsC(i) + nu0_HeII/etasC(i)^2*(etasC(k)-etasC(i))/6;
   wtsC(j) = wtsC(j) + nu0_HeII/etasC(j)^2*(etasC(k)-etasC(i))*2/3;
   wtsC(k) = wtsC(k) + nu0_HeII/etasC(k)^2*(etasC(k)-etasC(i))/6;
end

% compute the integrated radiation spectrum
chibar = sum(wtsA.*chiA) + sum(wtsB.*chiB) + sum(wtsC.*chiC);
chiA = chiA/chibar;
chiB = chiB/chibar;
chiC = chiC/chibar;
disp(sprintf('  chibar = %12e',chibar))

% shortcuts for radiation spectrum at ionization thresholds
chin1 = chifun(nu0_HI);
chin2 = chifun(nu0_HeI);
chin3 = chifun(nu0_HeII);


% loop over space
[nx,ny,nz] = size(Ef);
for k=1:nz, for j=1:ny, for i=1:nx
   % extract radiation values at this location
   Efval = Ef(i,j,k);
   E1val = E1(i,j,k);
   E2val = E2(i,j,k);
   E3val = E3(i,j,k);
   E1rat = E1val/Efval*chibar/chin1;
   E2rat = E2val/Efval*chibar/chin2;
   E3rat = E3val/Efval*chibar/chin3;

   % get E values at prescribed frequencies
   fEA = Efval.*chiA.*(E1rat).^S1A;
   fEB = Efval.*chiB.*(E1rat).^S1B.*(E2rat).^S2B;
   fEC = Efval.*chiC.*(E1rat).^S1C.*(E2rat).^S2C.*(E3rat).^S3C;

   % compute the integrated grey radiation
   Eg(i,j,k) = sum(wtsA.*fEA) + sum(wtsB.*fEB) + sum(wtsC.*fEC);
end, end, end


% end of script
