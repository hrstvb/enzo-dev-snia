clear

% initialize parameters
nnus = 21;  % must be odd
%nnus = 999;  % must be odd
disp(sprintf('\n  nnus = %i',nnus))

% set input radiation values
Efval = 2.714226485423713e-13;
E1val = 3.113596012320412e-34*Efval;
E2val = 4.329726899046675e-34*Efval;
E3val = 1.391988129146205e-34*Efval;


E1rat = E1val/Efval;
E2rat = E2val/Efval;
E3rat = E3val/Efval;


% set constants
hp = 6.6260693e-27;            % Planck's constant (ergs*s)
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
nu0_HI   = 13.6*ev2erg/hp;     % ionization threshold of HI (hz)
nu0_HeI  = 24.6*ev2erg/hp;     % ionization threshold of HeI (hz)
nu0_HeII = 54.4*ev2erg/hp;     % ionization threshold of HeII (hz)
Llimit = 0.1;                  % lower limit of int (shift away from 0)
Ulimit = 1 - eps^(0.5);        % upper limit of int (shift away from 1)


% set frequencies for evaluation
nbins = (nnus-1)/2;
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

disp(sprintf('  s*n* = %g %g %g %g %g %g\n',s1n1,s1n2,s1n3,s2n2,s2n3,s3n3))

%    species cross sections across intervals
fHIA   = sigHI(nusA);
fHIB   = sigHI(nusB);
fHIC   = sigHI(nusC);
fHeIB  = sigHeI(nusB);
fHeIC  = sigHeI(nusC);
fHeIIC = sigHeII(nusC);

disp(sprintf('  fH* = %g %g %g %g %g %g\n',...
   sum(fHIA),sum(fHIB),sum(fHIC),sum(fHeIB),sum(fHeIC),sum(fHeIIC)))

%    radiation spectrum across intervals
chiA = chifun(nusA);
chiB = chifun(nusB);
chiC = chifun(nusC);

disp(sprintf('  chi* = %g %g %g\n',sum(chiA),sum(chiB),sum(chiC)))

%    exponents for radiation terms
S3C = fHeIIC/s3n3;
S2C = fHeIC/s2n2 - s2n3*S3C;
S2B = fHeIB/s2n2;
S1C = fHIC/s1n1  - s1n3*S3C - s1n2*S2C;
S1B = fHIB/s1n1  - s1n2*S2B;
S1A = fHIA/s1n1;

disp(sprintf('  S* = %g %g %g %g %g %g\n',...
   sum(S3C),sum(S2B),sum(S1A),sum(S2C),sum(S1B),sum(S1C)))

% compute the quadrature weights across the intervals
wtsA = zeros(nnus,1);
wtsB = zeros(nnus,1);
wtsC = zeros(nnus,1);
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

disp(sprintf('  wts* = %g %g %g\n',sum(wtsA),sum(wtsB),sum(wtsC)))

% compute the integrated radiation spectrum
chibar = sum(wtsA.*chiA) + sum(wtsB.*chiB) + sum(wtsC.*chiC);

disp(sprintf('  chibar = %g',chibar))

%    shortcuts for radiation spectrum at ionization thresholds
chin1 = chifun(nu0_HI)/chibar;
chin2 = chifun(nu0_HeI)/chibar;
chin3 = chifun(nu0_HeII)/chibar;

disp(sprintf('  chin* = %g %g %g',chin1,chin2,chin3))

% get E values at prescribed frequencies
fEA = Efval.*chiA./chibar.*(E1rat./chin1).^S1A;
fEB = Efval.*chiB./chibar.*(E1rat./chin1).^S1B.*(E2rat./chin2).^S2B;
fEC = Efval.*chiC./chibar.*(E1rat./chin1).^S1C.*(E2rat./chin2).^S2C.*(E3rat./chin3).^S3C;

% compute the integrated photo-ionization terms
piHI   = sum(wtsA.*c./hp.*fHIA.*fEA./nusA) ...
       + sum(wtsB.*c./hp.*fHIB.*fEB./nusB) ...
       + sum(wtsC.*c./hp.*fHIC.*fEC./nusC);
piHeI  = sum(wtsB.*c./hp.*fHeIB.*fEB./nusB) ...
       + sum(wtsC.*c./hp.*fHeIC.*fEC./nusC);
piHeII = sum(wtsC.*c./hp.*fHeIIC.*fEC./nusC);

% compute the integrated photo-heating coefficients
GHI   = sum(wtsA.*c.*fHIA.*fEA.*(1-nu0_HI./nusA)) ...
      + sum(wtsB.*c.*fHIB.*fEB.*(1-nu0_HI./nusB)) ...
      + sum(wtsC.*c.*fHIC.*fEC.*(1-nu0_HI./nusC));
GHeI  = sum(wtsB.*c.*fHeIB.*fEB.*(1-nu0_HeI./nusB)) ...
      + sum(wtsC.*c.*fHeIC.*fEC.*(1-nu0_HeI./nusC));
GHeII = sum(wtsC.*c.*fHeIIC.*fEC.*(1-nu0_HeII./nusC));

% $$$ disp(sprintf('  E* = %g %g %g %g',Efval,E1rat,E2rat,E3rat))
% $$$ disp(sprintf('  fE* = %g %g %g',sum(fEA),sum(fEB),sum(fEC)))
% $$$ disp(sprintf('  pi* = %g %g %g',piHI,piHeI,piHeII))
% $$$ disp(sprintf('  G* = %g %g %g',GHI,GHeI,GHeII))


% compute these rates a bit differently
piHIfunA = GamfunHI(nusA);
piHIfunB = GamfunHI(nusB);
piHIfunC = GamfunHI(nusC);
piHeIfunB = GamfunHeI(nusB);
piHeIfunC = GamfunHeI(nusC);
piHeIIfunC = GamfunHeII(nusC);
GHIfunA = GfunHI(nusA);
GHIfunB = GfunHI(nusB);
GHIfunC = GfunHI(nusC);
GHeIfunB = GfunHeI(nusB);
GHeIfunC = GfunHeI(nusC);
GHeIIfunC = GfunHeII(nusC);
piHI2   = sum(wtsA.*piHIfunA)  + sum(wtsB.*piHIfunB) + sum(wtsC.*piHIfunC);
piHeI2  = sum(wtsB.*piHeIfunB) + sum(wtsC.*piHeIfunC);
piHeII2 = sum(wtsC.*piHeIIfunC);
GHI2    = sum(wtsA.*GHIfunA)  + sum(wtsB.*GHIfunB) + sum(wtsC.*GHIfunC);
GHeI2   = sum(wtsB.*GHeIfunB) + sum(wtsC.*GHeIfunC);
GHeII2  = sum(wtsC.*GHeIIfunC);


% compute errors in photo-ionization and photo-heating rates
disp('   ')
disp(sprintf(' piHI:   method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    piHI,piHI2,abs(piHI-piHI2)/abs(piHI)))
disp(sprintf(' piHeI:  method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    piHeI,piHeI2,abs(piHeI-piHeI2)/abs(piHeI)))
disp(sprintf(' piHeII: method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    piHeII,piHeII2,abs(piHeII-piHeII2)/abs(piHeII)))
disp(sprintf(' GHI:    method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    GHI,GHI2,abs(GHI-GHI2)/abs(GHI)))
disp(sprintf(' GHeI:   method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    GHeI,GHeI2,abs(GHeI-GHeI2)/abs(GHeI)))
disp(sprintf(' GHeII:  method1 = %17.10e,  method2 = %17.10e,  rdiff = %9.2e',...
    GHeII,GHeII2,abs(GHeII-GHeII2)/abs(GHeII)))
disp('   ')


% end of script
