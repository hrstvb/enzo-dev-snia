% load the radiation history file
% the first row should contain
%   lUnit, mUnit, tUnit, E1scale, E2scale, E3scale
% the subsequent rows should contain
%   time, Ef, E1, E2, E3, piHI, GHI
% where the E* are all in normalized units
load radhist.txt
lUnit = radhist(1,1);
mUnit = radhist(1,2);
tUnit = radhist(1,3);
E1scale = radhist(1,4);
E2scale = radhist(1,5);
E3scale = radhist(1,6);
Rhist = radhist(2:end,:);  clear radhist

% set derived unit quantities
dUnit = mUnit/lUnit/lUnit/lUnit;
rUnit = dUnit*lUnit*lUnit/tUnit/tUnit;
E1Unit = rUnit*E1scale;
E2Unit = rUnit*E2scale;
E3Unit = rUnit*E3scale;

% set some constants
hp = 6.6260693e-27;            % Planck's constant (ergs*s)
ev2erg = 1.60217653e-12;       % conversion constant from eV to ergs
c = 2.99792458e10;             % speed of light (cm/s)
nu0_HI   = 13.6*ev2erg/hp;     % ionization threshold of HI (hz)
nu0_HeI  = 24.6*ev2erg/hp;     % ionization threshold of HeI (hz)
nu0_HeII = 54.4*ev2erg/hp;     % ionization threshold of HeII (hz)

% set frequency-related arrays (independent of E)
nus = linspace(nu0_HI,2*nu0_HeII,10000);
chis = chifun(nus);
chibar = quadl('chifun',nu0_HI,10*nu0_HeII,1e-8);

% run through time steps
[Nt, tmp] = size(Rhist);
for tstep = 1:Nt
   
   % get the current values
   t    = Rhist(tstep,1);
   Ef   = Rhist(tstep,2);
   E1   = Rhist(tstep,3);
   E2   = Rhist(tstep,4);
   E3   = Rhist(tstep,5);
   piHI = Rhist(tstep,6);
   GHI  = Rhist(tstep,7);
   
   % get E(nu), Eot
   Enu = Erad2(nus,Ef,E1,E2,E3,rUnit);
   Eot = Ef*rUnit*chis/chibar;

   % get E and Ef values at frequency thresholds
   E1 = E1*E1Unit;
   E2 = E2*E2Unit;
   E3 = E3*E3Unit;
   Ef1 = chifun(nu0_HI)/chibar*Ef*rUnit;
   Ef2 = chifun(nu0_HeI)/chibar*Ef*rUnit;
   Ef3 = chifun(nu0_HeII)/chibar*Ef*rUnit;

   % create plot snapshot
   figure(1)
   loglog(nus,Enu,'r-',nu0_HI,E1,'ro',nu0_HeI,E2,'ro',nu0_HeII,E3,'ro',nus,Eot,'b-',nu0_HI,Ef1,'bo',nu0_HeI,Ef2,'bo',nu0_HeII,Ef3,'bo')
   title(sprintf('E(%.2e):   piHI = %.2e,   GHI = %.2e',t,piHI,GHI))
   pause(0.1)
   
end