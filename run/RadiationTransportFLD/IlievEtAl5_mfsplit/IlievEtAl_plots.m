function [] = IlievEtAl_plots(te,endonly)
%
% Usage:
%   [] = IlievEtAl_plots(te,endonly)
%
%   Inputs:
%     te = final time dump
%     endonly = parameter specifying to generate intermediate results (0)
%               or only the final results (1)
%
% This function post-processes the Enzo outputs to plot the various pictures
% from Test 5 of the 'Tests for Cosmological Radiative Code Comparison
% Project', by Iliev et al.
%
% Daniel R. Reynolds, 10/31/2008
%

% set some constants
Ngammadot = 5.0e48;    % ionization source strength [photons/sec]
aHII = 2.59e-13;       % recombination rate coefficient
mp = 1.67262171e-24;   % proton mass [g]
kb = 1.3806504e-16;    % Boltzmann constant
Myr = 3.15576e13;      % duration of a Megayear [sec]
nH = 1e-3;             % input hydrogen number density [cm^(-3)]
trec = 1/(aHII*nH);    % recombination time [sec]
tS = 0;

% initialize time-history outputs
%    col 1: time (t)
%    col 2: computed i-front radius
%    col 3: predicted i-front radius (rI)
%    col 4: stromgren sphere radius (rs)
%    col 5: final predicted stromgren sphere radius (rf)
rdata = zeros(te+1,5);
HIfraction = zeros(te+1,2);

% generate time-dependent fields
for tstep=0:te

   % directory and file name mangling
   s = sprintf('%i',tstep);
   if (tstep < 10)
      dir = [ 'DD000', s ];
      paramfile  = [ dir, '/data000' , s ];
   elseif (tstep < 100)
      dir = [ 'DD00', s ];
      paramfile  = [ dir, '/data00' , s ];
   elseif (tstep < 1000)
      dir = [ 'DD0', s ];
      paramfile  = [ dir, '/data0' , s ];
   elseif (tstep < 10000)
      dir = [ 'DD', s ];
      paramfile  = [ dir, '/data' , s ];
   else
      disp(sprintf('error: te = %i outside of range',te));
      break;
   end
   paramfile2 = [ paramfile, '.rtmodule' ];
   hdfile  = [ paramfile , '.cpu0000' ];
   
   % open parameter file for input
   fin = fopen(paramfile,'r');
   if fin < 0
      error([ 'Could not open ',paramfile,' for input']);
   end

   % examine parameter file for relevant fields
   %    left domain boundary
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'DomainLeftEdge'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 leftedge = str2num(tmp); frewind(fin); break
      end
   end
   xL = leftedge(1);  yL = leftedge(2);  zL = leftedge(3);
   %    right domain boundary
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'DomainRightEdge'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 rightedge = str2num(tmp); frewind(fin); break
      end
   end
   xR = rightedge(1);  yR = rightedge(2);  zR = rightedge(3);
   %    current time
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'InitialTime'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 t = str2num(tmp); frewind(fin); break
      end
   end
   %    gamma (ideal gas law parameter)
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'Gamma'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 gamma = str2num(tmp); frewind(fin); break
      end
   end
   %    density units
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'DensityUnits'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 dUnits = str2num(tmp); frewind(fin); break
      end
   end
   %    length units
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'MassUnits'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 mUnits = str2num(tmp); frewind(fin); break
      end
   end
   %    time units
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'TimeUnits'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 tUnits = str2num(tmp); frewind(fin); break
      end
   end
   lUnits = (mUnits / dUnits)^(1/3);
   fclose(fin);

   % open parameter2 file for input
   fin2 = fopen(paramfile2,'r');
   if fin2 < 0
      error([ 'Could not open ',paramfile2,' for input']);
   end

   % examine parameter file for relevant fields
   %    E1Units scaling factor
   for i=1:1000
      buffer = fgetl(fin2);
      if (findstr(buffer,'E1Units'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 E1Unit = str2num(tmp); frewind(fin2); break
      end
   end
   %    E2Units scaling factor
   for i=1:1000
      buffer = fgetl(fin2);
      if (findstr(buffer,'E2Units'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 E2Unit = str2num(tmp); frewind(fin2); break
      end
   end
   %    E3Units scaling factor
   for i=1:1000
      buffer = fgetl(fin2);
      if (findstr(buffer,'E3Units'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 E3Unit = str2num(tmp); frewind(fin2); break
      end
   end
   fclose(fin2);

   % convert to CGS units
   xL = xL * lUnits;  xR = xR * lUnits;  
   yL = yL * lUnits;  yR = yR * lUnits;  
   zL = zL * lUnits;  zR = zR * lUnits;  
   t  =  t * tUnits;
   
   % load hdf5 file, and convert data to CGS values
   rho = double(hdf5read(hdfile,'Grid00000001/Density')); 
   rho = rho*dUnits;
   HI = double(hdf5read(hdfile,'Grid00000001/HI_Density'));
   HI = HI*dUnits + mp*1e-10*ones(size(HI));
   HII = double(hdf5read(hdfile,'Grid00000001/HII_Density'));
   HII = HII*dUnits + mp*1e-10*ones(size(HII));;
   Ef = double(hdf5read(hdfile,'Grid00000001/FS_Radiation_Energy'));
   Ef = Ef*lUnits/tUnits*lUnits/tUnits*dUnits;
   E1 = double(hdf5read(hdfile,'Grid00000001/Radiation_Energy1'));
   E1 = E1*lUnits/tUnits*lUnits/tUnits*dUnits*E1Unit;
   E2 = double(hdf5read(hdfile,'Grid00000001/Radiation_Energy2'));
   E2 = E2*lUnits/tUnits*lUnits/tUnits*dUnits*E2Unit;
   E3 = double(hdf5read(hdfile,'Grid00000001/Radiation_Energy3'));
   E3 = E3*lUnits/tUnits*lUnits/tUnits*dUnits*E3Unit;
   energy = double(hdf5read(hdfile,'Grid00000001/Total_Energy'));
   energy = energy*lUnits/tUnits*lUnits/tUnits;
   vx = double(hdf5read(hdfile,'Grid00000001/x-velocity'));
   vx = vx*lUnits/tUnits;
   vy = double(hdf5read(hdfile,'Grid00000001/y-velocity'));
   vy = vy*lUnits/tUnits;
   vz = double(hdf5read(hdfile,'Grid00000001/z-velocity'));
   vz = vz*lUnits/tUnits;

   % compute grey radiation field over domain
   Eg = Egrey(Ef,E1,E2,E3);
   
   % get grid dimensions, compute volume element
   [Nx,Ny,Nz] = size(HI);
   dV = (xR-xL)*(yR-yL)*(zR-zL)/Nx/Ny/Nz;
   
   % integrate HII fraction over domain
   HIintegral = sum(sum(sum(HI./rho)))*dV;

   % compute volume of ionized region (50% or more HII)
   %   mirror this solution in each octant to get spherical region
   HIIvolume1 = sum(sum(sum(round(HII./rho))))*dV*8;
   HIIvolume2 = sum(sum(sum(HII./rho)))*dV*8;
%   HIIvolume = (HIIvolume1 + HIIvolume2)/2;
   HIIvolume = HIIvolume2;
   
   % compute radius of ionized region assuming spherical ionization region
   %    i.e. volume = 4/3*pi*r^3  ->  r = (3*volume/4/pi)^(1/3)
   radius = (3*HIIvolume/4/pi)^(1/3);

   % fill radius/time array [time, total neutral fraction]
   HIfraction(tstep+1,1) = t/trec;
   HIfraction(tstep+1,2) = HIintegral/xR/yR/zR;

   % compute Stromgren radius, predicted HII region radius
   %   pressure/density/temperature behind I-front
   Icells = radius/xR*Nx;      % number of cells behind front
   pI = 0;  dI = 0;  Ti = 0;  Te = 0;  ncells = 0;
   for k=1:Nz,  for j=1:Ny,  for i=1:Nx
      if (sqrt(i*i+j*j+k*k)<=Icells) 
%	 pI = pI + energy(i,j,k)*rho(i,j,k)*(gamma-1) ...
%	         - 0.5*(vx(i,j,k)^2 + vy(i,j,k)^2 + vz(i,j,k)^2);
	 pI = pI + energy(i,j,k)*rho(i,j,k)*(gamma-1);
	 dI = dI + rho(i,j,k);
	 Ti = Ti + energy(i,j,k)*rho(i,j,k)*(gamma-1)/(2*rho(i,j,k) - HI(i,j,k))*mp/kb;
	 ncells = ncells + 1;
      else
	 Te = Te + energy(i,j,k)*rho(i,j,k)*(gamma-1)/(2*rho(i,j,k) - HI(i,j,k))*mp/kb;
      end
   end, end, end
   Ti = Ti/ncells;                   % average temperature behind i-front
   Te = Te/(Nx*Ny*Nz - ncells);      % average temperature beyond i-front
   rs0 = (3*Ngammadot/4/pi/aHII/nH/nH)^(1/3);
   Ti = 1e4;      % temperature inside HII region
   Te = 1e2;      % temperature outside HII region
   rf  = (2*Ti/Te)^(2/3)*rs0;  % predicted final HII region radius
   cspeed = sqrt(pI/dI);  % average sound speed behind I-front

   % if radius <= rs0, use static solution, otherwise use dynamic (with 't'
   % set as the time since reaching the Stromgren radius.
   % (in both of these, accomodate for the fact that Matlab screws up
   %  fractional powers of negative numbers)
   if (radius < rs0)     % use static solution
      arg = 1-exp(-t/trec);
      if (arg < 0)
	 ranal = -rs0*abs(arg)^(1/3);
      else
	 ranal = rs0*arg^(1/3);
      end
   else     % use dynamic solution
      if (tS == 0), tS = t; end  % store time at which rs0 is reached
      tdiff = t-tS;
      arg = 1 + 7*cspeed*tdiff/4/rs0;
      if (arg < 0) 
	 ranal = -rs0*abs(arg)^(4/7);
      else
	 ranal = rs0*arg^(4/7);
      end
   end
   
   % fill radius/time array 
   rdata(tstep+1,:) = [t/trec, radius, ranal, rs0, rf];
   
   
   % generate continuous plots (if desired)
   Myr0 = 0;  Myr10 = 0;  Myr30 = 0;  Myr100 = 0;  
   Myr200 = 0;  Myr500 = 0;  Myr1000 = 0;
   if (tstep == 0),   Myr0    = 1; t = 0;        MyrStub = '0Myr';    end
   if (tstep == 1),   Myr10   = 1; t = Myr*10;   MyrStub = '10Myr';   end
   if (tstep == 3),   Myr30   = 1; t = Myr*30;   MyrStub = '30Myr';   end
   if (tstep == 10),  Myr100  = 1; t = Myr*100;  MyrStub = '100Myr';  end
   if (tstep == 20),  Myr200  = 1; t = Myr*200;  MyrStub = '200Myr';  end
   if (tstep == 50),  Myr500  = 1; t = Myr*500;  MyrStub = '500Myr';  end
   if (tstep == 100), Myr1000 = 1; t = Myr*1000; MyrStub = '1000Myr'; end
   if ((endonly == 0) || Myr0 || Myr10 || Myr30 || Myr100 || Myr200 || Myr500 || Myr1000)
      
      % i-front radius vs time
      h = figure(1);
      plot(rdata(1:tstep+1,1), rdata(1:tstep+1,2)/rs0,'b-',...
	   rdata(1:tstep+1,1), rdata(1:tstep+1,3)/rs0,'r--',...
	   rdata(1:tstep+1,1), rf/rs0*ones(tstep+1,1),'k:','LineWidth',2);
      title(sprintf('Propagation of HII Region, t = %g Myr',t/Myr),'FontSize',16);
      xlabel('t/t_{rec}','FontSize',14);  ylabel('r_I/r_S','FontSize',14);
      hl = legend('computed','analytical','final',2);  
      set(hl,'FontSize',14);  set(gca,'FontSize',14);
      axis([0 4.5 0 3.0]);  grid on
      if (tstep == te) 
% $$$ 	 saveas(h,'rad_vs_time.fig');
% $$$ 	 saveas(h,'rad_vs_time.eps','epsc');
	 saveas(h,'rad_vs_time.png');
	 save -ascii radiusdata.txt rdata
      end
      
      % i-front velocity vs time
      h = figure(2);
      times = (rdata(1:tstep,1) + rdata(2:tstep+1,1))/2;
      velocity = (rdata(2:tstep+1,2) - rdata(1:tstep,2)) ...
	       ./(rdata(2:tstep+1,1) - rdata(1:tstep,1))./rs0;
      vel_anal = (rdata(2:tstep+1,3) - rdata(1:tstep,3)) ...
	       ./(rdata(2:tstep+1,1) - rdata(1:tstep,1))./rs0;
      for k=2:tstep-1
	 if (vel_anal(k) > (vel_anal(k-1)+vel_anal(k+1))),
	    vel_anal(k) = (vel_anal(k-1)+vel_anal(k+1))/2;
	 end
      end
      plot(times,velocity,'b-',...
	   times,vel_anal,'r--','LineWidth',2);
      title(sprintf('Velocity of HII Region, t=%g Myr',t/Myr),'FontSize',16);
      xlabel('t/t_{rec}','FontSize',14);  ylabel('v_i / (r_S/t_{rec})','FontSize',14);
      hl = legend('computed','analytical');  
      set(hl,'FontSize',14);  set(gca,'FontSize',14);
%      axis([0, 4.5, -0.2, 1.2]);  grid on
      if (tstep == te) 
% $$$ 	 saveas(h,'vel_vs_time.fig');
% $$$ 	 saveas(h,'vel_vs_time.eps','epsc');
	 saveas(h,'vel_vs_time.png');
	 vdata =  [times, velocity, vel_anal];
	 save -ascii velocitydata.txt vdata;  clear vdata;
      end
      
      % plot neutral fraction history
      h = figure(3);
      plot(HIfraction(1:tstep+1,1),HIfraction(1:tstep+1,2),'b-','LineWidth',2);
      title(sprintf('Total Neutral Fraction, t = %g Myr',t/Myr),'FontSize',16);
      xlabel('t/t_{rec}','FontSize',14);
      ylabel('1-x','FontSize',14);
      set(gca,'FontSize',14);  axis([0 4.5 0.5 1.1]);
      grid on
      if (tstep == te) 
% $$$ 	 saveas(h,'frac_vs_time.fig');
% $$$ 	 saveas(h,'frac_vs_time.eps','epsc');
	 saveas(h,'frac_vs_time.png');
	 save -ascii nfracdata.txt HIfraction
      end

      
      % plot HI fraction contour through z=0
      xspan = linspace(0,1,Nx)';
      yspan = linspace(0,1,Ny)';
      h = figure(4);
      surf(yspan,xspan,log10(HI(:,:,1)./rho(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log HI fraction, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'HIcontour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'HIcontour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'HIcontour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot FS radiation density contour through z=0
      h = figure(5);
      surf(yspan,xspan,log10(Ef(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log FS radiation density, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'Efcontour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'Efcontour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'Efcontour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot radiation density 1 contour through z=0
      h = figure(6);
      surf(yspan,xspan,log10(E1(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log radiation 1 density, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'E1contour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E1contour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E1contour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot radiation density 2 contour through z=0
      h = figure(7);
      surf(yspan,xspan,log10(E2(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log radiation 2 density, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'E2contour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E2contour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E2contour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot radiation density 3 contour through z=0
      h = figure(8);
      surf(yspan,xspan,log10(E3(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log radiation 3 density, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'E3contour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E3contour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E3contour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot grey radiation density contour through z=0
      h = figure(9);
      surf(yspan,xspan,log10(Eg(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log grey radiation density, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'Egcontour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'Egcontour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'Egcontour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot pressure contour through z=0
      Press = (gamma-1)*rho.*(energy-0.5.*(vx.^2+vy.^2+vz.^2));
      h = figure(10);
      surf(yspan,xspan,log10(Press(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log Pressure, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'PressContour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'PressContour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'PressContour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot temperature contour through z=0
      Temp = Press./(2*rho - HI)*mp/kb;
      h = figure(11);
      surf(yspan,xspan,log10(Temp(:,:,1)));
      shading interp, colorbar, axis equal, view(0,90)
      title(sprintf('log Temp, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);
% $$$       fname = [ 'TempContour_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'TempContour_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'TempContour_' , MyrStub , '.png' ];
      saveas(h,fname);

      % compute spherically-averaged profiles for HI, HII fractions,
      % temperature, radiation, pressure, kinetic energy, internal energy,
      % total energy, number density 
      Nradii = 100;
      Hradii  = linspace(0,sqrt(3),Nradii);
      rad_idx = zeros(Nx,Ny,Nz);
      for k=1:Nz
	 zloc = (k-1/2)/Nz;
	 for j=1:Ny
	    yloc = (j-1/2)/Ny;
	    for i=1:Nx
	       xloc = (i-1/2)/Nx;
	       rad_idx(i,j,k) = max([1,ceil(sqrt(xloc^2+yloc^2+zloc^2)/sqrt(3)*Nradii)]);
	    end
	 end
      end
      Hcount  = eps*ones(size(Hradii));
      HIprof  = zeros(size(Hradii));
      HIIprof = zeros(size(Hradii));
      Tempprof = zeros(size(Hradii));
      Efprof = zeros(size(Hradii));
      E1prof = zeros(size(Hradii));
      E2prof = zeros(size(Hradii));
      E3prof = zeros(size(Hradii));
      Egprof = zeros(size(Hradii));
      KEProf = zeros(size(Hradii));
      TEProf = zeros(size(Hradii));
      IEProf = zeros(size(Hradii));
      PressProf = zeros(size(Hradii));
      nHProf = zeros(size(Hradii));
      rsquared = linspace(1,Nradii,Nradii).^2;
      for k=1:Nz
	 for j=1:Ny
	    for i=1:Nx
	       idx = rad_idx(i,j,k);
	       HIprof(idx) = HIprof(idx) + HI(i,j,k)/rho(i,j,k);
	       HIIprof(idx) = HIIprof(idx) + HII(i,j,k)/rho(i,j,k);
	       Tempprof(idx) = Tempprof(idx) + Temp(i,j,k);
	       Efprof(idx) = Efprof(idx) + Ef(i,j,k);
	       E1prof(idx) = E1prof(idx) + E1(i,j,k);
	       E2prof(idx) = E2prof(idx) + E2(i,j,k);
	       E3prof(idx) = E3prof(idx) + E3(i,j,k);
	       Egprof(idx) = Egprof(idx) + Eg(i,j,k);
	       PressProf(idx) = PressProf(idx) + Press(i,j,k);
	       KEProf(idx) = KEProf(idx) + 0.5*(vx(i,j,k)^2 + vy(i,j,k)^2 + vz(i,j,k)^2);
	       TEProf(idx) = TEProf(idx) + energy(i,j,k);
	       IEProf(idx) = IEProf(idx) + energy(i,j,k) - 0.5*(vx(i,j,k)^2 + vy(i,j,k)^2 + vz(i,j,k)^2);
	       nHProf(idx) = nHProf(idx) + rho(i,j,k)/mp;
	       Hcount(idx) = Hcount(idx) + 1;
	    end
	 end
      end
      HIprof = HIprof./Hcount;
      HIIprof = HIIprof./Hcount;
      Tempprof = Tempprof./Hcount;
      Efprof = Efprof./Hcount;
      E1prof = E1prof./Hcount;
      E2prof = E2prof./Hcount;
      E3prof = E3prof./Hcount;
      Egprof = Egprof./Hcount;
      PressProf = PressProf./Hcount;
      KEProf = KEProf./Hcount;
      TEProf = TEProf./Hcount;
      IEProf = IEProf./Hcount;
      nHProf = nHProf./Hcount;
      
      % plot averaged chemistry profiles
      h = figure(12);
%      semilogy(Hradii,HIprof,'b-',Hradii,HIIprof,'r--','LineWidth',2);
      plot(Hradii,log10(HIprof),'b-',Hradii,log10(HIIprof),'r--','LineWidth',2);
      hl = legend('HI','HII');  
      set(hl,'FontSize',14);  set(gca,'FontSize',14);
%      xlabel('r/L_{box}','FontSize',14);  ylabel('x, 1-x','FontSize',14);
      xlabel('r/L_{box}','FontSize',14);  ylabel('log(x), log(1-x)','FontSize',14);
      title(sprintf('HI, HII profiles, t = %g Myr',t/Myr),'FontSize',16);
%      axis([0, 1.2, 1e-7, 10]);  grid on
      axis([0, 1.2, -7, 1]);  grid on
% $$$       fname = [ 'profiles_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'profiles_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'profiles_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged temperature profile
      h = figure(13);
      plot(Hradii,log10(Tempprof),'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('log10(T) [K]','FontSize',14);
      title(sprintf('Temperature profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  axis([0, 1.2, 2.0, 5.0]);  grid on
% $$$       fname = [ 'TempProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'TempProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'TempProfile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged FS Radiation profile
      h = figure(14);
      semilogy(Hradii,Efprof,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('E','FontSize',14);
      title(sprintf('FS Radiation profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'EfProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'EfProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'EfProfile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged Radiation 1 profile
      h = figure(15);
      semilogy(Hradii,E1prof,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('E','FontSize',14);
      title(sprintf('Radiation 1 profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'E1Profile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E1Profile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E1Profile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged Radiation 2 profile
      h = figure(16);
      semilogy(Hradii,E2prof,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('E','FontSize',14);
      title(sprintf('Radiation 2 profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'E2Profile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E2Profile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E2Profile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged Radiation 3 profile
      h = figure(17);
      semilogy(Hradii,E3prof,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('E','FontSize',14);
      title(sprintf('Radiation 3 profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'E3Profile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'E3Profile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'E3Profile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged grey Radiation profile
      h = figure(18);
      semilogy(Hradii,Egprof,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('E','FontSize',14);
      title(sprintf('Grey radiation profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'EgProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'EgProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'EgProfile_' , MyrStub , '.png' ];
      saveas(h,fname);

      % plot averaged Pressure profile
      h = figure(19);
      semilogy(Hradii,PressProf,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('p','FontSize',14);
      title(sprintf('Pressure profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  
      %axis([0 1.2 1e-35 1e-10])
      ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'PressProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'PressProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'PressProfile_' , MyrStub , '.png' ];
      saveas(h,fname);


      % plot averaged Energy profiles
      h = figure(20);
      plot(Hradii,log10(TEProf),'b-',Hradii,log10(KEProf),'r--',Hradii,log10(IEProf),'k:','LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('log10(energy)','FontSize',14);
      hl = legend('Total Energy','Kinetic Energy','Internal Energy',4);
      title(sprintf('Energy profiles, t = %g Myr',t/Myr),'FontSize',16);
      set(hl,'FontSize',14);  set(gca,'FontSize',14);  
      axis([0 1.2 7 14])
      %ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'EnergyProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'EnergyProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'EnergyProfile_' , MyrStub , '.png' ];
      saveas(h,fname);


      % plot averaged number density profile
      h = figure(21);
      semilogy(Hradii,nHProf,'LineWidth',2);
      xlabel('r/L_{box}','FontSize',14);  ylabel('n','FontSize',14);
      title(sprintf('Number density profile, t = %g Myr',t/Myr),'FontSize',16);
      set(gca,'FontSize',14);  
      axis([0 1.2 1e-4 1e-2])
      %ax = axis; ax(1) = 0; ax(2) = 1.2; axis(ax);
      grid on
% $$$       fname = [ 'nProfile_' , MyrStub , '.fig' ];
% $$$       saveas(h,fname);
% $$$       fname = [ 'nProfile_' , MyrStub , '.eps' ];
% $$$       saveas(h,fname,'epsc');
      fname = [ 'nProfile_' , MyrStub , '.png' ];
      saveas(h,fname);

   end   % end plotting stuff
   
   if (endonly == 0) 
      disp(sprintf('finished plots for time step %i, press any key to continue',tstep))
      pause
   end
   
end  % end time steps



% end of file
