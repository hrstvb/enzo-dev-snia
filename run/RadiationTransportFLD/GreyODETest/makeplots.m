function [] = makeplots(te)
%
% Usage:
%   [] = makeplots(te)
%
%   Inputs:
%     te = final time dump
%
% This function post-processes the Enzo outputs to:
%  (a) check that the fields remain spatially homogeneous, and
%  (b) plot the time-histories of each variable
%
% Daniel R. Reynolds
%

% initialize time-history outputs (each row is a snapshot)
% cols: [t, rho, ne, Eg, nHI, nHII, etot]
tdata = zeros(te+1,7);

% generate time-dependent fields
for tstep=0:te

   % directory and file name mangling
   s = sprintf('%i',tstep);
   if (tstep < 10)
      paramfile  = [ 'DD000', s , '/data000' , s ];
   elseif (tstep < 100)
      paramfile  = [ 'DD00', s , '/data00' , s ];
   elseif (tstep < 1000)
      paramfile  = [ 'DD0', s , '/data0' , s ];
   elseif (tstep < 10000)
      paramfile  = [ 'DD', s , '/data' , s ];
   else
      disp(sprintf('error: te = %i outside of range',te));
      break;
   end
   hdfile  = [ paramfile , '.cpu0000' ];
   
   % open parameter file for input
   fin = fopen(paramfile,'r');
   if fin < 0
      error([ 'Could not open ',fname,' for input']);
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
      if (findstr(buffer,'LengthUnits'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 lUnits = str2num(tmp); frewind(fin); break
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
   fclose(fin);

   % convert to CGS units
   xL = xL * lUnits;  xR = xR * lUnits;  
   yL = yL * lUnits;  yR = yR * lUnits;  
   zL = zL * lUnits;  zR = zR * lUnits;  
   t  =  t * tUnits;
   
   % load hdf5 file, and convert data to CGS values
   rho = double(hdf5read(hdfile,'Grid00000001/Density')); 
   ne = double(hdf5read(hdfile,'Grid00000001/Electron_Density')); 
   nHI = double(hdf5read(hdfile,'Grid00000001/HI_Density'));
   nHII = double(hdf5read(hdfile,'Grid00000001/HII_Density'));
   Eg = double(hdf5read(hdfile,'Grid00000001/Grey_Radiation_Energy'));
   etot = double(hdf5read(hdfile,'Grid00000001/Total_Energy'));
   ne = ne./rho + 1e-9;
   nHI = nHI./rho + 1e-9;
   nHII = nHII./rho + 1e-9;
   rho = rho*dUnits;
   Eg = Eg*lUnits/tUnits*lUnits/tUnits*dUnits;
   etot = etot*lUnits/tUnits*lUnits/tUnits;
   
   % store values in tdata array
   % cols: [t, rho, ne, Eg, nHI, nHII, etot]
   tdata(tstep+1,:) = [t, rho(1,1,1), ne(1,1,1), Eg(1,1,1), nHI(1,1,1), ...
	  nHII(1,1,1), etot(1,1,1)];
   
   % get grid dimensions
   [Nx,Ny,Nz] = size(rho);

   % check that the domain is spatially homogeneous
   utol = 1e-8;
   %   density
   umax = max(max(max(rho)));
   uavg = sum(sum(sum(rho)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i density field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

   %   electron density
   umax = max(max(max(ne)));
   uavg = sum(sum(sum(ne)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i electron density field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

   %   radiation 
   umax = max(max(max(Eg)));
   uavg = sum(sum(sum(Eg)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i radiation field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

   %   nHI density
   umax = max(max(max(nHI)));
   uavg = sum(sum(sum(nHI)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i nHI density field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

   %   nHII density
   umax = max(max(max(nHII)));
   uavg = sum(sum(sum(nHII)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i nHII density field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

   %   total energy
   umax = max(max(max(etot)));
   uavg = sum(sum(sum(etot)))/Nx/Ny/Nz;
   if (abs(umax - uavg)/abs(umax) > utol) 
      disp(sprintf('Warning: at tstep %i energy field is inhomogeneous (%g vs %g)',tstep,umax,uavg));
      exit
   end

end  % end time steps

% generate plots

%   density vs time
h = figure(1);
semilogy(tdata(:,1),tdata(:,2));
title('Density evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('rho','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'rho.png');

%   electron density vs time
h = figure(2);
semilogy(tdata(:,1),tdata(:,3));
title('Electron density evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('ne','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'ne.png');
      
%   radiation vs time
h = figure(3);
semilogy(tdata(:,1),tdata(:,4));
title('Radiation evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('Eg','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'Eg.png');
      
%   nHI density vs time
h = figure(4);
semilogy(tdata(:,1),tdata(:,5));
title('HI fraction evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('nHI','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'nHI.png');
      
%   nHII density vs time
h = figure(5);
semilogy(tdata(:,1),tdata(:,6));
title('HII fraction evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('nHII','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'nHII.png');
      
%   energy vs time
h = figure(6);
semilogy(tdata(:,1),tdata(:,7));
title('Total energy evolution','FontSize',14);
xlabel('tstep','FontSize',12);  ylabel('etot','FontSize',12);
set(gca,'FontSize',12); grid on
saveas(h,'etot.png');
      


% end of file
