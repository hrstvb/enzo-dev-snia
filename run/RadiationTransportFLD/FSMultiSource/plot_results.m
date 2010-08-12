function [] = plot_results(tend)
%
% usage:  plot_results(tend)
%

% loop over time steps
for it = 0:5:tend

   % set subdirectory to load files from
   if (it < 10) 
      dir = [ 'DD000' , num2str(it) ];
      stub = [ 'data000' , num2str(it) ];
   elseif (it < 100) 
      dir = [ 'DD00' , num2str(it) ];
      stub = [ 'data00' , num2str(it) ];
   elseif (it < 1000)
      dir = [ 'DD0' , num2str(it) ];
      stub = [ 'data0' , num2str(it) ];
   elseif (it < 10000)
      dir = [ 'DD' , num2str(it) ];
      stub = [ 'data' , num2str(it) ];
   else
      error([ 'Time step ',it,' is out of range']);
   end

   % get Units from parameter file
   pfile = [ dir, '/', stub ];
   fin = fopen(pfile,'r');
   if (fin < 0)
      error([ 'Could not open ',pfile,' for input']);
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
   %    right domain boundary
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'DomainRightEdge'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 rightedge = str2num(tmp); frewind(fin); break
      end
   end
   %    time
   for i=1:1000
      buffer = fgetl(fin);
      if (findstr(buffer,'InitialTime'))
	 [text,buffer] = strtok(buffer,'=');
	 [text,tmp] = strtok(buffer);
	 t = str2num(tmp); frewind(fin); break
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

   % set radiation units, box bounds, time
   rUnits = (lUnits/tUnits)*dUnits;
   leftedge = leftedge*lUnits;
   rightedge = rightedge*lUnits;
   t  =  t*tUnits;

   % set hdf5 name
   hfile = [ dir, '/', stub, '.cpu0000' ];
   gridname = [ 'Grid00000001' ];
   
   % read radiation field and re-scale
   name = [ gridname, '/FS_Radiation' ];
   tmp = double(hdf5read(hfile,name));
   E = tmp*rUnits;
   clear tmp
   
   % plot 2D projection
   [Nx,Ny,Nz] = size(E);
   xspan = linspace(leftedge(1),rightedge(1),Nx);
   yspan = linspace(leftedge(2),rightedge(2),Ny);
%   figure(1), contour(yspan,xspan,log10(sum(abs(E),3)));
   figure(1), contour(yspan,xspan,sum(abs(E),3));
   a = axis; a(1:4) = [leftedge(1), rightedge(1), leftedge(2), rightedge(2)];
   colorbar, axis equal
   title(sprintf('log E_{fs}, t = %g sec',t),'FontSize',14);

   % plot 1D projection
%   figure(2), plot(xspan,log10(sum(sum(E,3),2)));
   figure(2), plot(xspan,sum(sum(E,3),2));
   xlabel('y'), title(sprintf('log E_{fs}, t = %g sec',t),'FontSize',14);
   
   clear E;
   pause(0.5)
   
end  % time step loop

% end of file
