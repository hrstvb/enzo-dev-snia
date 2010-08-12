function [] = plot_results()
%
% usage:  plot_results()
%

% loop over time steps
tend = 10;
for it = 0:tend

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
   rUnits = dUnits*(lUnits/tUnits)^2;
   leftedge = leftedge*lUnits;
   rightedge = rightedge*lUnits;
   t = t*tUnits;

   % set hdf5 name
   hfile = [ dir, '/', stub, '.cpu0000' ];
   gridname = [ 'Grid00000001' ];
   
   % read radiation field and re-scale
   name = [ gridname, '/FS_Radiation' ];
   tmp = double(hdf5read(hfile,name));
   E = tmp*rUnits;
   clear tmp
   
   % plot 2D slice at zlevel = 0
   [Nx,Ny,Nz] = size(E);
   xspan = linspace(leftedge(1),rightedge(1),Nx)'./(rightedge(1)-leftedge(1));
   yspan = linspace(leftedge(2),rightedge(2),Ny)'./(rightedge(2)-leftedge(2));
   figure(1), surf(yspan,xspan,log10(E(:,:,1)));
   shading interp, colorbar, axis equal, view(0,90)
%   figure(1), contour(yspan,xspan,log10(E(:,:,1)));
%   colorbar, axis equal, view(0,90)
   xlabel('y'), ylabel('x')
   title(sprintf('log E_{fs}, t = %g Tf',it/10),'FontSize',14);

   % plot 1D profile
   figure(2), plot(xspan,log10(E(:,1,1)));
   axis([0 1 -30 -20]);
   xlabel('x'), title(sprintf('log E_{fs}, t = %g Tf',it/10),'FontSize',14);

   % plot 1D profile
   figure(3), plot(xspan,E(:,1,1));
   axis([0 1 0 1e-20]);
   xlabel('x'), title(sprintf('E_{fs}, t = %g Tf',it/10),'FontSize',14);

   % plot scaled 1D profile
   figure(4), plot(xspan,E(:,1,1).*xspan.*xspan);
%   axis([0 1 0 5e-22]);
   xlabel('x'), title(sprintf('E_{fs}*r^2, t = %g Tf',it/10),'FontSize',14);

   % clean up
   clear E;

   % save plots from mid-way through run
   if (it == 5) 
      h = figure(1);
      saveas(h,'Econtour.fig');
      saveas(h,'Econtour.eps','epsc');
      saveas(h,'Econtour.png');
   
      h = figure(2);
      saveas(h,'Eprofile.fig');
      saveas(h,'Eprofile.eps','epsc');
      saveas(h,'Eprofile.png');

      h = figure(3);
      saveas(h,'Eprofile2.fig');
      saveas(h,'Eprofile2.eps','epsc');
      saveas(h,'Eprofile2.png');
end
   
   % hold on a second
   pause(1)
   
end  % time step loop

% end of file
