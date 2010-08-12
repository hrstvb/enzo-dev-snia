function [A,stencil,nx,ny,nz] = read_HypreStructMat(matfile)

% Usage: [A,stencil,nx,ny,nz] = read_HypreStructMat('matfile')
% 
% reads in the file 'matfile' and determines both the stencil components on
% this variable, as well as the data array that determines the associated
% matrix values.  The matrix data is then converted to a standard matrix
% format, in which unknowns are ordered with the x-dimension first, then y,
% then z; e.g. a 2x2x2 problem will index the matrix entries as 
%       (1,1,1) -> 1
%       (2,1,1) -> 2
%       (1,2,1) -> 3
%       (2,2,1) -> 4
%       (1,1,2) -> 5
%       (2,1,2) -> 6
%       (1,2,2) -> 7
%       (2,2,2) -> 8
%
% Note: this routine assumes (for now) that the matrix is neither symmetric
% nor constant-coefficient.
%
% Outputs:
%       A - sparse matrix containing the specified HYPRE Struct matrix
%       stencil - matrix containing the stencil structure
%       nx - x-dimensional size of grid
%       ny - y-dimensional size of grid
%       nz - z-dimensional size of grid
%
% Daniel R. Reynolds
% reynolds@smu.edu

% open file
fid = fopen(matfile,'r');

% get the grid dimensions
line = fgets(fid);  % should be 'StructMatrix'
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Symmetric: 0'
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'ConstantCoefficient: 0'
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Grid:'
ndim = fscanf(fid,'%i',1);  
line = fgets(fid);  % rest of line
if ((ndim < 1) | (ndim > 3))
   disp(sprintf('error: illegal ndim = %i',ndim));
   return
end
line = fgets(fid);  % should contain the number 1 (don't know what it's for)
extents = fscanf(fid,'0:  (%i, %i, %i) x (%i, %i, %i)',6);  
line = fgets(fid);  % rest of line
nx = extents(4)-extents(1)+1;
ny = extents(5)-extents(2)+1;
nz = extents(6)-extents(3)+1;
if (ndim < 3) 
   if (nz ~= 1)
      disp(sprintf('error, ndim = %i and nz = %i do not match',ndim,nz))
      return
   end
end
if (ndim == 1) 
   if (ny ~= 1)
      disp(sprintf('error, ndim = 1 and ny = %i do not match',ny))
      return
   end
end
      
% get the stencil size
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Stencil:'
stsize = fscanf(fid,'%i',1);
line = fgets(fid);  % rest of line

% get the stencil couplings
stencil = zeros(stsize,ndim);
for i=1:stsize
   strow = fscanf(fid,'%i: %i %i %i',4);
   line = fgets(fid);  % rest of line
   stencil(i,1) = strow(2);
   stencil(i,2) = strow(3);
   stencil(i,3) = strow(4);
end


% read the data and put into sparse matrix
irows = zeros(nx*ny*nz*stsize,1);
icols = zeros(nx*ny*nz*stsize,1);
vals  = zeros(nx*ny*nz*stsize,1);
line = fgets(fid);  % should be blank
line = fgets(fid);  % should be 'Data:'
for i=1:nx*ny*nz*stsize
   matrow = fscanf(fid,'%i: (%i, %i, %i; %i) %g',6);
   line = fgets(fid);  % rest of line
   ix = matrow(2)-extents(1)+1;  % shift to start at 1
   iy = matrow(3)-extents(2)+1;  % shift to start at 1
   iz = matrow(4)-extents(3)+1;  % shift to start at 1
   is = matrow(5)+1;
   jx = ix + stencil(is,1);
   jy = iy + stencil(is,2);
   jz = iz + stencil(is,3);
   if (jx == 0), jx = nx; end
   if (jy == 0), jy = ny; end
   if (jz == 0), jz = nz; end
   if (jx == nx+1), jx = 1; end
   if (jy == ny+1), jy = 1; end
   if (jz == nz+1), jz = 1; end
   
   irows(i) = ((iz-1)*ny + (iy-1))*nx + ix;
   icols(i) = ((jz-1)*ny + (jy-1))*nx + jx;
   vals(i) = matrow(6);
end
A = sparse(irows,icols,vals,nx*ny*nz,nx*ny*nz,nx*ny*nz*stsize);


% close hypre matrix file
fclose(fid);

% if matrix was symmetric, propagate to full matrix (not implemented yet)

% if matrix was constant-coefficient, propagate to full matrix 
% (not implemented yet)


% end of routine