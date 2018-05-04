%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
---------------------------------------------------------------------------
                backprojection(proj,param)
---------------------------------------------------------------------------
    DESCRIPTION:
    This function makes the simple backprojection of 2D images, in order to
    reconstruct a certain number of slices.
 
    The function calculates for each voxel, which detector pixel is
    associated with that voxel in the specific projection. That is done for
    all angles specified.
    
    The geometry is for DBT with half cone-beam. All parameters are set in 
    "ParameterSettings" code. 
 
    INPUT:

    - proj = 2D projection images 
    - param = Parameter of all geometry

    OUTPUT:

    - data3d = Volume data.

    Reference: Patent US5872828

    -----------------------------------------------------------------------
    Copyright (C) <2018>  <Rodrigo de Barros Vimieiro>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
% =========================================================================
%% Backprojection Code
function data3d = backprojection(proj,param)

global gpuprocess

% Stack of reconstructed slices
if(gpuprocess == 0)
    data3d = zeros(param.ny, param.nx, param.nz,'single');
else
    data3d = zeros(param.ny, param.nx, param.nz,'single','gpuArray');
end

% Object Coordinate sytem in (mm)
[xCoord,yCoord] = meshgrid(param.xs,param.ys);
if(gpuprocess == 1)
    xCoord = gpuArray(xCoord);
    yCoord = gpuArray(yCoord);
end

% Get parameters from struct
DSR = param.DSR;
DDR = param.DDR;
tubeAngle = deg2rad(param.tubeDeg);
numVPixels = param.nv;
numUPixels = param.nu;
numProjs = param.nProj;
zCoords = param.zs;     % Z object coordinates

sliceRange = param.sliceRange;
iROI = param.iROI;
jROI = param.jROI;

% For each projection
for p=1:numProjs
    
    % Get specif angle for the backprojection
    teta = tubeAngle(p);
    
    % Get specif projection from stack
    proj_tmp = proj(:,:,p);
    
    % For each slice
    for nz=sliceRange(1):sliceRange(end)

        % Calculates a relation of voxel and detector coordinates
        pvCoord = yCoord(iROI,jROI)+((zCoords(nz).*((DSR.*sin(teta))+ yCoord(iROI,jROI)))./...
                        ((DSR.*cos(teta))+ DDR -zCoords(nz)));
                
        puCoord = (xCoord(iROI,jROI).*((DSR.*cos(teta))+DDR))./ ...
                 ((DSR.*cos(teta))+DDR-zCoords(nz));
                
        % Coordinate in Pixel of detector plane axis origin
        u0 = numUPixels;
        v0 = numVPixels/2;
        
        % Represent detector plane axis (Xi,Yi) in the image plane axis (i,j) (covert mm to Pixels)
        puCoord = -(puCoord./param.du) + u0;
        pvCoord =  (pvCoord./param.dv) + v0;      
        
        % Interpolation of the pixel coordinates of the "Projection" at the
        % calculated pixel coordinates for the slices
        slice_tmp = interp2(proj_tmp,puCoord,pvCoord,'linear',0);       
        
        % Acumulate current slice with slice from previus projection
        data3d(iROI,jROI,nz) = data3d(iROI,jROI,nz) + slice_tmp;
        
    end % Loop end slices 
    
end % Loop end Projections
if(gpuprocess == 1)
    data3d = gather(data3d);
end
data3d(:) = data3d(:)./ numProjs;
end