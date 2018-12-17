%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 backprojection(proj,param,projNumber)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function makes the simple backprojection of 2D images, in order 
%     to reconstruct a certain number of slices.
%  
%     The function calculates for each voxel, which detector pixel is
%     associated with that voxel in the specific projection. That is done 
%     for all angles specified.
%     
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - proj = 2D projection images 
%     - param = Parameter of all geometry
%     - projNumber = Vector with projections numbers to be processed
% 
%     OUTPUT:
% 
%     - data3d = Volume data.
% 
%     Reference: Patent US5872828
% 
%     ---------------------------------------------------------------------
%     Copyright (C) <2018>  <Rodrigo de Barros Vimieiro>
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
% =========================================================================
%% Backprojection Code
function data3d = backprojection(proj,param,projNumber)

global gpuprocess

% Get parameters from struct
DSR = param.DSR;
DDR = param.DDR;
tubeAngle = deg2rad(param.tubeDeg);
numVPixels = param.nv;
numUPixels = param.nu;
numSlices = param.nz;
numProjs = param.nProj;
zCoords = param.zs;     % Z object coordinates

% Test if there's specific angles
if(isempty(projNumber))
    projNumber = 1:numProjs;
else
    if(max(projNumber(:)) <= numProjs)
        numProjs = size(projNumber,2);
    else
        error('Projection number exceeds the maximum for the equipment.')
    end
end


% Stack of reconstructed slices
if(gpuprocess == 1)
    data3d = zeros(param.ny, param.nx, param.nz,'single','gpuArray');
else
    data3d = zeros(param.ny, param.nx, param.nz,'single');
end

% Object Coordinate sytem in (mm)
[xCoord,yCoord] = meshgrid(param.xs,param.ys);
if(gpuprocess == 1)
    xCoord = gpuArray(xCoord);
    yCoord = gpuArray(yCoord);
end

% For each projection
for p=1:numProjs
    
    % Get specific projection number
    projN = projNumber(p);
    
    % Get specific angle for the backprojection
    teta = tubeAngle(projN);
    
    % Get specif projection from stack
    proj_tmp = proj(:,:,p);
    
    % For each slice
    for nz=1:numSlices

        % Calculates a relation of voxel and detector coordinates
        pvCoord = yCoord+((zCoords(nz).*((DSR.*sin(teta))+ yCoord))./...
                        ((DSR.*cos(teta))+ DDR -zCoords(nz)));
                
        puCoord = (xCoord.*((DSR.*cos(teta))+DDR))./ ...
                 ((DSR.*cos(teta))+DDR-zCoords(nz));
                
        % Coordinate in Pixel of detector plane axis origin
        u0 = numUPixels;
        v0 = numVPixels/2;
        
        % Represent detector plane axis (Xi,Yi) in the image plane axis (i,j) (covert mm to Pixels)
        puCoord = -(puCoord./param.du) + u0;
        pvCoord =  (pvCoord./param.dv) + v0;      
        
        % Interpolation of the pixel coordinates of the "Projection" at the
        % calculated pixel coordinates for the slices
        slice_tmp = interp2(proj_tmp,puCoord,pvCoord,'linear', 0);       
        
        % Acumulate current slice with slice from previus projection
        data3d(:,:,nz) = data3d(:,:,nz) + slice_tmp;
        
    end % Loop end slices 
    
end % Loop end Projections
if(gpuprocess == 1)
    data3d = gather(data3d);
end
data3d(:) = data3d(:)./ numProjs;
end