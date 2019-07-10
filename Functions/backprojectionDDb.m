%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 backprojectionDDb(proj,param,projNumber)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function reconstruct the 3D volume from projections, based on 
%     the Branchless Distance-Driven principle. 
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - proj = 2D projections for each angle 
%     - param = Parameter of all geometry
%     - projNumber = Vector with projections numbers to be processed
%
%     OUTPUT:
% 
%     - data3d = reconstructed volume.
% 
%     Reference: 
%     - Branchless Distance Driven Projection and Backprojection,
%     Samit Basu and Bruno De Man (2006)
%     - GPU Acceleration of Branchless Distance Driven Projection and 
%     Backprojection, Liu et al (2016)
%     - GPU-Based Branchless Distance-Driven Projection and Backprojection,
%     Liu et al (2017)
%     - A GPU Implementation of Distance-Driven Computed Tomography, 
%     Ryan D. Wagner (2017)
%     ---------------------------------------------------------------------
%     Copyright (C) <2019>  <Rodrigo de Barros Vimieiro>
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
%% 3-D Backprojection Branchless Distance Driven Code

function data3d = backprojectionDDb(proj,param,projNumber)

% Map detector and object boudaries
param.us = double((param.nu:-1:0)*param.du);
param.vs = double((-(param.nv)/2:1:(param.nv)/2)*param.dv);
param.xs = double((param.nx:-1:0)*param.dx);
param.ys = double((-(param.ny)/2:1:(param.ny)/2)*param.dy);
param.zs = double((0:1:param.nz-1)*param.dz + param.DAG + (param.dz/2));

% Detector boudaries coordinate sytem in (mm)
[detX,detY] = meshgrid(param.us,param.vs);
detZ = zeros(size(detX),'double');

% Object boudaries coordinate sytem in (mm)
[objX,objY] = meshgrid(param.xs,param.ys);
objZ = double(param.zs);     % Z object coordinates

% X-ray tube initial position
tubeX = double(0);
tubeY = double(0);
tubeZ = double(param.DSD);

% Iso-center position
isoY = double(0);
isoZ = double(param.DDR);

% Projection angles
tubeAngle = double(deg2rad(param.tubeDeg));
detAngle = double(deg2rad(param.detectorDeg));

% Number of detector and voxels for each direction
nDetX = double(param.nu);
nDetY = double(param.nv);
nPixX = double(param.nx);
nPixY = double(param.ny);

deltaDetX = double(param.du); 
deltaDetY = double(param.dv); 
deltaObjX = double(param.dx); 
deltaObjY = double(param.dy);

nSlices = double(param.nz);
nProjs = double(param.nProj);

% Voxel size on Z
pixsZ = double(param.dz); 

% Test if there's specific angles
if(isempty(projNumber))
    projNumber = 1:nProjs;
else
    if(max(projNumber(:)) <= nProjs)
        nProjs = size(projNumber,2);
    else
        error('Projection number exceeds the maximum for the equipment.')
    end
end

% Stack of reconstructed slices
data3d = zeros(param.ny, param.nx, param.nz,'double');

% Integration of each 2D projection
% (S.1. Integration - Liu et al (2017))
projI = integrate2D(proj);

% For each projection
for p=1:nProjs
       
    % Get specific projection number
    projN = projNumber(p);
    
    % Get specif tube angle for the projection
    theta = tubeAngle(projN);
    
    % Get specif detector angle for the projection
    phi = detAngle(projN);
    
    % Tubre rotation
    rtubeY = ((tubeY - isoY)*cos(theta)-(tubeZ - isoZ)*sin(theta) )+isoY;
    rtubeZ = ((tubeY - isoY)*sin(theta)+(tubeZ - isoZ)*cos(theta) )+isoZ;

    % Detector rotation
    rdetY = ((detY - isoY).*cos(phi)-(detZ - isoZ).*sin(phi) )+isoY;
    rdetZ = ((detY - isoY).*sin(phi)+(detZ - isoZ).*cos(phi) )+isoZ;
    
    
    % For each slice
    for nz=1:nSlices
        
        % Temporary variable for each slice
        slice = zeros(nPixY,nPixX,'double');
                  
        % Map detector onto XY plane in the specifc slice   
        [detmX,detmY] = mapp2xyZ(tubeX,rtubeY,rtubeZ,detX,rdetY,rdetZ,objZ(nz));
        
        % Bilinear interpolation of the integrated proj(proj coords) on
        % image slice coords
        % (S.2. Interpolation - Liu et al (2017))
        sliceI = interp2(detmX,detmY,projI(:,:,p),objX,objY,'linear',0);
                
        for xPix=2:nPixX+1
    
            % x-ray angle in X coord
            gamma = atan((objX(1,xPix-1)+(deltaObjX/2)-tubeX)/(rtubeZ - objZ(nz))); 
            
            for yPix=2:nPixY+1                            
                
                % x-ray angle in Y coord
                alpha = atan((objY(yPix-1,1)+(deltaObjY/2)-rtubeY)/(rtubeZ - objZ(nz)));              
                
                % S.3. Differentiation - Eq. 24 - Liu et al (2017)
                slice(yPix-1,xPix-1) = slice(yPix-1,xPix-1)...
                                         +((sliceI(yPix,xPix)...
                                         -sliceI(yPix,xPix-1)...
                                         -sliceI(yPix-1,xPix)...
                                         +sliceI(yPix-1,xPix-1))...
                                         *(deltaDetX.*deltaDetY.*pixsZ./(cos(alpha)*cos(gamma)))...
                                         *(1/(deltaObjX*deltaObjY)));
                                                                              
            end % Loop end Y coord
        end % Loop end X coord
      
        % Accumulate for each proj angle
        data3d(:,:,nz) = data3d(:,:,nz) + slice;
        
    end % Loop end slices     
end % Loop end Projections

data3d = data3d ./ nProjs;

end

%% Integrate function
function [Ppj] = integrate2D (imageStack)

%         2-D image
%
%     | -> J
%     v -------------
%     I |           |
%       |           |
%       |           |
%       |           |
%       -------------

ni_pixel = size(imageStack,1)+1;
nj_pixel = size(imageStack,2)+1;
n_stack = size(imageStack,3);

p_v = zeros(ni_pixel,nj_pixel,n_stack,'double');    % Pad with zeros
p_v(2:end,2:end,:) = imageStack;            % Load slices

P_i = zeros(ni_pixel,1,n_stack,'double');
P_j = zeros(1,nj_pixel,n_stack,'double');

Ppj = zeros(ni_pixel,nj_pixel,n_stack,'double');

% Integrate on J direction
for pj=2:nj_pixel    
   P_i = P_i + p_v(:,pj,:);
   Ppj(:,pj,:) = P_i;
end

% Integrate on I direction
for pi=2:ni_pixel    
   P_j = P_j + Ppj(pi,:,:);
   Ppj(pi,:,:) = P_j;
end

end

%% Map on XY plane (on specific slice)
function [x,y] = mapp2xyZ(x1,y1,z1,x2,y2,z2,z)
    x = ((x2-x1).*(z-z2)-(x2.*z1)+(x2.*z2))./(-z1+z2);
    y = ((y2-y1).*(z-z2)-(y2.*z1)+(y2.*z2))./(-z1+z2);
end