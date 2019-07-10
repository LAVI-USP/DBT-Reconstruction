%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projectionDDb(data3d,param,projNumber)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates the volume projection based on the
%     Branchless Distance-Driven principle.  
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - data3d = 3D volume for projection 
%     - param = Parameter of all geometry
%     - projNumber = Vector with projections numbers to be processed
%
%     OUTPUT:
% 
%     - proj = projections for each angle.
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
%% 3-D Projection Branchless Distance Driven Code
function proj = projectionDDb(data3d,param,projNumber)

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

% Stack of projections
proj = zeros(param.nv, param.nu, nProjs,'double');

% Integration of 2D slices over the whole volume 
% (S.1. Integration. - Liu et al (2017))
data3dI = integrate2D(data3d);

% For each projection
for p=1:nProjs
    
    % Temporary variable to acumulate all projection from the slices
    proj_tmp = zeros(nDetY,nDetX,'double');
    
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

        % Map detector onto XY plane in the specifc slice   
        [detmX,detmY] = mapp2xyZ(tubeX,rtubeY,rtubeZ,detX,rdetY,rdetZ,objZ(nz));
        
        % Bilinear interpolation of the integrated image(slice coords) on
        % mapped detector coords
        % (S.2. Interpolation - Liu et al (2017))
        proj_data3dI = interp2(objX,objY,data3dI(:,:,nz),detmX,detmY,'linear',0);
        
        %figure
        %plot(detmX,detmY,'k*','MarkerSize',6)
        %hold on
        %plot(objX,objY,'k.','MarkerSize',6)
        %hold off
        %legend('Detector Boundaries','Pixel Boundaries')
        
        % These two lines find the max value of the integrated image, in 
        % order to replicate it in the remaining rows
        [~,row_max] = max(proj_data3dI(:,end));   
        proj_data3dI(row_max+1:end,:) = repmat(proj_data3dI(row_max,:),param.nv+1-row_max,1);
        
        for xDet=2:nDetX+1
              
            % x-ray angle in X coord
            gamma = atan((detX(1,xDet-1)+(deltaDetX/2)-tubeX)/rtubeZ); 
            
            for yDet=2:nDetY+1
                
                % Detector X coord overlap calculation
                overLapX = detmX(yDet,xDet-1)-detmX(yDet,xDet);
                
                % Detector Y coord overlap calculation
                overLapY = detmY(yDet,1)-detmY(yDet-1,1);
                
                 % x-ray angle in Y coord
                alpha = atan((rdetY(yDet-1,1)+(deltaDetY/2)-rtubeY)/rtubeZ);              
                             
                % S.3. Differentiation - Eq. 24 - Liu et al (2017)
                proj_tmp(yDet-1,xDet-1) = proj_tmp(yDet-1,xDet-1)...
                                         +((proj_data3dI(yDet,xDet)...
                                         -proj_data3dI(yDet,xDet-1)...
                                         -proj_data3dI(yDet-1,xDet)...
                                         +proj_data3dI(yDet-1,xDet-1))...
                                         *(deltaObjX*deltaObjX.*pixsZ./(cos(alpha)*cos(gamma)*overLapX*overLapY)));                                                                              
            end % Loop end Y coord
        end % Loop end X coord
    
    end % Loop end slices
   
    proj(:,:,p) = proj_tmp;

end % Loop end Projections
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