%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projection(data3d,param,projNumber)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates for each detector pixel, which voxel is
%     associated with that pixel in the specific projection. That is done 
%     for all angles specified.
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
%     - proj = projection for each angle.
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
%% Projection Code
function proj = projection(data3d,param,projNumber)

global animation

% Get parameters from struct
DSR = param.DSR;
DDR = param.DDR;
tubeAngle = deg2rad(param.tubeDeg);
numVPixels = param.nv;
numUPixels = param.nu;
numYVoxels = param.ny;
numXVoxels = param.nx;
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

% Stack of projections
proj = zeros(param.nv, param.nu, numProjs,'single');


% Detector Coordinate sytem in (mm)
[uCoord,vCoord] = meshgrid(param.us,param.vs);

% Object Coordinate sytem in (mm) (just for 3D visualization)
[xCoord,yCoord] = meshgrid(param.xs,param.ys);
[   ~  ,zCoord] = meshgrid(param.ys,param.zs);

% For each projection
for p=1:numProjs
    
    % Get specific projection number
    projN = projNumber(p);
    
    % Get specific tube angle for the projection
    teta = tubeAngle(projN);
    
    % Temporary projection variable to acumulate all projection from the slices
    proj_tmp = zeros(numVPixels,numUPixels,'single');
    
    % For each slice
    for nz=1:numSlices
          
        % Calculates a relation of detector and voxel coordinates
        pyCoord = ((vCoord.*((DSR.*cos(teta))+DDR-zCoords(nz)))-(zCoords(nz).*DSR.*sin(teta)))./...
                                            (DSR.*cos(teta)+DDR);
        
        pxCoord = (uCoord.*((DSR.*cos(teta))+DDR-zCoords(nz)))./...
                             ((DSR.*cos(teta))+DDR);          
        
        % 3D Visualization
        if(nz==numSlices && animation == 1)
            figure(1)
            % Calculus of projection on the detector for visualization
            pvCoord = yCoord+((zCoords(nz).*((DSR.*sin(teta))+ yCoord))./...
                                 ((DSR.*cos(teta))-zCoords(nz)));

            puCoord = (xCoord.*DSR.*cos(teta))./ ...
                   ((DSR.*cos(teta))-zCoords(nz));
               
            % Draw 3D animation   
            draw3d(xCoord,yCoord,zCoord,puCoord,pvCoord,param,projN,teta);
        end              
             
        % Coordinate in Pixel of slice plane axis origin
        x0 = numXVoxels;
        y0 = numYVoxels/2;
        
        % Represent slice plane axis (X,Y) in the image plane axis (i,j) (covert mm to Pixels)
        pxCoord = -(pxCoord./param.dx) + x0;
        pyCoord =  (pyCoord./param.dy) + y0;  
      
        % Interpolation of the pixel coordinates of the "Slice" at the calculated pixel coordinates for the detector
        proj_tmp = proj_tmp + interp2(data3d(:,:,nz),pxCoord,pyCoord,'linear',0);
              
    end % Loop end slices
    
    proj(:,:,p) = proj_tmp;
   
    if(animation == 1)
        % 2D Visualization
        figure(1)
        subplot(2,2,2)
        imshow(imrotate(data3d(:,:,round(numSlices/2)),90),[])
        title('Slice');axis on;
        subplot(2,2,4)
        imshow(imrotate(proj(:,:,p),90),[]); title(['Proj ',num2str(projN)]);axis on;
    end 
    
end % Loop end Projections
end