%% Author: Rodrigo de Barros Vimieiro
% Date: July, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projectionDD(data3d,param)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates the volume projection based on the
%     Distance-Driven principle. It works by calculating the overlap in X 
%     and Y axis of the volume and the detector boundaries.  
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - data3d = 3D volume for projection 
%     - param = Parameter of all geometry
%
%     OUTPUT:
% 
%     - proj = projections for each angle.
% 
%     Reference: Three-Dimensional Digital Tomosynthesis - Yulia 
%     Levakhina (2014), Cap 3.6 and 3.7.
%
%     Original Paper: De Man, Bruno, and Samit Basu. "Distance-driven 
%     projection and backprojection in three dimensions." Physics in 
%     Medicine & Biology (2004).
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
%% 3-D Distance Driven Projection Code
function proj = projectionDD(data3d,param)

% Stack of projections
proj = zeros(param.nv, param.nu, param.nProj,'single');


% Map detector and object boudaries
param.us = (param.nu:-1:0)*param.du;
param.vs = (-(param.nv)/2:1:(param.nv)/2)*param.dv;
param.xs = (param.nx:-1:0)*param.dx;
param.ys = (-(param.ny)/2:1:(param.ny)/2)*param.dy;
param.zs = (0:1:param.nz)*param.dz + param.DAG + (param.dz/2);

% Detector boudaries coordinate sytem in (mm)
[detX,detY] = meshgrid(param.us,param.vs);
detZ = zeros(size(detX));

% Object boudaries coordinate sytem in (mm)
[objX,objY] = meshgrid(param.xs,param.ys);
objZ = param.zs;     % Z object coordinates

% X-ray tube initial position
tubeX = 0;
tubeY = 0;
tubeZ = param.DSD;

% Iso-center position
isoY = 0;
isoZ = param.DDR;

% Projection angles
tubeAngle = deg2rad(param.tubeDeg);
detAngle = deg2rad(param.detectorDeg);

% Number of detector and voxels for each direction
nDetX = param.nu;
nDetY = param.nv;
nPixX = param.nx;
nPixY = param.ny;

nSlices = param.nz;
nProjs = param.nProj;

% Voxel size on Z
pixsZ = param.dz; 

% For each projection
for p=1:nProjs
    
    % Temporary variable to acumulate all projection from the slices
    proj_tmp = zeros(nDetY,nDetX,'single');
    
    % Get specif tube angle for the projection
    theta = tubeAngle(p);
    
    % Get specif detector angle for the projection
    phi = detAngle(p);
    
    % Tubre rotation
    rtubeY = ((tubeY - isoY)*cos(theta)-(tubeZ - isoZ)*sin(theta) )+isoY;
    rtubeZ = ((tubeY - isoY)*sin(theta)+(tubeZ - isoZ)*cos(theta) )+isoZ;

    % Detector rotation
    rdetY = ((detY - isoY).*cos(phi)-(detZ - isoZ).*sin(phi) )+isoY;
    rdetZ = ((detY - isoY).*sin(phi)+(detZ - isoZ).*cos(phi) )+isoZ;
    
     
    % Map detector onto XY plane(Inside proj loop in case detector rotates)
    [detmX,detmY] = mapp2xy(tubeX,rtubeY,rtubeZ,detX,rdetY,rdetZ);
    
    % Z coord does not change in X direction, so we can take only one
    % collumn of Y mapped detecotr boundaries
    detmY = detmY(:,1);        
   
    % Pixel start index and increment
    detIstart = 1;
    detIinc = 1;    
    % Mapped detector length
    deltaDetmY = detmY(detIstart+detIinc)-detmY(detIstart);
        
    % For each slice
    for nz=1:nSlices
        
        % Flip X (Img coord is reverse to Global)
        slice = data3d(:,end:-1:1,nz);
          
        % Map slice onto XY plane
        [pixmX,pixmY] = mapp2xy(tubeX,rtubeY,rtubeZ,objX,objY,objZ(nz)); 
        
        %  - Z coord does not change in one slice, so we can take just one 
        %    line of mapped coordinates per slice
        %  - Flip X (Img coord is reverse to Global)
        pixmX = pixmX(1,end:-1:1);
        pixmY = pixmY(:,1);       
        
        % Pixel start index and increment
        pixIstart = 1;
        pixIinc = 1; 
        % Mapped pixel length
        deltaPixmX = pixmX(pixIstart+pixIinc)-pixmX(pixIstart);
        deltaPixmY = pixmY(pixIstart+pixIinc)-pixmY(pixIstart);
        
        % Start pixel and detector indices
        detIndY = detIstart;
        pixIndY = pixIstart;                 
        
        % Case 1
        % Find first detector overlap maped with pixel maped on Y
        if(detmY(detIndY)-pixmY(pixIstart)<-deltaDetmY)
            while(detmY(detIndY)-pixmY(pixIstart)<-deltaDetmY)            
                detIndY = detIndY + detIinc;            
            end
        else
        % Case 2    
        % Find first pixel overlap maped with detector maped on Y            
            if(detmY(detIstart)-pixmY(pixIndY)>deltaPixmY)
                while(detmY(detIstart)-pixmY(pixIndY)>deltaPixmY)            
                    pixIndY = pixIndY + pixIinc;            
                end
            end
        end
        
        % Get the left coordinate of the first overlap on Y axis
        % Try the following lines for a better understanding
        % % ---------------------------------------------------------------
        % % figure
        % % plot(detmY(:),zeros(1,size(detmY,1)),'r.','MarkerSize',6)
        % % hold on
        % % plot(pixmY(:),zeros(1,size(pixmY,1)),'b.','MarkerSize',6)
        % % hold off
        % % legend('Detector Boundaries','Pixel Boundaries')
        % % title('Overlap on Y axis')
        % % ---------------------------------------------------------------
        if( detmY(detIndY) < pixmY(pixIndY) )
            moving_left_boundaryY = pixmY(pixIndY);            
        else
            moving_left_boundaryY = detmY(detIndY);
        end       
        
        % Loop over Y intersections
        while((detIndY<=nDetY)&&(pixIndY<=nPixY))
             
            alpha = atan((detmY(detIndY)+(deltaDetmY/2)-rtubeY)/rtubeZ);
            
            % Case A, when you jump to the next detector boundarie but stay
            % in the same pixel
            if(detmY(detIndY+1)<=pixmY(pixIndY+1))
                overLapY = (detmY(detIndY+1)- moving_left_boundaryY)...
                    /deltaDetmY; % Normalized overlap Calculation
            else
            % Case B, when you jump to the next pixel boundarie but stay 
            % in the same detector
                overLapY = (pixmY(pixIndY+1)- moving_left_boundaryY)...
                    /deltaDetmY; % Normalized overlap Calculation
            end
                                 
            %% X overlap    
            detIndX = detIstart;
            pixIndX = pixIstart; 
            
            % Get row/coll of X, which correspond to that Y overlap det
            detmXrow = detmX(detIndY,end:-1:1);
            
            % Mapped detecor length on X
            deltaDetmX = detmXrow(detIstart+detIinc)-detmXrow(detIstart); 
            
            % Case 1
            % Find first detector overlap maped with pixel maped on X 
            if(detmXrow(detIndX)-pixmX(pixIstart)<-deltaDetmX)
                while(detmXrow(detIndX)-pixmX(pixIstart)<-deltaDetmX)            
                    detIndX = detIndX + detIinc;            
                end
            else
            % Case 2
            % Find first pixel overlap maped with detector maped on X            
                if(detmXrow(detIstart)-pixmX(pixIndX)>deltaPixmX)
                    while(detmXrow(detIstart)-pixmX(pixIndY)>deltaPixmX)            
                        pixIndX = pixIndX + pixIinc;            
                    end
                end
            end

            % Get the left coordinate of the first overlap on X axis
            % Try the following lines for a better understanding
            % % ---------------------------------------------------------------
            % % figure
            % % plot(detmXrow(:),zeros(1,size(detmXrow,1)),'r.','MarkerSize',6)
            % % hold on
            % % plot(pixmX(:),zeros(1,size(pixmX,1)),'b.','MarkerSize',6)
            % % hold off
            % % legend('Detector Boundaries','Pixel Boundaries')
            % % title('Overlap on X axis')
            % % ---------------------------------------------------------------

            if( detmXrow(detIndX) < pixmX(pixIndX) )
                moving_left_boundaryX = pixmX(pixIndX);            
            else
                moving_left_boundaryX = detmXrow(detIndX);
            end
        
                
            % Loop over X intersections
            while((detIndX<=nDetX)&&(pixIndX<=nPixX))
                 
                gamma = atan((detmXrow(detIndX)+(deltaDetmX/2)-tubeX)/rtubeZ);                
                
                % Case A, when you jump to the next detector boundarie but stay
                % in the same pixel
                if(detmXrow(detIndX+1)<=pixmX(pixIndX+1))
                    
                    overLapX = (detmXrow(detIndX+1)- moving_left_boundaryX)...
                    /deltaDetmX; % Normalized overlap Calculation
                                            
                    proj_tmp((detIndX-1)*nDetY+detIndY) = proj_tmp((detIndX-1)*nDetY+detIndY) ...
                    + overLapX * overLapY * slice((pixIndX-1)*nPixY+pixIndY )...
                    * pixsZ/(cos(alpha)*cos(gamma));               
                
                    detIndX = detIndX + detIinc;
                    moving_left_boundaryX = detmXrow(detIndX);
                    
                else
                % Case B, when you jump to the next pixel boundarie but stay 
                % in the same detector
                
                    overLapX = (pixmX(pixIndX+1)- moving_left_boundaryX)...
                    /deltaDetmX; % Normalized overlap Calculation
                                            
                    proj_tmp((detIndX-1)*nDetY+detIndY) = proj_tmp((detIndX-1)*nDetY+detIndY) ...
                    + overLapX * overLapY * slice((pixIndX-1)*nPixY+pixIndY )...
                    * pixsZ/(cos(alpha)*cos(gamma));
                
                    pixIndX = pixIndX + pixIinc;
                    moving_left_boundaryX = pixmX(pixIndX);
                    
                end
            end
           
            %% Back to Y overlap
            
            % Case A, when you jump to the next detector boundarie but stay
            % in the same pixel
            if(detmY(detIndY+1)<=pixmY(pixIndY+1))
                detIndY = detIndY + detIinc;
                moving_left_boundaryY = detmY(detIndY);
            else
            % Case B, when you jump to the next pixel boundarie but stay 
            % in the same detector
                pixIndY = pixIndY + pixIinc;
                moving_left_boundaryY = pixmY(pixIndY);
            end
            
        end % Y Overlap loop
                     
    end % Loop end slices
    
    proj(:,:,p) = fliplr(proj_tmp);    
    
end % Loop end Projections

end
%% Map on XY plane
function [x,y] = mapp2xy(x1,y1,z1,x2,y2,z2)
    x = x2 + z2 .* (x2-x1)./(z1-z2);
    y = y2 + z2 .* (y2-y1)./(z1-z2);
end