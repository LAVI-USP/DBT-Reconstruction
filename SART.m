%% Author: Rodrigo de Barros Vimieiro
% Date: May, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                       SART(proj,nIter,parameter)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function reconstruct iteratively the volume through 
%     Simultaneous Algebraic Reconstruction Technique (SART) method.
%     
%     The geometry is for DBT with half cone-beam. All parameters are set in 
%     "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - proj  = 2D projection images 
%     - nIter = Specific iterations to save volume data 
%     - parameter = Parameter of all geometry
% 
%     OUTPUT:
% 
%     - allreconData3d{1,...} = Cell with Volumes of each specific iteration.
%     - allreconData3d{2,...} = Cell with informations of each specific iteration.
% 
%     Reference: Three-Dimensional Digital Tomosynthesis - Yulia Levakhina (2014)
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
%% Recon Code -  Iterative reconstruction: SART
function allreconData3d = SART(proj,nIter,parameter)

global showinfo animation

numProjs = parameter.nProj;

info.startDateAndTime = char(datetime('now','Format','MM-dd-yyyy''  ''HH:mm:ss'));
info.reconMeth = 'SART';

highestValue = (2^parameter.bitDepth) - 1;

allIterTime = zeros(nIter(end),1);   % Time data for each iteration
allreconData3d = cell(1,size(nIter,2)); % Recon data for each iteration

% Initial estimated data 
reconData3d = zeros(parameter.ny, parameter.nx, parameter.nz,'single');

% Pre calculation of Projection normalization
tempvar = animation; animation = 0; 
proj_norm = projection(ones(parameter.ny, parameter.nx, parameter.nz, 'single'),parameter, []);
animation = tempvar; clear tempvar;

% Pre calculation of Backprojection normalization
vol_norm = backprojection(ones(parameter.nv, parameter.nu, parameter.nProj, 'single'), parameter, []);

if(showinfo)
    fprintf('----------------\nStarting SART Iterations... \n\n')
end

% Start Iterations
for iter = 1:nIter(end)
    tStart = tic;
    % For each projection
    for p=1:numProjs
        
        % Error between raw data and projection of estimated data 
        proj_diff = proj(:,:,p) - projection(reconData3d,parameter,p);  

        proj_diff = proj_diff ./ proj_norm(:,:,p); % Projection normalization
        proj_diff(isnan(proj_diff)) = 0;
        proj_diff(isinf(proj_diff)) = 0;

        upt_term = backprojection(proj_diff,parameter,p);
        upt_term = upt_term ./ vol_norm; % Volume normalization
        upt_term(isnan(upt_term)) = 0;
        upt_term(isinf(upt_term)) = 0;

        % Updates the previous estimation 
        reconData3d = reconData3d + upt_term; 
    end
        
    % Truncate to highest and minimum value
    reconData3d(reconData3d>highestValue) = highestValue;
    reconData3d(reconData3d<0) = 0;

    allIterTime(iter,1) = toc(tStart);
    
    % Save data 
    indIter = find(nIter == iter);
    if(indIter~=0)
        allreconData3d{1,indIter} = reconData3d(parameter.iROI,parameter.jROI,parameter.sliceRange);
        info.IterationNumber = num2str(iter);
        info.IterationTime = num2str(sum(allIterTime(:)));
        allreconData3d{2,indIter} = info;      
    end
    
    if(showinfo) 
        fprintf('Iteration %d Time: %.2f\n\n',iter,allIterTime(iter,1)); 
    end
    
end
end
