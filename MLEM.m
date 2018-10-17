%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                         MLEM(proj,nIter,parameter)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function reconstruct iteratively the volume through 
%     Maximum-Likelihood Expectation-Maximization (MLEM) method.
%     
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
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
%     Reference: Medical Image Reconstruction A Conceptual Tutorial - 
%     Gengsheng Lawrence Zeng (2010)
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
%% Recon Code -  Iterative reconstruction: MLEM
function allreconData3d = MLEM(proj,nIter,parameter)

global showinfo

info.startDateAndTime = char(datetime('now','Format','MM-dd-yyyy''  ''HH:mm:ss'));
info.reconMeth = 'MLEM';

highestValue = (2^parameter.bitDepth) - 1;

allIterTime = zeros(nIter(end),1);   % Time data for each iteration
allreconData3d = cell(1,size(nIter,2)); % Recon data for each iteration

% Initial estimated volume data 
reconData3d = zeros(parameter.ny, parameter.nx, parameter.nz,'single');
reconData3d(:) = 1;

% Volume normalization
vol_norm = backprojection(ones(parameter.nv, parameter.nu, parameter.nProj, 'single'), parameter);

if(showinfo)
    fprintf('----------------\nStarting MLEM Iterations... \n\n')
end

% Start Iterations
for iter = 1:nIter(end)
    
    tStart = tic;
    
    % Error ratio between raw data and projection of estimated data  
    proj_ratio = proj./projection(reconData3d,parameter);
    proj_ratio(isnan(proj_ratio)) = 0;
    proj_ratio(isinf(proj_ratio)) = 0;

    upt_term = backprojection(proj_ratio, parameter);
    upt_term = upt_term ./vol_norm; % Volume normalization
    upt_term(isnan(upt_term)) = 0;  
    upt_term(isinf(upt_term)) = 0;

    reconData3d = reconData3d.*upt_term; % Updates the previous estimation 
        
    % Truncate to highest value
    reconData3d(reconData3d>highestValue) = highestValue;
    
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
    
end% Loop end iterations
end
