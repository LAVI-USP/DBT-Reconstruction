%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
---------------------------------------------------------------------------
                MLEM(proj,nIter,parameter,showinfo)
---------------------------------------------------------------------------
    DESCRIPTION:
    This function makes a iterative reconstruction method. MLEM stands for
    Maximum likelihood Expectation Maximization.
    
    The geometry is for DBT with half cone-beam. All parameters are set in 
    "ParameterSettings" code. 
 
    INPUT:

    - proj  = 2D projection images 
    - nIter = Specific iterations to save volume data 
    - param = Parameter of all geometry
    - showinfo = Show information on Command Window

    OUTPUT:

    - allreconData3d = Cell with Volumes of each specific iteration.
    - allIterTime = time of all iterations

    Reference: Medical Image Reconstruction A Conceptual Tutorial - 
    Gengsheng Lawrence Zeng (2010)

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
%% Recon Code -  Iterative reconstruction: MLEM
function [allreconData3d,allIterTime] = MLEM(proj,nIter,parameter,showinfo)

highestValue = (2^parameter.bitDepth) - 1;

allIterTime = zeros(nIter(end),1);   % Time data for each iteration
allreconData3d = cell(1,size(nIter,2)); % Recon data for each iteration

% Inicial estimation (Set to maximum)
reconData3d = zeros(parameter.ny, parameter.nx, parameter.nz,'single');
reconData3d(parameter.iROI,parameter.jROI,parameter.sliceRange) = highestValue;

% Volume normalization
NorVol = backprojection(ones(parameter.nv, parameter.nu, parameter.nProj, 'single'), parameter);

if(showinfo)
    fprintf('Starting Iterations... \n\n\n')
end

% Start Iterations
for iter = 1:nIter(end)
    
    tStart = tic;
    
    proj_ratio = proj./projection(reconData3d,parameter,0);
    proj_ratio(isnan(proj_ratio)) = 0;
    proj_ratio(isinf(proj_ratio)) = 0;

    img_ratio = backprojection(proj_ratio, parameter)./NorVol;
    img_ratio(isnan(img_ratio)) = 0;
    img_ratio(isinf(img_ratio)) = 0;

    reconData3d = reconData3d.*img_ratio;
    
    allIterTime(iter,1) = toc(tStart);
    
    % Truncate to highest value
    reconData3d(reconData3d>highestValue) = highestValue;
    
    % Save data 
    indIter = find(nIter == iter);
    if(indIter~=0)
        allreconData3d{1,indIter} = reconData3d;
    end
    
    if(showinfo) 
        fprintf('Itaration %d Time: %.2f\n\n',iter,toc(tStart));  
    end
    
end% Loop end iterations
end