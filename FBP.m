%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 FBP(proj,filterType,cutoff,parameter)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function makes the Filtered Backprojection of 2D images, 
%     in order to reconstruct a certain number of slices.
%     
%     The geometry is for DBT with half cone-beam. All parameters are set 
%     in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - proj = 2D projection images 
%     - filterType = Filter type to be used
%     - cutoff = Cut off frequency up to Nyquist
%     - parameter = Parameter of all geometry
% 
%     OUTPUT:
% 
%     - reconData3d = Volume data reconstructed with FBP.
%     - time = [Filtering Time , Backproject Time]
% 
%     Reference: Jiang Hsieh's book (second edition)
%     Reference: Fessler Book -> (http://web.eecs.umich.edu/~fessler/book/c-tomo.pdf)
% 
%     -----------------------------------------------------------------------
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
%% Recon Code - Analytical reconstruction: FBP
function [reconData3d,time] = FBP(proj,filterType,cutoff,parameter)

global showinfo

if(showinfo)
    fprintf('----------------\nStarting FBP... \n\n')
end

tStart = tic;
if(strcmp(filterType,'FBP'))
    proj = filterProj(proj,parameter,cutoff);  % Filter projections
end
time = toc(tStart);
tStart = tic;

% Make the Backprojection
reconData3d = backprojection(proj,parameter);

% Store only interested slices and ROIs
reconData3d = reconData3d(parameter.iROI,parameter.jROI,parameter.sliceRange);

time = [time,toc(tStart)];
end
