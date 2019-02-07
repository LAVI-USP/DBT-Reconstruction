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
%     - reconData3d{1} = Volume data reconstructed with FBP.
%     - reconData3d{2} = Struct informations like - Filter and Bp Time
% 
%     Reference: Jiang Hsieh's book (second edition)
%     Reference: Fessler, J. A.: Fundamentals of CT Reconstruction in 2D 
%     and 3D. In: Brahme, A. (eds.) Comprehensive Biomedical Physics, 
%     1st ed., vol. 2, pp 263-295. Elsevier, Netherlands (2014).
%     doi:10.1016/B978-0-444-53632-7.00212-4.
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
function reconData3d = FBP(proj,filterType,cutoff,parameter)

global showinfo

info.startDateAndTime = char(datetime('now','Format','MM-dd-yyyy''  ''HH:mm:ss'));
info.reconMeth = filterType;

if(showinfo)
    fprintf('----------------\nStarting FBP... \n\n')
end

if(strcmp(filterType,'FBP'))
    tStart = tic;
    proj = filterProj(proj,parameter,cutoff);  % Filter projections
    info.FilterTime = num2str(toc(tStart));
    info.CutOffFreq = num2str(cutoff);
end

tStart = tic;
% Make the Backprojection
reconData3d = backprojection(proj,parameter,[]);
info.BackprojectTime = num2str(toc(tStart));

% Store only interested slices and ROIs
reconData3d = {reconData3d(parameter.iROI,parameter.jROI,parameter.sliceRange)};
reconData3d{2,1} = info;

end
