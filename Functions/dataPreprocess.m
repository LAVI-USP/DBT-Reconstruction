%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                    dataPreprocess(data,parameter)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function pre process the projection images.
%  
%  
%     INPUT:
% 
%     - data = Projections.
%     - parameter = Parameters.
% 
%     OUTPUT:
% 
%     - data_mod = Modified image.
%     - parameter_mod = Modified parameters.
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
%
%}
% =========================================================================
function [data_mod,parameter_mod] = dataPreprocess(data,parameter)

parameter_mod = parameter;

%% Segmentation
Gap = 20;
meanSlideW = 251;

maxValue = single(max(max(data(:,:,5))));

vertProj = sum(abs(maxValue-data(:,:,5)));	% Horizontal Profile
[~,Ind] = max(diff(movmean(vertProj,meanSlideW)));		% Smooth the signal and takes its max positive derivative
Ind = Ind-meanSlideW-Gap;	% Subtract the ind found to a gap

data_mod = data(:,Ind:end,:);

% Modifies parameters based on segmentation
parameter_mod.nu = size(data_mod,2);    % Number of pixels (columns)
parameter_mod.nv = size(data_mod,1);    % number of pixels (rows) 
parameter_mod.su = parameter_mod.nu.*parameter_mod.du;	
parameter_mod.sv = parameter_mod.nv.*parameter_mod.dv; 
parameter_mod.us = (parameter_mod.nu-1:-1:0)*parameter_mod.du;
parameter_mod.vs = (-(parameter_mod.nv-1)/2:1:(parameter_mod.nv-1)/2)*parameter_mod.dv;

%% Transform intensity image in attenuation coefficients
data_mod = -log(data_mod./single(2^parameter_mod.bitDepth-1));

if(numel(find(data_mod==inf)) ~= 0)
    error('Numeric error. Please, double-check if you input data has zeros on it')
end

end
