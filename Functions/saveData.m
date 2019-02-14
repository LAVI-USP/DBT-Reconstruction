%% Author: Rodrigo de Barros Vimieiro
% Date: Jun, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%    saveData(dataRecon3d,parameter,title,folder,infoDicom)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function save reconstructed data.
%  
%  
%     INPUT:
% 
%     - dataRecon3d = Reconstructed data
%     - parameter = Parameter of all geometry
%     - title = Title of the data (Dicom, Shepp-Logan, VCT)
%     - infoDicom = Dicom header from projections
%     - folder = Folder to save the data
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
%% Save data
function saveData(dataRecon3d,varargin)

global gpuprocess

optionals = {[],[],[],[]}; % placeholder for inputs

numInputs = nargin - 1; % Minus number of required inputs
inputVar = 1;

while numInputs > 0
    if ~isempty(varargin{inputVar})
        optionals{inputVar} = varargin{inputVar};
    end
    inputVar = inputVar + 1;
    numInputs = numInputs - 1;
end

parameter = optionals{1};
title = optionals{2};
folder = optionals{3};
infoDicom = optionals{4};

if(isempty(folder)) 
    folder='output';
else
    mkdir([folder,filesep,title]);
end


% Some more informations to save
for k=1:size(dataRecon3d,2)
    info = dataRecon3d{2,k};
    if(gpuprocess)
        info.GPU = 'True';
    else
        info.GPU = 'False';
    end
    info.Title = title;
    dataRecon3d{2,k} = info;
end


filestring = [folder,filesep,title,filesep,'Recon',info.reconMeth,int2str(gpuprocess)];

if(~isempty(infoDicom))
    writeDicom(dataRecon3d,infoDicom,filestring,info.reconMeth,parameter);
end

save(filestring,'dataRecon3d','-v7.3')

end
