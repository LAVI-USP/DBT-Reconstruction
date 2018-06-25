%% Author: Rodrigo de Barros Vimieiro
% Date: Jun, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%    saveData(dataRecon3d,parameter,answer,infoDicom)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function save reconstructed data.
%  
%  
%     INPUT:
% 
%     - dataRecon3d = Reconstructed data
%     - parameter = Parameter of all geometry
%     - answer = Type of data (Dicom, Shepp-Logan, VCT)
%     - infoDicom = Dicom header from projections
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
function saveData(dataRecon3d,parameter,answer,infoDicom)

global gpuprocess

% Some more informations to save
for k=1:size(dataRecon3d,2)
    info = dataRecon3d{2,k};
    if(gpuprocess)
        info.GPU = 'True';
    else
        info.GPU = 'False';
    end
    info.Title = answer;
    dataRecon3d{2,k} = info;
end


filestring = ['res',filesep,answer,filesep,'Recon',info.reconMeth,int2str(gpuprocess)];

if(~strcmp(answer,'Shepp-Logan'))
    writeDicom(dataRecon3d,infoDicom,filestring,info.reconMeth,parameter);
end
save(filestring,'dataRecon3d','-v7.3')

end