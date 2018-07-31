%% Author: Rodrigo de Barros Vimieiro
% Date: Jun, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%    writeDicom(dataRecon3d,infoDicomRecon,filestring,reconmeth,parameter)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function write a set of Dicom images with respective headers.
%  
%  
%     INPUT:
% 
%     - dataRecon3d = Reconstructed data
%     - infoDicomProjs = Dicom header from projections
%     - filestring = Path to save
%     - reconmeth = Reconstructed method
%     - parameter = Parameter of all geometry
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
%% Write Dicom images
function writeDicom(dataRecon3d,infoDicomRecon,filestring,reconmeth,parameter)

% Get reconstructed dicom headers
% uiwait(msgbox('Select reconstructed image Dicom file to get headers.','Dicom','Warn'));
% [filename,path] = uigetfile('*.dcm');
% infoDicomRecon = dicominfo([path,filename]);

% Weather is Hologic or GE
dc=0;
if(strcmp(parameter.type,'vct')||strcmp(parameter.type,'hologic'))
    dc = -1;
end

% Adjust data to save, based on the recon method
if(strcmp(reconmeth,'BP')||strcmp(reconmeth,'MLEM')||strcmp(reconmeth,'SART'))
    for k=1:size(dataRecon3d,2)
        dataRecon3d{k} = (2^14-1).*mat2gray(dataRecon3d{k});
        dataRecon3d{k} = uint16(dataRecon3d{k});
    end
else
    if(strcmp(reconmeth,'FBP'))
        for k=1:size(dataRecon3d,2)
            dataRecon3d{k} = dataRecon3d{k} + abs(min(dataRecon3d{k}(:)));
            dataRecon3d{k} = (2^14-1).*mat2gray(dataRecon3d{k});
            dataRecon3d{k} = uint16(dataRecon3d{k});
        end
        
    end
end
  
% Write Dicom files in respective folder
for k=1:size(dataRecon3d,2) 
    filedir = [filestring,filesep,int2str(k)];
    mkdir(filedir)
    for i=1:size(dataRecon3d{k},3)
        fileimg = [filedir,filesep,num2str(i+dc),'.dcm'];
        dicomwrite(dataRecon3d{k}(:,:,i), fileimg, infoDicomRecon(:,1), 'CreateMode', 'copy');
    end
end

end