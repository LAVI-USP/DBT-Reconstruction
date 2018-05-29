%% Author: Rodrigo de Barros Vimieiro
% Date: May, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                         readDicom(imgdir)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function read a set of Dicom images from a directory.
%  
%  
%     INPUT:
% 
%     - imgdir = Directory for Dicom files 
% 
%     OUTPUT:
% 
%     - dataDicom = Stack of Dicom images.
%     - infoDicom = Stack of Dicom headers.
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
%% Read Dicom images from files
function [dataDicom,infoDicom] = readDicom(imgdir)

img_list = dir([imgdir,'/*.dcm']);       % List dicom files
img_count = size(img_list,1);

for i=1:img_count
    imgAux = single(dicomread([imgdir, filesep ,  img_list(i,:).name]));  % Read Dicom
    infoAux = dicominfo([imgdir, filesep ,  img_list(i,:).name]);         % Read Dicom headers
    dataDicom(:,:,infoAux.InstanceNumber) = imgAux;
    infoDicom(:,infoAux.InstanceNumber) = infoAux;
end
end