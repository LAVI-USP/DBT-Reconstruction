%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
---------------------------------------------------------------------------
                        readDicom(imgdir)
---------------------------------------------------------------------------
    DESCRIPTION:
    This function read a set of Dicom images from a directory.
 
 
    INPUT:

    - imgdir = Directory for Dicom files 

    OUTPUT:

    - dataDicom = Stack of Dicom images.
    - infoDicom = Stack of Dicom headers.

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
%% Read Dicom images from files
function [dataDicom,infoDicom] = readDicom(imgdir)

if isunix
    img_list = ls([imgdir,'/*.dcm']);       % List dicom files
    img_list = strsplit(img_list,'\n')';    % Split files
    img_list(end,:) = [];                   % Get ride off last one
else
    img_list = ls([imgdir,'\*.dcm']);
end

img_count = size(img_list,1);

for i=1:img_count
    filename = char(strtrim(img_list(i,:)));    % Remove white space from string
    if isunix
        imgAux = single(dicomread(filename));   % Read Dicom       
        infoAux = dicominfo(filename);          % Read Dicom headers
    else
        imgAux = single(dicomread([imgdir, '\' ,  filename]));  % Read Dicom
        infoAux = dicominfo([imgdir, '\' ,  filename]);         % Read Dicom headers
    end
    dataDicom(:,:,infoAux.InstanceNumber) = imgAux(:,:,1);
    infoDicom(:,infoAux.InstanceNumber) = infoAux;
end

end