%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{

    DESCRIPTION:

    The goal of this software is to present an open source reconstruction 
    toolbox, which features the three basic types of DBT reconstruction, 
    with the possibility of different acquisition geometries. This toolbox 
    is intended for academic usage and it aims at expanding the academic 
    research in this field.  Since it is an open source toolbox, 
    researchers are welcome to contribute to the current version of the 
    software.

    Department of Electrical and Computer Engineering, 
    São Carlos School of Engineering, University of São Paulo, 
    São Carlos, Brazil

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
%%                      Reconstruction Code                              %%
close all;clear;clc

% Global parameters
showinfo = uint8(1);        % Show projection animation
saveinfo = uint8(1);        % Save reconstructed volume

%% GUI - Data decision

answer = questdlg('Load Dicom file or create a virtual Phantom?', ...
	'Data decision', ...
	'Dicom','Phantom','Dicom');
% Handle response
switch answer
    case 'Dicom'
        disp(' Waiting for Dicom files.')
        data = 1;
    case 'Phantom'
        disp(' Creating a virtual phantom.')
        data = 2;
    otherwise
        display('Cancelled by user.');
    	return;    
end


%% Load components

addpath('Functions');
addpath('Parameters');

if(~exist('res','dir'))
    mkdir('res')
    saveinfo = 1;
end
addpath('res');



if(data == 1)   %           ** Dicom data **
    
    ParameterSettings_GE
           
    % Load Projection Data    
    uiwait(msgbox('Select the path of projection Dicom files.','Dicom','Warn'));
    path_User = userpath;
    path_ProjData = uigetdir(path_User(1:end-1)); 
    
    if(path_ProjData == 0)
        display(' Cancelled by user.');
    	return;
    else
        userpath(path_ProjData)       
        
        [dataProj,infoDicom] = readDicom(path_ProjData);
        parameter.bitDepth = infoDicom(:,1).BitDepth;
        % Pre process projections
        [dataProj,parameter] = dataPreprocess(dataProj,parameter);
    end
    
else    %                   ** Phantom data **
    
    ParameterSettings_Phantom;

    % Create Shepp-Logan phantom
    data3d = single(phantom3d('Modified Shepp-Logan', parameter.nz));   
    data3d(data3d<0) = eps;

    % Make the Projections
    if(showinfo || saveinfo)
        % Create a figure of screen size
        figureScreenSize()
        dataProj = projection(data3d,parameter,showinfo);
        if(saveinfo)
            save('res/proj.mat','dataProj')
        else
            load proj.mat
        end
    else
        load proj.mat
    end
end
display(' Starting reconstruction...');
 %% Set specific recon parameters

nIter = [5,6];       % Iteration to be saved (In case of MLEM)
filterType = uint8(1);      % Filter type (In case of FBP) -> 0 - No filter  1 - Filter

%% Reconstruction methods

%                       ** Uncomment to use **
[dataRecon3d,time] = FBP(dataProj,filterType,parameter);
if(saveinfo)
    save('res/ReconFBP.mat','dataRecon3d')
end

%                       ** Uncomment to use **
% [dataRecon3d,time] = SART(dataProj,nIter,parameter,showinfo);
% if(saveinfo)
%     save('res/ReconSART.mat','dataRecon3d','-v7.3')
%     xlswrite('res/SART_IterValues',time)   
% end

%                       ** Uncomment to use **
% [dataRecon3d,time] = MLEM(dataProj,nIter,parameter,showinfo);
% if(saveinfo)
%     save('res/ReconMLEM.mat','dataRecon3d','-v7.3')
%     xlswrite('res/MLEM_IterValues',time)   
% end

display(' Finished');
display(' Please, check "res" folder for results');

%% Show info

if(showinfo && data == 2)
    % Create a figure of screen size
    figureScreenSize()
    % 3D Backprojection Visualization
    for k=1:parameter.nz
        subplot(1,2,1)
        imshow(data3d(:,:,k))
        title(['Slice Orig ',num2str(k)]);axis on;
        subplot(1,2,2)
        imshow(dataRecon3d(:,:,k),[])    % [.09 .24] limit for MLEM-15
        title(['Slice Recon ',num2str(k)]);axis on;
        pause(0.08)
    end
end