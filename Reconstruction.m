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

%% Global parameters
global showinfo saveinfo animation gpuprocess

showinfo = uint8(1);        % Show projection animation
saveinfo = uint8(1);        % Save reconstructed volume
animation = uint8(1);       % Graphical animation
gpuprocess = uint8(0);      % Pprocessing on GPU

%% GUI - Data decision

answer = questdlg('Load Dicom file or create a virtual Phantom?', ...
	'Data decision', ...
	'Dicom','Phantom','Dicom');
% Handle response
switch answer
    case 'Dicom'
        fprintf('Waiting for Dicom files \n')
        data = 1;
    case 'Phantom'
        fprintf('Creating a virtual phantom \n')
        data = 2;
    otherwise
        fprintf('Cancelled by user \n');
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
           
    if(~exist('res/Dicom','dir'))
        mkdir('res/Dicom')
    end
    
    % Load Projection Data    
    uiwait(msgbox('Select the path of projection Dicom files.','Dicom','Warn'));
    path_User = userpath;
    path_ProjData = uigetdir(path_User(1:end-1)); 
    
    if(path_ProjData == 0)
        fprintf('Cancelled by user \n');
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
    
    if(~exist('res/Phantom','dir'))
        mkdir('res/Phantom')
    end
    addpath('res/Phantom');

    % Create Shepp-Logan phantom
    data3d = single(phantom3d('Modified Shepp-Logan', parameter.nz));   
    data3d(data3d<0) = eps;

    % Make the Projections
    if(animation || saveinfo)
        % Create a figure of screen size
        figureScreenSize()
        dataProj = projection(data3d,parameter);
        if(saveinfo)
            save('res/Phantom/proj.mat','dataProj')
        else
            load proj.mat
        end
    else
        load proj.mat
    end
    
end

if(gpuprocess)
    fprintf(2,'Processing on GPU \n');
end
fprintf('Starting reconstruction \n');

%% Set specific recon parameters

nIter = [4,5];          % Iteration to be saved (In case of MLEM or SART)
filterType = 'FBP';     % Filter type: 'BP', 'FBP'

%% Reconstruction methods

%                       ** Uncomment to use **
[dataRecon3d,time] = FBP(dataProj,filterType,parameter);
if(saveinfo)
    filestring = ['res/',answer,'/Recon',filterType,int2str(gpuprocess)];
    save(filestring,'dataRecon3d')
    xlswrite(filestring,time)
end

%                       ** Uncomment to use **
% [dataRecon3d,time] = SART(dataProj,nIter,parameter);
% if(saveinfo)
%     filestring = ['res/',answer,'/ReconSART',int2str(gpuprocess)];
%     save(filestring,'dataRecon3d','-v7.3')
%     xlswrite(filestring,time)   
% end

%                       ** Uncomment to use **
% [dataRecon3d,time] = MLEM(dataProj,nIter,parameter);
% if(saveinfo)
%     filestring = ['res/',answer,'/ReconMLEM',int2str(gpuprocess)];
%     save(filestring,'dataRecon3d','-v7.3')
%     xlswrite(filestring,time)  
% end

fprintf('Finished \n');
fprintf('Please, check "res" folder for results \n');

%% Show info

if(animation && data == 2)
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