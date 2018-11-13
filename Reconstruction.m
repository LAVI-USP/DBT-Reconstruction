%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% 
%     DESCRIPTION:
% 
%     The goal of this software is to present an open source reconstruction 
%     toolbox, which features the three basic types of DBT reconstruction, 
%     with the possibility of different acquisition geometries. This 
%     toolbox is intended for academic usage and it aims at expanding the 
%     academic research in this field.  Since it is an open source toolbox, 
%     researchers are welcome to contribute to the current version of the 
%     software.
% 
%     Department of Electrical and Computer Engineering, 
%     São Carlos School of Engineering, University of São Paulo, 
%     São Carlos, Brazil
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
%%                      Reconstruction Code                              %%
close all;clear;clc

%% Global parameters
global showinfo saveinfo animation gpuprocess noise

showinfo = uint8(1);        % Show projection animation
saveinfo = uint8(1);        % Save reconstructed volume
animation = uint8(1);       % Graphical animation
gpuprocess = uint8(0);      % Processing on GPU
noise = uint8(0);           % Add noise to phantom projections

%% GUI - Data decision

answer = questdlg('Load Dicom file or create a Shepp-Logan Phantom?', ...
	'Data decision', ...
	'Dicom','Shepp-Logan','Dicom');
% Handle response
switch answer
    case 'Dicom'
        answer = questdlg('Clinical or VCT?', ...
            'Data decision', ...
            'Clinical','VCT','Clinical');
        switch answer
            case 'Clinical'
                fprintf('Waiting for Clinical images \n')
                data = 1;
            case 'VCT'
                fprintf('Waiting for VCT images \n')
                data = 2;
            otherwise
                fprintf('Cancelled by user \n');
                return
        end
    case 'Shepp-Logan'
        fprintf('Creating a Shepp-Logan phantom \n')
        data = 3;
    otherwise
        fprintf('Cancelled by user \n');
    	return;    
end

%% Load components

addpath(genpath('Functions'));
addpath(genpath('Parameters'));

if(~exist('output','dir'))
    mkdir('output')
    saveinfo = 1;
end
addpath('res');

if(data == 1 || data == 2)   %   ** Dicom data **
    
    if(data==1)
        ParameterSettings_GE
        if(~exist(['res',filesep,'Clinical'],'dir'))
            mkdir(['res',filesep,'Clinical'])
        end
    else
        ParameterSettings_VCT
        if(~exist(['res',filesep,'VCT'],'dir'))
            mkdir(['res',filesep,'VCT'])
        end
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
        
        [dataProj,infoDicom] = readDicom(path_ProjData,parameter);
        parameter.bitDepth = infoDicom(:,1).BitDepth;
        % Pre process projections
        [dataProj,parameter] = dataPreprocess(dataProj,parameter);
    end
    
else    %                   ** Shepp-Logan data **
    
    ParameterSettings_Phantom;
    
    if(~exist(['res',filesep,'Shepp-Logan'],'dir'))
        mkdir(['res',filesep,'Shepp-Logan'])
    end
    addpath(['res',filesep,'Shepp-Logan']);

    % Create Shepp-Logan phantom
    data3d = single(phantom3d('Modified Shepp-Logan', parameter.nz));   
    data3d(data3d<0) = eps;

    % Make the Projections
    if(animation || saveinfo)
        % Create a figure of screen size
        figureScreenSize()
        dataProj = projection(data3d,parameter);
        if(noise)
            dataProj = insertNoise(dataProj,parameter); % Insert noise in projs
        end
        if(saveinfo)
            save(['res',filesep,'Shepp-Logan',filesep,'proj.mat'],'dataProj')
            infoDicom = [];
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

nIter = [2,3,4,8];           % Iteration to be saved (MLEM or SART)
filterType = 'FBP';          % Filter type: 'BP', 'FBP'
cutoff = 0.75;               % Percentage until cut off frequency (FBP)

%% Reconstruction methods

%                       ## Uncomment to use ##
dataRecon3d = FBP(dataProj,filterType,cutoff,parameter);
if(saveinfo)
    saveData(dataRecon3d,parameter,answer,[],infoDicom);
end

%                       ## Uncomment to use ##
% dataRecon3d = MLEM(dataProj,nIter,parameter);
% if(saveinfo)
%     saveData(dataRecon3d,parameter,answer,[],infoDicom);
% end

%                       ## Uncomment to use ##
% dataRecon3d = SART(dataProj,nIter,parameter);
% if(saveinfo)
%     saveData(dataRecon3d,parameter,answer,[],infoDicom);
% end

fprintf('Finished \n');
fprintf('Please, check "res" folder for results \n');

%% Show info

if(animation && data == 3)
    % Create a figure of screen size
    figureScreenSize()
    % 3D Backprojection Visualization
    for k=1:parameter.nz
        subplot(1,2,1)
        imshow(data3d(:,:,k))
        title(['Slice Orig ',num2str(k)]);axis on;
        subplot(1,2,2)
        imshow(dataRecon3d{1,end}(:,:,k),[])   
        title(['Slice Recon ',num2str(k)]);axis on;
        pause(0.08)
    end
end