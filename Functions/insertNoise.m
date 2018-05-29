%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projNoisy = insertNoise(proj,param)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function insert noise in the projection images, for virtual
%     phantoms, considering quantum and electronic noise.
%     
%     INPUT:
% 
%     - proj = 2D projection images 
%     - param = Parameter of all geometry
% 
%     OUTPUT:
% 
%     - projNoisy = 2D projections with noise.
% 
%     Reference: ---
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
%% Insert noise on projections
function projNoisy = insertNoise(proj,param)

if(strcmp(param.type,'shepplogan'))
    highestValue = (2^14) - 1;
    
    projNoisy = mat2gray(proj);
    projNoisy = 1-projNoisy;
    projNoisy = projNoisy .* highestValue;
    
    % Decrease projection value due to radiation dose
    projNoisy = projNoisy ./ param.nProj;
    
else
    if(strcmp(param.type,'vct'))       
        projNoisy = proj .* 4.84; 
    end
end
    

for p=1:param.nProj  
    
    % Scaled Possion noise with normal distribution
    poissonNoise = sqrt(param.lambda.*projNoisy(:,:,p)) .* randn([param.nv,param.nu]);

    % Electronic noise
    electronicNoise = sqrt(param.sigmaE) .* randn([param.nv,param.nu]);

    % Build image with noise model
    imageNoisePoissonGaussian = projNoisy(:,:,p) + poissonNoise + electronicNoise;

    projNoisy(:,:,p) = imageNoisePoissonGaussian + param.tau;
end

if(strcmp(param.type,'shepplogan'))
    projNoisy(projNoisy>highestValue) = highestValue;
    projNoisy(projNoisy<0) = 0;
    
    projNoisy = projNoisy ./ highestValue;
    projNoisy = 1-projNoisy;
    projNoisy = projNoisy .* max(proj(:));
end

% figure;imshow(proj(:,:,p),[])
% figure;imshow(projNoisy(:,:,p),[])

end