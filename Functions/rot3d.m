%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                      rot3d(X,Y,Z, theta, aroundAxis)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     Counter-clockwise rotation around axis with positive angles.
%  
%  
%     INPUT:
% 
%     - X,Y,Z = Vectors with coordinates to be rotated
%     - theta = Rotation angle
%     - aroundAxis = Around which axis
% 
%     OUTPUT:
% 
%     - Xr,Yr,Zr = rotated coordinates
% 
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
%% 3D Rotation Matrix Code
function [Xr,Yr,Zr] = rot3d(X,Y,Z, theta, aroundAxis)

[m,n] = size(Z);

Xr = zeros(m,n);
Yr = Xr;
Zr = Xr;

M = ones(4,numel(Z));

k=1;
for i=1:1:m
    for j=1:1:n
        M(1,k) = X(i,j);
        M(2,k) = Y(i,j); 
        M(3,k) = Z(i,j); 
        k=k+1;
    end
end

switch aroundAxis
    case 'x'
        % Rotation matrix 3D around X axis
        mRotate = [  1       0           0      0
                     0   cos(theta)  -sin(theta)  0
                     0   sin(theta)   cos(theta)  0
                     0       0           0      1];
    case 'y'
        % Rotation matrix 3D around Y axis
        mRotate = [ cos(theta)   0  sin(theta)   0
                     0          1       0      0
                   -sin(theta)   0  cos(theta)   0
                     0          0       0      1];
    case 'z'
        % Rotation matrix 3D around Z axis
        mRotate = [ cos(theta)   -sin(theta)   0  0
                    sin(theta)    cos(theta)   0  0
                       0            0        1  0
                       0            0        0  1];       
end
Mrotated = (mRotate*M);

k=1;
for i=1:1:m
    for j=1:1:n
        Xr(i,j) = Mrotated(1,k);
        Yr(i,j) = Mrotated(2,k);
        Zr(i,j) = Mrotated(3,k);
        k=k+1;
    end
end

end