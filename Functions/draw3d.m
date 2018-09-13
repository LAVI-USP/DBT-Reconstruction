%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% 
%     DESCRIPTION:
%     Draw x-ray source, axis, projection plane and object plane.
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
%% Draw 3D animation
function draw3d(xCoord,yCoord,zCoord,puCoord,pvCoord,param,projnum,theta)

sb = subplot(2,2,[1 3]);    % Create a subplot and merge 1 and 3
cla(sb)                     % Clear current subplot window 

% Detector Coordinate sytem in (mm)
[uCoord,vCoord] = meshgrid(param.us,param.vs);

axis([ param.us(1,end),param.us(1,1)+10,param.vs(1,1)-10,param.vs(1,end)+10,-5,1.1*param.DSD]);
xlabel(['x(mm)'],'FontSize', 20);% / ', num2str(param.nx),'(Voxel Unity)'])
ylabel(['y(mm)'],'FontSize', 20);% / ', num2str(param.ny),'(Voxel Unity)'])
zlabel(['z(mm)'],'FontSize', 20);%;% / ', num2str(param.nz),'(Voxel Unity)'])
title(num2str(projnum))
view([94.6,27.4])
grid on
hold on

% Detector plane vertices coordinates
vetU = [uCoord(1,1),uCoord(end,1),uCoord(end,end),uCoord(1,end)];
vetV = [vCoord(1,1),vCoord(end,1),vCoord(end,end),vCoord(1,end)];

% Projection plane vertices coordinates
vetPu = [puCoord(1,1),puCoord(end,1),puCoord(end,end),puCoord(1,end)];
vetPv = [pvCoord(1,1),pvCoord(end,1),pvCoord(end,end),pvCoord(1,end)];

% Object volume vertices coordinates
vetObj =  [ xCoord(1,1)     yCoord(1,1)       zCoord(1,1)
            xCoord(end,1)   yCoord(end,1)     zCoord(1,1)
            xCoord(end,end) yCoord(end,end)   zCoord(1,1)
            xCoord(1,end)   yCoord(1,end)     zCoord(1,1)
            xCoord(1,1)     yCoord(1,1)       zCoord(end,1)
            xCoord(end,1)   yCoord(end,1)     zCoord(end,1)
            xCoord(end,end) yCoord(end,end)   zCoord(end,1)
            xCoord(1,end)   yCoord(1,end)     zCoord(end,1)];

% Link of vertices for each face
vetor_faces =  [1 2 3 4
                5 6 7 8
                1 2 6 5
                3 4 8 7
                1 4 8 5
                2 6 7 3];
            
%% Draw Everything

% Draw the world's coordinate system to use as reference
%   - Colors: x in red, y in green, z in blue.
%drawAxes([0,0,0,0,0,0], [0.15*param.su, 0.15*param.su, 0.1*param.sz], 2, ['r','g','b']);

% Draw the sensor coordinate system to use as reference
%   - Colors: x in red, y in green, z in blue.
%drawAxes([param.us(1),param.vs(1),0,0,0,0], [-0.15*param.su, 0.15*param.su, 0], 2, ['r','g','b']);

% Rotate x-ray source
xSource = 0;
ySource = double(param.DSR.*sin(-theta));
zSource = double(param.DSR.*cos(theta)+param.DDR);
% Draw x-ray source
plot3(xSource,ySource,zSource,'r.','markersize',45)
% Draw text
text(xSource,ySource-.15*abs(ySource),1.04*zSource,'X-Ray Source','FontSize', 16)

% Draw Detector plane
patch(vetU,vetV,[0,0,0,0],'FaceColor','none','EdgeColor','black','linewidth',2,'LineStyle','--')

% Draw Projection plane
patch(vetPu,vetPv,[0,0,0,0],'yellow','linewidth',1)

% Draw Object volume
patch('Vertices',vetObj,'Faces',vetor_faces,'FaceColor', 'blue','FaceAlpha',0.6)

% Draw lines from the X-ray source to the projection plane passing 
% through the upper slice of the object
plot3([xSource,xCoord(1,1),puCoord(1,1)],[ySource,yCoord(1,1),pvCoord(1,1)],[zSource,zCoord(end,1),0],'--')
plot3([xSource,xCoord(end,1),puCoord(end,1)],[ySource,yCoord(end,1),pvCoord(end,1)],[zSource,zCoord(end,1),0],'--')
plot3([xSource,xCoord(1,end),puCoord(1,end)],[ySource,yCoord(1,end),pvCoord(1,end)],[zSource,zCoord(end,1),0],'--')
plot3([xSource,xCoord(end,end),puCoord(end,end)],[ySource,yCoord(end,end),pvCoord(end,end)],[zSource,zCoord(end,1),0],'--')

drawnow

end