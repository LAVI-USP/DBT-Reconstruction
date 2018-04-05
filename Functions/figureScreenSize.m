%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{

    DESCRIPTION:
    Create a figure of 80% of sreen size.

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
%% Create a figure
function  figureScreenSize()
screen_size = get(0,'Screensize');
w = screen_size(3);
h = screen_size(4);
fw = round(w*0.8);
fh = round(h*0.8);
fig_opt = [(w-fw)/2 (h-fh)/2 fw fh];
figure('position', fig_opt)
end
