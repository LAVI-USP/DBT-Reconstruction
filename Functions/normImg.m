%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
---------------------------------------------------------------------------
                            normImg(data)
---------------------------------------------------------------------------
    DESCRIPTION:

    Normalize image between 0 and 1.

    INPUT:

    - data = Data to be normalized 

    OUTPUT:

    - dataNorm = Normalized data.
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
function dataNorm = normImg(data)

minData = min(data(:));

dataNorm = data - minData;

maxData = max(dataNorm(:));

dataNorm = dataNorm ./ maxData;

end
