/*
%% Author: Rodrigo de Barros Vimieiro
% Date: September, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This is the header function
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
%% 3-D Distance Driven Projection Code
*/

//typedef double user_dataType;
//typedef float user_dataType;
//#define user_mexdataType mxDOUBLE_CLASS
//#define user_mexdataType mxSINGLE_CLASS

#include "mex.h"
#define _USE_MATH_DEFINES
#include <math.h>

void mapp2xy(float* const pXmapp, float* const pYmapp, float tubeX, float tubeY, float tubeZ, float* const pXcoord, float* const pYcoord, float* const pZcoord, const int nXelem, const int nYelem);
