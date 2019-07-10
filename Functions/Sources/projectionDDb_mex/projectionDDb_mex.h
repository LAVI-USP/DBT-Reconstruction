/*
%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This is the header function
%     ---------------------------------------------------------------------
%     Copyright (C) <2019>  <Rodrigo de Barros Vimieiro>
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
%% 3-D Projection Branchless Distance Driven Code (CPU-Multithread)
*/

//typedef double user_dataType;
//typedef float user_dataType;
//#define user_mexdataType mxDOUBLE_CLASS
//#define user_mexdataType mxSINGLE_CLASS

#include "mex.h"
#include <omp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>

// Global variable
int nThreads;

void projectionDDb(double* const pProj,
	double* const pVolume,
	double* const pTubeAngle,
	double* const pDetAngle,
	const  int nProj,
	const  int nPixX,
	const  int nPixY,
	const  int nSlices,
	const  int nDetX,
	const  int nDetY,
	const double dx,
	const double dy,
	const double dz,
	const double du,
	const double dv,
	const double DSD,
	const double DDR,
	const double DAG);

void mapBoudaries(double* pBound,
	const int nElem,
	const double valueLeftBound,
	const double sizeElem,
	const double offset);

void mapDet2Slice(double* const pXmapp,
	double* const pYmapp,
	double tubeX,
	double tubeY,
	double tubeZ,
	double * const pXcoord,
	double * const pYcoord,
	double * const pZcoord,
	double ZSlicecoord,
	const int nXelem,
	const int nYelem);

void bilinear_interpolation(double* projI,
	double* pVolume,
	double* pDetmX,
	double* pDetmY,
	const int nDetXMap,
	const int nDetYMap,
	const int nPixXMap,
	const int nPixYMap,
	const double pixelSize,
	const unsigned int nz);

void differentiation(double* pProj,
	double* projI,
	double* const pDetmX,
	double* const pDetmY,
	double tubeX,
	double rtubeY,
	double rtubeZ,
	double* const pDetX,
	double* const pRdetY,
	double* const pRdetZ,
	const int nDetX,
	const int nDetY,
	const int nDetXMap,
	const int nDetYMap,
	const double du,
	const double dv,
	const double dx,
	const double dy,
	const double dz,
	const unsigned int p);
