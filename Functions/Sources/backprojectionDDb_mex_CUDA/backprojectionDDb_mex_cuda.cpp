/*
%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 backprojectionDDb_mex_CUDA(proj,param)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function reconstruct the 3D volume from projections, based on
%     the Distance-Driven principle. It works by calculating the overlap
%     in X and Y axis of the volume and the detector boundaries.
%     The geometry is for DBT with half cone-beam. All parameters are set
%     in "ParameterSettings" code.
%
%     INPUT:
%
%     - proj = 2D projections for each angle
%     - param = Parameter of all geometry
%
%     OUTPUT:
%
%     - data3d = reconstructed volume.
%
%     Reference:
%     - Branchless Distance Driven Projection and Backprojection,
%     Samit Basu and Bruno De Man (2006)
%     - GPU Acceleration of Branchless Distance Driven Projection and
%     Backprojection, Liu et al (2016)
%     - GPU-Based Branchless Distance-Driven Projection and Backprojection,
%     Liu et al (2017)
%     - A GPU Implementation of Distance-Driven Computed Tomography,
%     Ryan D. Wagner (2017)
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
%% 3-D Distance Driven Back-projection Code
*/

#include "backprojectionDDb_mex_cuda.h"

/****************************************************************************
* function: mexFunction() - Matlab interface function.						*
* INPUTS:																	*
*   nlhs - Number of expected output mxArrays, specified as an integer.		*
*   plhs[] - Array of pointers to the expected mxArray output arguments.	*
*   nrhs - Number of input mxArrays, specified as an integer.				*
*   prhs[] - Array of pointers to the mxArray input arguments.				*
****************************************************************************/

void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]) {


	/* Check for proper number of arguments */
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin",
			"backprojection_mex requires two input arguments.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
			"backprojection_mex requires one output argument.");
	}

	/* Check if the input is of proper type */
	if (!mxIsDouble(prhs[0])) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument has to be double.");
	}
	if (!mxIsStruct(prhs[1])) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"Second argument has to be a struct.");
	}

	double* const pProj = (double*)(mxGetPr(prhs[0])); // This variable has to come as single/double
	double* const pTubeAngleD = mxGetPr(mxGetField(prhs[1], 0, "tubeDeg"));
	double* const pDetAngleD = mxGetPr(mxGetField(prhs[1], 0, "detectorDeg"));


	const int unsigned nProj = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nProj"));

	const int unsigned nPixX = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nx"));
	const int unsigned nPixY = (const int)mxGetScalar(mxGetField(prhs[1], 0, "ny"));
	const int unsigned nSlices = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nz"));
	const int unsigned nDetX = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nu"));
	const int unsigned nDetY = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nv"));

	const double dx = (const double)mxGetScalar(mxGetField(prhs[1], 0, "dx"));
	const double dy = (const double)mxGetScalar(mxGetField(prhs[1], 0, "dy"));
	const double dz = (const double)mxGetScalar(mxGetField(prhs[1], 0, "dz"));
	const double du = (const double)mxGetScalar(mxGetField(prhs[1], 0, "du"));
	const double dv = (const double)mxGetScalar(mxGetField(prhs[1], 0, "dv"));

	const double DSD = (const double)mxGetScalar(mxGetField(prhs[1], 0, "DSD"));
	const double DDR = (const double)mxGetScalar(mxGetField(prhs[1], 0, "DDR"));
	const double DAG = (const double)mxGetScalar(mxGetField(prhs[1], 0, "DAG"));

	mwSize NRow = mxGetM(prhs[0]);
	mwSize NCol = mxGetN(prhs[0]);
	mwSize NElemVol = mxGetNumberOfElements(prhs[0]);

	/* Check if the input is in proper size */
	if (NRow != nDetY) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument needs to have the same number of rows as in the configuration file.");
	}

	if (NCol / nProj != nDetX) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument needs to have the same number of columns as in the configuration file.");
	}

	// Cast double vectors to single 
	double* const pTubeAngle = (double*)mxMalloc(nProj * sizeof(double));
	double* const pDetAngle = (double*)mxMalloc(nProj * sizeof(double));

	for (unsigned int n = 0; n < nProj; n++) {
		pTubeAngle[n] = (double)pTubeAngleD[n];
		pDetAngle[n] = (double)pDetAngleD[n];
	}

	//mexPrintf("Nx:%d Ny:%d Nz:%d \nNu:%d Nv:%d \nDx:%.2f Dy:%.2f Dz:%.2f \nDu:%.2f Dv:%.2f \nDSD:%.2f DDR:%.2f \n", nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR);


	// Output memory allocation
	mwSize dims[3] = { nPixY,nPixX,nSlices };

	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	double* const pVolume = (double*)mxGetData(plhs[0]);

	backprojectionDDb(pVolume, pProj, pTubeAngle, pDetAngle, nProj, nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR, DAG);

	mxFree(pTubeAngle);
	mxFree(pDetAngle);

	return;

}