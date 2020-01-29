/*
%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projectionDDb_mex(data3d,param)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates the volume projection based on the
%     Branchless Distance-Driven principle.
%     The geometry is for DBT with half cone-beam. All parameters are set
%     in "ParameterSettings" code.
%
%     INPUT:
%
%     - data3d = 3D volume for projection
%     - param = Parameter of all geometry
%
%     OUTPUT:
%
%     - proj = projections for each angle.
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
%% 3-D Projection Branchless Distance Driven Code (CPU-Multithread)
*/

#include "projectionDDb_mex.h"


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
			"projection_mex requires two input arguments.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:nargout",
			"projection_mex requires one output argument.");
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

	double* const pVolume = (double*)(mxGetPr(prhs[0])); // This variable has to come as single/double
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
	if (NRow != nPixY) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument needs to have the same number of rows as in the configuration file.");
	}

	if (NCol / nSlices != nPixX) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument needs to have the same number of columns as in the configuration file.");
	}

	// Cast single vectors to double 
	double* const pTubeAngle = (double*)mxMalloc(nProj * sizeof(double));
	double* const pDetAngle = (double*)mxMalloc(nProj * sizeof(double));

	for (unsigned int n = 0; n < nProj; n++) {
		pTubeAngle[n] = (double)pTubeAngleD[n];
		pDetAngle[n] = (double)pDetAngleD[n];
	}

	//mexPrintf("Nx:%d Ny:%d Nz:%d \nNu:%d Nv:%d \nDx:%.2f Dy:%.2f Dz:%.2f \nDu:%.2f Dv:%.2f \nDSD:%.2f DDR:%.2f \n", nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR);


	// Output memory allocation
	mwSize dims[3] = { nDetY,nDetX,nProj };

	plhs[0] = mxCreateNumericArray(3, dims, mxDOUBLE_CLASS, mxREAL);
	double* const pProj = (double*)mxGetData(plhs[0]);

	projectionDDb(pProj, pVolume, pTubeAngle, pDetAngle, nProj, nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR, DAG);

	mxFree(pTubeAngle);
	mxFree(pDetAngle);

	return;

}


// CPU Branchless Distance-driven projection
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
	const double DAG)
{

	nThreads = omp_get_max_threads();

	mexPrintf("CPU running with maximum number of threads: %d \n", nThreads);
	
	clock_t time;

	time = clock();
	
	// Number of mapped detector and pixels
	const int nDetXMap = nDetX + 1;
	const int nDetYMap = nDetY + 1;
	const int nPixXMap = nPixX + 1;
	const int nPixYMap = nPixY + 1;


	// Allocate memory
	double* const projI = (double*)mxMalloc(nDetXMap*nDetYMap * sizeof(double));


	// Pointer for projections coordinates
	// Allocate memory for projections coordinates
	double* const pDetX = (double*)mxMalloc(nDetXMap * sizeof(double));
	double* const pDetY = (double*)mxMalloc(nDetYMap * sizeof(double));
	double* const pDetZ = (double*)mxMalloc(nDetYMap * sizeof(double));
	double* const pObjX = (double*)mxMalloc(nPixXMap * sizeof(double));
	double* const pObjY = (double*)mxMalloc(nPixYMap * sizeof(double));
	double* const pObjZ = (double*)mxMalloc(nSlices * sizeof(double));


	// Pointer for mapped coordinates
	// Allocate memory for mapped coordinates
	double* const pDetmY = (double*)mxMalloc(nDetYMap * sizeof(double));
	double* const pDetmX = (double*)mxMalloc(nDetYMap * nDetXMap * sizeof(double));



	// Pointer for rotated detector coords
	// Allocate memory for rotated detector coords
	double* const pRdetY = (double*)mxMalloc(nDetYMap * sizeof(double));
	double* const pRdetZ = (double*)mxMalloc(nDetYMap * sizeof(double));


	// Map detector and object boudaries
	mapBoudaries(pDetX, nDetXMap, (double)nDetX, -du, 0.0);

	mapBoudaries(pDetY, nDetYMap, nDetY / 2.0, dv, 0.0);

	mapBoudaries(pDetZ, nDetYMap, 0.0, 0.0, 0.0);

	mapBoudaries(pObjX, nPixXMap, (double)nPixX, -dx, 0.0);

	mapBoudaries(pObjY, nPixYMap, nPixY / 2.0, dy, 0.0);

	mapBoudaries(pObjZ, nSlices, 0.0, dz, DAG + (dz / 2.0));


	// X - ray tube initial position
	double tubeX = 0;
	double tubeY = 0;
	double tubeZ = DSD;

	// Iso - center position
	double isoY = 0;
	double isoZ = DDR;


	// Allocate memory for temp projection variable
	double* const pProjt = (double*)mxMalloc(nDetY *nDetX * nProj * sizeof(double));
	double* const pVolumet = (double*)mxMalloc(nPixYMap *nPixXMap * nSlices * sizeof(double));

	
	// Initiate variables value with 0
	for (int p = 0; p < nProj; p++)
		for (int x = 0; x < nDetX; x++)
			for (int y = 0; y < nDetY; y++)
				pProjt[(p*nDetY *nDetX) + (x*nDetY) + y] = 0;
	

	// Integration of 2D slices over the whole volume
	// (S.1.Integration. - Liu et al(2017))
	/*
		
		2 - D image

		 -->J
		|	-------------
		v	|			|
		I	|			|
			|			|
			|			|
			-------------
	*/
	
	
	// Initialize first column and row with zeros
	for (int nz = 0; nz < nSlices; nz++) {

		for (int y = 0; y < nPixYMap; y++)
			pVolumet[(nz*nPixYMap *nPixXMap) + y] = 0;
		
		for (int x = 1; x < nPixXMap; x++)
			pVolumet[(nz*nPixYMap *nPixXMap) + (x*nPixYMap)] = 0;
	}
		

	// Integrate on I direction
	for (int nz = 0; nz < nSlices; nz++) 
		#pragma omp parallel for num_threads(nThreads) schedule(static)
		for (int x = 0; x < nPixX; x++){
			double sum = 0;
			for (int y = 0; y < nPixY; y++){
				sum += pVolume[(nz*nPixY *nPixX) + (x*nPixY) + y];
				pVolumet[(nz*nPixYMap *nPixXMap) + ((x+1)*nPixYMap) + y + 1] = sum;
			}
		}

	// Integrate on J direction
	for (int nz = 0; nz < nSlices; nz++) 	
		#pragma omp parallel for num_threads(nThreads) schedule(static)
		for (int y = 1; y < nPixYMap; y++) 
			for (int x = 2; x < nPixXMap; x++) 
				pVolumet[(nz*nPixYMap *nPixXMap) + (x*nPixYMap) + y] += pVolumet[(nz*nPixYMap *nPixXMap) + ((x - 1)*nPixYMap) + y];	



	// For each projection
	for (int p = 0; p < nProj; p++) {
		
		// Get specif tube angle for the projection
		double theta = pTubeAngle[p] * M_PI / 180.0;

		// Get specif detector angle for the projection
		double phi = pDetAngle[p] * M_PI / 180.0;

		//mexPrintf("Tube angle:%f Det angle:%f\n", theta, phi);

		// Tube rotation
		double rtubeY = ((tubeY - isoY)*cos(theta) - (tubeZ - isoZ)*sin(theta)) + isoY;
		double rtubeZ = ((tubeY - isoY)*sin(theta) + (tubeZ - isoZ)*cos(theta)) + isoZ;

		//mexPrintf("R tube Y:%f R tube Z:%f\n", rtubeY, rtubeZ);

		// Detector rotation
		for (int y = 0; y < nDetYMap; y++) {
			pRdetY[y] = ((pDetY[y] - isoY)*cos(phi) - (pDetZ[y] - isoZ)*sin(phi)) + isoY;
			pRdetZ[y] = ((pDetY[y] - isoY)*sin(phi) + (pDetZ[y] - isoZ)*cos(phi)) + isoZ;
		}


		// For each slice
		for (int nz = 0; nz < nSlices; nz++) {
			
			/*

			Map detector onto XY plane(Inside proj loop in case detector rotates)

			*** Note: Matlab has linear indexing as a column of elements, i.e, the elements are actually stored in memory as queued columns.
	
			*/

			// Map slice onto XY plane
			mapDet2Slice(pDetmX, pDetmY, tubeX, rtubeY, rtubeZ, pDetX, pRdetY, pRdetZ, pObjZ[nz], nDetXMap, nDetYMap);

			/*
			S.2. Interpolation - Liu et al (2017)
			*/

			bilinear_interpolation(projI, pVolumet, pDetmX, pDetmY, nDetXMap, nDetYMap, nPixXMap, nPixYMap, dx, nz);

			/*
			S.3. Differentiation - Eq. 24 - Liu et al (2017)
			*/

			differentiation(pProjt, projI, pDetmX, pDetmY, tubeX, rtubeY, rtubeZ, pDetX, pRdetY, pRdetZ, nDetX, nDetY, nDetXMap, nDetYMap, du, dv, dx, dy, dz, p);

		} // Loop end slices

	} // Loop end Projections

	
	// Copy pProjt to pProj
	for (int p = 0; p < nProj; p++)
		for (int x = 0; x < nDetX; x++)
			for (int y = 0; y < nDetY; y++)
				pProj[(p*nDetX*nDetY) + (x*nDetY) + y] = pProjt[(p*nDetX*nDetY) + (x*nDetY) + y];
	
	
	// Interp
	//for (int x = 0; x < nDetXMap; x++)
	//	for (int y = 0; y < nDetYMap; y++) {
	//		pProj[(0 * nDetXMap*nDetYMap) + (x*nDetYMap) + y] = projI[(x*nDetYMap) + y];
	//	}

	// Volume			
	//for (int p = 0; p < nSlices; p++)
	//	for (int x = 0; x < nPixXMap; x++)
	//		for (int y = 0; y < nPixYMap+1; y++) {
	//			pProj[(p * nPixXMap*nPixYMap) + (x*nPixYMap) + y] = pVolumet[(p * nPixXMap*nPixYMap) + (x*nPixYMap) + y];
	//		}
	
	mxFree(pProjt);
	mxFree(projI);
	mxFree(pVolumet);
	mxFree(pDetX);
	mxFree(pDetY);
	mxFree(pDetZ);
	mxFree(pObjX);
	mxFree(pObjY);
	mxFree(pObjZ);
	mxFree(pDetmY);
	mxFree(pDetmX);
	mxFree(pRdetY);
	mxFree(pRdetZ);
	
	time = clock() - time;
	mexPrintf("Processing time: %f seconds.\n", ((double)time) / CLOCKS_PER_SEC);

	return;	
}


// Make boundaries of detector and slices
void mapBoudaries(double* pBound, 
	const int nElem, 
	const double valueLeftBound, 
	const double sizeElem, 
	const double offset){

	for (int k = 0; k < nElem; k++)
		pBound[k] = (k - valueLeftBound) * sizeElem + offset;

	return;
}

// Map on XY plane
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
	const int nYelem){

	#pragma omp parallel num_threads(nThreads) 
	{
		int ind;
		#pragma omp for schedule(static)
			for (int x = 0; x < nXelem; x++)
				for (int y = 0; y < nYelem; y++) {

					ind = (x*nYelem) + y;
					pXmapp[ind] = ((pXcoord[x] - tubeX)*(ZSlicecoord - pZcoord[y]) - (pXcoord[x] * tubeZ) + (pXcoord[x] * pZcoord[y])) / (-tubeZ + pZcoord[y]);

					if (x == 0)
						pYmapp[y] = ((pYcoord[y] - tubeY)*(ZSlicecoord - pZcoord[y]) - (pYcoord[y] * tubeZ) + (pYcoord[y] * pZcoord[y])) / (-tubeZ + pZcoord[y]);
				}
	}

	return;
}

// Bilinear interpolation
void bilinear_interpolation(double* projI, 
	double* pVolume, 
	double* pDetmX, 
	double* pDetmY, 
	const int nDetXMap, 
	const int nDetYMap,
	const int nPixXMap, 
	const int nPixYMap, 
	const double pixelSize,
	const unsigned int nz) {


	/*

	S.2. Interpolation - Liu et al (2017)

	Reference:
	- https://en.wikipedia.org/wiki/Bilinear_interpolation
	- https://stackoverflow.com/questions/21128731/bilinear-interpolation-in-c-c-and-cuda

	Note: *** We are using the Unit Square equation ***

	alpha = X - X1
	1-alpha = X2 - X

	beta = Y - Y1
	1-beta = Y2 - Y

	----> Xcoord 
	|				0___________
	v				|_d00_|_d10_|
	Ycoord			|_d01_|_d11_|

	*/

	// These are the boundaries of the slices, ranging from (0-nPixX)
	const double PixXbound = (double) (nPixXMap - 1);
	const double PixYbound = (double) (nPixYMap - 1);

	#pragma omp parallel for num_threads(nThreads) schedule(static) 
	for (int x = 0; x < nDetXMap; x++)
		for (int y = 0; y < nDetYMap; y++) {

			// Adjust the mapped coordinates to cross the range of (0-nPixX).*dx 
			// Divide by pixelSize to get a unitary pixel size
			double xNormData = PixXbound - pDetmX[x * nDetYMap + y] / pixelSize;
			signed int xData = (signed int)floor(xNormData);
			double alpha = xNormData - xData;

			// Adjust the mapped coordinates to cross the range of (0-nPixY).*dy  
			// Divide by pixelSize to get a unitary pixel size
			double yNormData = (PixYbound / 2.0) + (pDetmY[y] / pixelSize);
			signed int yData = (signed int)floor(yNormData);
			double beta = yNormData - yData;

			double d00, d01, d10, d11;
			if (((xNormData) >= 0) && ((xNormData) <= PixXbound) && ((yNormData) >= 0) && ((yNormData) <= PixYbound))	d00 = pVolume[(nz*nPixYMap*nPixXMap) + (xData*nPixYMap + yData)];					else    d00 = 0.0;
			if (((xData + 1) > 0) && ((xData + 1) <= PixXbound) && ((yNormData) >= 0) && ((yNormData) <= PixYbound))	d10 = pVolume[(nz*nPixYMap*nPixXMap) + ((xData + 1)*nPixYMap + yData)]; 			else    d10 = 0.0;
			if (((xNormData) >= 0) && ((xNormData) <= PixXbound) && ((yData + 1) > 0) && ((yData + 1) <= PixYbound))	d01 = pVolume[(nz*nPixYMap*nPixXMap) + (xData*nPixYMap + yData + 1)]; 				else    d01 = 0.0;
			if (((xData + 1) > 0) && ((xData + 1) <= PixXbound) && ((yData + 1) > 0) && ((yData + 1) <= PixYbound))		d11 = pVolume[(nz*nPixYMap*nPixXMap) + ((xData + 1)*nPixYMap + yData + 1)];			else    d11 = 0.0;

			double result_temp1 = alpha * d10 + (-d00 * alpha + d00);
			double result_temp2 = alpha * d11 + (-d01 * alpha + d01);

			projI[x * nDetYMap + y] = beta * result_temp2 + (-result_temp1 * beta + result_temp1);

		}
		
	return;
}

// Differentiation
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
	const unsigned int p) {


	/*

	S.3. Differentiation - Eq. 24 - Liu et al (2017)

	Detector integral projection
	___________
	|_A_|_B_|___|
	|_C_|_D_|___|
	|___|___|___|


	(y,x)
	________________
	|_A_|__B__|_____|
	|_C_|(0,0)|(0,1)|
	|___|(1,0)|(1,1)|


	Coordinates on intergal projection:

	A = x * nDetYMap + y
	B = ((x+1) * nDetYMap) + y
	C = x * nDetYMap + y + 1
	D = ((x+1) * nDetYMap) + y + 1

	*/

	#pragma omp parallel num_threads(nThreads)
	{
		unsigned int coordA;
		unsigned int coordB;
		unsigned int coordC;
		unsigned int coordD;

		double detMapSizeX;
		double detMapSizeY;
		double gamma;
		double alpha;

		double dA, dB, dC, dD;

		#pragma omp for schedule(static)
		for (int x = 0; x < nDetX; x++)
			for (int y = 0; y < nDetY; y++) {

				coordA = x * nDetYMap + y;
				coordB = ((x + 1) * nDetYMap) + y;
				coordC = coordA + 1;
				coordD = coordB + 1;


				// Detector X coord overlap calculation
				detMapSizeX = pDetmX[coordA] - pDetmX[coordB];

				// Detector Y coord overlap calculation
				detMapSizeY = pDetmY[y + 1] - pDetmY[y];

				// x - ray angle in X coord
				gamma = atan((pDetX[x] + (du / 2) - tubeX) / rtubeZ);

				// x - ray angle in Y coord
				alpha = atan((pRdetY[y] + (dv / 2) - rtubeY) / rtubeZ);


				dA = projI[coordA];
				dB = projI[coordB];
				dC = projI[coordC];
				dD = projI[coordD];

				// Treat border of interpolated integral detector
				if (dC == 0 && dD == 0) {
					dC = dA;
					dD = dB;
				}


				// S.3.Differentiation - Eq. 24 - Liu et al(2017)
				pProj[(nDetX*nDetY*p) + (x * nDetY) + y] += ((dD - dC - dB + dA)*(dx*dy*dz / (cos(alpha)*cos(gamma)*detMapSizeX*detMapSizeY)));
			}
	}
	return;
}