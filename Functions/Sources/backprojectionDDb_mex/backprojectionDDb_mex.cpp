/*
%% Author: Rodrigo de Barros Vimieiro
% Date: March, 2019
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 backprojectionDDb_mex(proj,param)
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
%% 3-D Back-projection Branchless Distance Driven Code (CPU-Multithread)
*/

#include "backprojectionDDb_mex.h"

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

}


// CPU Branchless Distance-driven back-projection
void backprojectionDDb(double* const pVolume,
	double* const pProj,
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

	mexPrintf("CPU running with maximum number of thrreads: %d \n", nThreads);
	
	clock_t time;

	time = clock();
	
	// Number of mapped detector and pixels
	const int nDetXMap = nDetX + 1;
	const int nDetYMap = nDetY + 1;
	const int nPixXMap = nPixX + 1;
	const int nPixYMap = nPixY + 1;


	// Allocate memory
	double* const pSliceI = (double*)mxMalloc(nPixXMap*nPixYMap * sizeof(double));


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
	double* const pProjt = (double*)mxMalloc(nDetYMap *nDetXMap * nProj * sizeof(double));
	double* const pVolumet = (double*)mxMalloc(nPixY *nPixX * nSlices * sizeof(double));

	
	// Initiate variables value with 0
	for (int nz = 0; nz < nSlices; nz++)
		for (int x = 0; x < nPixX; x++)
			for (int y = 0; y < nPixY; y++)
				pVolumet[(nz*nPixY *nPixX) + (x*nPixY) + y] = 0;
	

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
	for (int p = 0; p < nProj; p++) {

		for (int y = 0; y < nDetYMap; y++)
			pProjt[(p*nDetYMap *nDetXMap) + y] = 0;
		
		for (int x = 1; x < nDetXMap; x++)
			pProjt[(p*nDetYMap *nDetXMap) + (x*nDetYMap)] = 0;
	}
		

	// Integrate on I direction
	for (int p = 0; p < nProj; p++) 
		#pragma omp parallel for num_threads(nThreads) schedule(runtime)
		for (int x = 0; x < nDetX; x++){
			double sum = 0;
			for (int y = 0; y < nDetY; y++){
				sum += pProj[(p*nDetY *nDetX) + (x*nDetY) + y];
				pProjt[(p*nDetYMap *nDetXMap) + ((x+1)*nDetYMap) + y + 1] = sum;
			}
		}
	
	// Integrate on J direction
	for (int p = 0; p < nProj; p++) 	
		#pragma omp parallel for num_threads(nThreads) schedule(runtime)
		for (int y = 1; y < nDetYMap; y++) 
			for (int x = 2; x < nDetXMap; x++) 
				pProjt[(p*nDetYMap *nDetXMap) + (x*nDetYMap) + y] += pProjt[(p*nDetYMap *nDetXMap) + ((x - 1)*nDetYMap) + y];

	
	double* pDetmX_tmp = pDetmX + (nDetYMap * (nDetXMap - 2));

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
								  
			bilinear_interpolation(pSliceI, pProjt, pObjX, pObjY, pDetmX_tmp, pDetmY, nPixXMap, nPixYMap, nDetXMap, nDetYMap, nDetX, nDetY, p);

			/*
			S.3. Differentiation - Eq. 24 - Liu et al (2017)
			*/
						  
			differentiation(pVolumet, pSliceI, tubeX, rtubeY, rtubeZ, pObjX, pObjY, pObjZ, nPixX, nPixY, nPixXMap, nPixYMap, du, dv, dx, dy, dz, nz);

		} // Loop end slices

	} // Loop end Projections

	
	// Copy pProjt to pProj
	for (int nz = 0; nz < nSlices; nz++)
		for (int x = 0; x < nPixX; x++)
			for (int y = 0; y < nPixY; y++)
				pVolume[(nz*nPixX*nPixY) + (x*nPixY) + y] = pVolumet[(nz*nPixX*nPixY) + (x*nPixY) + y] / (double) nProj;
	
	
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
	mxFree(pSliceI);
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
		#pragma omp for schedule(runtime)
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
void bilinear_interpolation(double* pSliceI,
	double* pProj,
	double* pObjX,
	double* pObjY,
	double* pDetmX,
	double* pDetmY,
	const int nPixXMap,
	const int nPixYMap,
	const int nDetXMap,
	const int nDetYMap,
	const int nDetX,
	const int nDetY,
	const unsigned int np){

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

	Matlab code --
	objX = nDetX - (objX ./ detmX(1, end - 1));
	objY = (objY ./ detmX(1, end - 1)) - (detmY(1) / detmX(1, end - 1));

	*/

	#pragma omp parallel num_threads(nThreads) 
	{
		double result_temp1, result_temp2;

		double xNormData, yNormData;
		signed int xData, yData;
		double  alpha, beta;

		double d00, d01, d10, d11;

		#pragma omp for schedule(runtime) 
		for (int x = 0; x < nPixXMap; x++)
			for (int y = 0; y < nPixYMap; y++) {

				// Adjust the mapped coordinates to cross the range of (0-nDetX).*duMap 
				//  Divide by pixelSize to get a unitary pixel size
				const double xNormData = nDetX - pObjX[x] / pDetmX[0];
				const signed int xData = (signed int) floor(xNormData);
				const double  alpha = xNormData - xData;

				// Adjust the mapped coordinates to cross the range of (0-nDetY).*dyMap  
				//  Divide by pixelSize to get a unitary pixel size
				const double yNormData = (pObjY[y] / pDetmX[0]) - (pDetmY[0] / pDetmX[0]);
				const signed int yData = (signed int) floor(yNormData);
				const double  beta = yNormData - yData;

				double d00, d01, d10, d11;
				if (((xNormData) >= 0) && ((xNormData) <= nDetX) && ((yNormData) >= 0) && ((yNormData) <= nDetY))	d00 = pProj[(np*nDetYMap*nDetXMap) + (xData*nDetYMap + yData)];						else    d00 = 0.0;
				if (((xData + 1) > 0) && ((xData + 1) <= nDetX) && ((yNormData) >= 0) && ((yNormData) <= nDetY))	d10 = pProj[(np*nDetYMap*nDetXMap) + ((xData + 1)*nDetYMap + yData)]; 				else    d10 = 0.0;
				if (((xNormData) >= 0) && ((xNormData) <= nDetX) && ((yData + 1) > 0) && ((yData + 1) <= nDetY))	d01 = pProj[(np*nDetYMap*nDetXMap) + (xData*nDetYMap + yData + 1)]; 				else    d01 = 0.0;
				if (((xData + 1) > 0) && ((xData + 1) <= nDetX) && ((yData + 1) > 0) && ((yData + 1) <= nDetY))		d11 = pProj[(np*nDetYMap*nDetXMap) + ((xData + 1)*nDetYMap + yData + 1)];			else    d11 = 0.0;

				double result_temp1 = alpha * d10 + (-d00 * alpha + d00);
				double result_temp2 = alpha * d11 + (-d01 * alpha + d01);

				pSliceI[x * nPixYMap + y] = beta * result_temp2 + (-result_temp1 * beta + result_temp1);
			}
	}

	return;
}

// Differentiation
void differentiation(double* pVolume,
	double* pSliceI,
	double tubeX, 
	double rtubeY, 
	double rtubeZ,
	double* const pObjX, 
	double* const pObjY, 
	double* const pObjZ, 
	const int nPixX,
	const int nPixY,
	const int nPixXMap,
	const int nPixYMap,
	const double du, 
	const double dv, 
	const double dx, 
	const double dy, 
	const double dz, 
	const unsigned int nz){

	/*

	S.3. Differentiation - Eq. 24 - Liu et al (2017)

	Slice integral projection
	___________
	|_A_|_B_|___|
	|_C_|_D_|___|
	|___|___|___|


	(y,x)
	________________
	|_A_|__B__|_____|
	|_C_|(0,0)|(0,1)|
	|___|(1,0)|(1,1)|


	Coordinates on intergal slice:

	A = x * nPixYMap + y
	B = ((x+1) * nPixYMap) + y
	C = x * nPixYMap + y + 1
	D = ((x+1) * nPixYMap) + y + 1

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

		#pragma omp for schedule(runtime)
		for (int x = 0; x < nPixX; x++)
			for (int y = 0; y < nPixY; y++) {

				coordA = x * nPixYMap + y;
				coordB = ((x + 1) * nPixYMap) + y;
				coordC = coordA + 1;
				coordD = coordB + 1;

				// x - ray angle in X coord
				gamma = atan((pObjX[x] + (dx / 2) - tubeX) / (rtubeZ - pObjZ[nz]));
				
				// x - ray angle in Y coord
				alpha = atan((pObjY[y] + (dy / 2) - rtubeY) / (rtubeZ - pObjZ[nz]));
				

				dA = pSliceI[coordA];
				dB = pSliceI[coordB];
				dC = pSliceI[coordC];
				dD = pSliceI[coordD];

				// Treat border of interpolated integral detector
				if (dC == 0 && dD == 0) {
					dC = dA;
					dD = dB;
				}


				// S.3.Differentiation - Eq. 24 - Liu et al(2017)
				pVolume[(nPixX*nPixY*nz) + (x * nPixY) + y] += ((dD - dC - dB + dA)*(du*dv*dz / (cos(alpha)*cos(gamma)*dx*dy)));
			}
	}
	return;
}