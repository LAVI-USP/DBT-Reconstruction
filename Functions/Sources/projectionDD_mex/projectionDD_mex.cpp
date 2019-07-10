/*
%% Author: Rodrigo de Barros Vimieiro
% Date: September, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 projectionDD_mex(data3d,param)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function calculates the volume projection based on the
%     Distance-Driven principle. It works by calculating the overlap in X
%     and Y axis of the volume and the detector boundaries.
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
%     Reference: Three-Dimensional Digital Tomosynthesis - Yulia
%     Levakhina (2014), Cap 3.6 and 3.7.
%
%     Original Paper: De Man, Bruno, and Samit Basu. "Distance-driven
%     projection and backprojection in three dimensions." Physics in
%     Medicine & Biology (2004).
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
%% 3-D Distance Driven Projection Code
*/

#include "projectionDD_mex.h"

/****************************************************************************
* function: mexFunction() - Matlab interface function.						*
* INPUTS:																	*
*   nlhs - Number of expected output mxArrays, specified as an integer.		*
*   plhs[] - Array of pointers to the expected mxArray output arguments.	*
*   nrhs - Number of input mxArrays, specified as an integer.				*
*   prhs[] - Array of pointers to the mxArray input arguments.				*
****************************************************************************/

void mexFunction(int nlhs, mxArray* plhs[],
	int nrhs, const mxArray* prhs[]) {


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
	if (!mxIsSingle(prhs[0])) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"First argument has to be single.");
	}
	if (!mxIsStruct(prhs[1])) {
		mexErrMsgIdAndTxt("MATLAB:mexcpp:typeargin",
			"Second argument has to be a struct.");
	}


	float* pVolume = (float*)mxGetPr(prhs[0]);
	double* pTubeAngleD = mxGetPr(mxGetField(prhs[1], 0, "tubeDeg"));
	double* pDetAngleD = mxGetPr(mxGetField(prhs[1], 0, "detectorDeg"));

	const int nProj = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nProj"));

	const int nPixX = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nx"));
	const int nPixY = (const int)mxGetScalar(mxGetField(prhs[1], 0, "ny"));
	const int nSlices = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nz"));
	const int nDetX = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nu"));
	const int nDetY = (const int)mxGetScalar(mxGetField(prhs[1], 0, "nv"));

	const float dx = (const float)mxGetScalar(mxGetField(prhs[1], 0, "dx"));
	const float dy = (const float)mxGetScalar(mxGetField(prhs[1], 0, "dy"));
	const float dz = (const float)mxGetScalar(mxGetField(prhs[1], 0, "dz"));
	const float du = (const float)mxGetScalar(mxGetField(prhs[1], 0, "du"));
	const float dv = (const float)mxGetScalar(mxGetField(prhs[1], 0, "dv"));

	const float DSD = (const float)mxGetScalar(mxGetField(prhs[1], 0, "DSD"));
	const float DDR = (const float)mxGetScalar(mxGetField(prhs[1], 0, "DDR"));
	const float DAG = (const float)mxGetScalar(mxGetField(prhs[1], 0, "DAG"));

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

	// Cast double vectors to single 
	float* const pTubeAngle = (float*)mxMalloc(nProj * sizeof(float));
	float* const pDetAngle = (float*)mxMalloc(nProj * sizeof(float));

	for (int n = 0; n < nProj; n++) {
		pTubeAngle[n] = (float)pTubeAngleD[n];
		pDetAngle[n] = (float)pDetAngleD[n];
	}

	//mexPrintf("Nx:%d Ny:%d Nz:%d \nNu:%d Nv:%d \nDx:%.2f Dy:%.2f Dz:%.2f \nDu:%.2f Dv:%.2f \nDSD:%.2f DDR:%.2f \n", nPixX, nPixY, nSlices, nDetX, nDetY, dx, dy, dz, du, dv, DSD, DDR);


	const int nDetXMap = nDetX + 1;
	const int nDetYMap = nDetY + 1;
	const int nPixXMap = nPixX + 1;
	const int nPixYMap = nPixY + 1;

	// Memory allocation for projections coordinates
	float* const pDetX = (float*)mxMalloc(nDetXMap * sizeof(float));
	float* const pDetY = (float*)mxMalloc(nDetYMap * sizeof(float));
	float* const pDetZ = (float*)mxMalloc(nDetYMap * sizeof(float));
	float* const pObjX = (float*)mxMalloc(nPixXMap * sizeof(float));
	float* const pObjY = (float*)mxMalloc(nPixYMap * sizeof(float));
	float* const pObjZ = (float*)mxMalloc(nSlices * sizeof(float));

	// Memory allocation for mapped coordinates
	float* const pDetmY = (float*)mxMalloc(nDetYMap * sizeof(float));
	float* const pDetmX = (float*)mxMalloc(nDetYMap * nDetXMap * sizeof(float));
	float* const pPixmX = (float*)mxMalloc(nPixXMap * sizeof(float));
	float* const pPixmY = (float*)mxMalloc(nPixYMap * sizeof(float));


	// Memory allocation for rotated detecto coords
	float* const pRdetY = (float*)mxMalloc(nDetYMap * sizeof(float));
	float* const pRdetZ = (float*)mxMalloc(nDetYMap * sizeof(float));

	// Memory allocation for slice
	float* const pSlice = (float*)mxMalloc(nPixX * nPixY * sizeof(float));


	// Map detector and object boudaries
	for (int k = nDetX, ind = 0; k >= 0; k--, ind++)					pDetX[ind] = k * du;
	for (int k = -nDetY / 2, ind = 0; k <= nDetY / 2; k++, ind++)		pDetY[ind] = k * dv;
	for (int k = 0; k <= nDetY; k++)								pDetZ[k] = 0;
	for (int k = nPixX, ind = 0; k >= 0; k--, ind++)				pObjX[ind] = k * dx;
	for (int k = -nPixY / 2, ind = 0; k <= nPixY / 2; k++, ind++)	pObjY[ind] = k * dy;
	for (int k = 0, ind = 0; k <= nSlices; k++, ind++)				pObjZ[ind] = k * dz + DAG + dz / 2;


	// X - ray tube initial position
	float tubeX = 0;
	float tubeY = 0;
	float tubeZ = DSD;

	// Iso - center position
	float isoY = 0;
	float isoZ = DDR;

	// Allocate memory for temp projection variable
	float* const pProjt = (float*)mxMalloc(nDetY * nDetX * nProj * sizeof(float));

	// Initiate temp proj values with zeros
	for (int p = 0; p < nProj; p++)
		for (int u = 0; u < nDetX; u++)
			for (int v = 0; v < nDetY; v++)
				pProjt[(p * nDetY * nDetX) + (u * nDetY) + v] = 0;

	// For each projection
	for (int p = 0; p < nProj; p++) {

		// Get specif tube angle for the projection
		float theta = pTubeAngle[p] * M_PI / 180.f;

		// Get specif detector angle for the projection
		float phi = pDetAngle[p] * M_PI / 180.f;


		// Tubre rotation
		float rtubeY = ((tubeY - isoY) * (float)cos(theta) - (tubeZ - isoZ) * (float)sin(theta)) + isoY;
		float rtubeZ = ((tubeY - isoY) * (float)sin(theta) + (tubeZ - isoZ) * (float)cos(theta)) + isoZ;

		// Detector rotation
		for (int v = 0; v < nDetYMap; v++) {
			pRdetY[v] = ((pDetY[v] - isoY) * (float)cos(phi) - (pDetZ[v] - isoZ) * (float)sin(phi)) + isoY;
			pRdetZ[v] = ((pDetY[v] - isoY) * (float)sin(phi) + (pDetZ[v] - isoZ) * (float)cos(phi)) + isoZ;
		}


		// Map detector onto XY plane(Inside proj loop in case detector rotates)
		mapp2xy(pDetmX, pDetmY, tubeX, rtubeY, rtubeZ, pDetX, pRdetY, pRdetZ, nDetXMap, nDetYMap);


		// Pixel start index and increment
		int detIstart = 0;
		int detIinc = 1;

		// Mapped detector length
		float deltaDetmY = pDetmY[detIstart + detIinc] - pDetmY[detIstart];


		// For each slice
		for (int z = 0; z < nSlices; z++) {

			// Flip X(Img coord is reverse to Global)
			for (int x = 0, x_inv = nPixX - 1; x < nPixX; x++, x_inv--)
				for (int y = 0; y < nPixY; y++)
					pSlice[(x * nPixY) + y] = pVolume[(z * nPixX * nPixY) + (x_inv * nPixY) + y];

			// Tmp Z coords value for Y direction
			float* const pObjZt = (float*)mxMalloc(nPixYMap * sizeof(float));

			// Tmp Pixel X mapped coords
			float* const pPixmXt = (float*)mxMalloc(nPixYMap * nPixXMap * sizeof(float));

			// Get specif Z coord value for each slice
			for (int k = 0; k < nPixYMap; k++)	pObjZt[k] = pObjZ[z];

			// Map slice onto XY plane
			mapp2xy(pPixmXt, pPixmY, tubeX, rtubeY, rtubeZ, pObjX, pObjY, pObjZt, nPixXMap, nPixYMap);

			// Flip X(Img coord is reverse to Global)
			for (int x = 0, x_inv = nPixXMap - 1; x < nPixXMap; x++, x_inv--)
				pPixmX[x] = pPixmXt[x_inv * nPixYMap];

			// Free temp variables
			mxFree(pObjZt);
			mxFree(pPixmXt);

			// Pixel start index and increment
			int pixIstart = 0;
			int pixIinc = 1;

			// Mapped pixel length
			float deltaPixmX = pPixmX[pixIstart + pixIinc] - pPixmX[pixIstart];
			float deltaPixmY = pPixmY[pixIstart + pixIinc] - pPixmY[pixIstart];

			// Start pixel and detector indices
			int detIndY = detIstart;
			int pixIndY = pixIstart;

			// Case 1
			// Find first detector overlap maped with pixel maped on Y
			if (pDetmY[detIndY] - pPixmY[pixIstart] < -deltaDetmY)
				while (pDetmY[detIndY] - pPixmY[pixIstart] < -deltaDetmY)
					detIndY = detIndY + detIinc;

			else
				// Case 2
				// Find first pixel overlap maped with detector maped on Y
				if (pDetmY[detIstart] - pPixmY[pixIndY] > deltaPixmY)
					while (pDetmY[detIstart] - pPixmY[pixIndY] > deltaPixmY)
						pixIndY = pixIndY + pixIinc;

			float moving_left_boundaryY = NULL;

			// Get the left coordinate of the first overlap on Y axis
			if (pDetmY[detIndY] < pPixmY[pixIndY])
				moving_left_boundaryY = pPixmY[pixIndY];
			else
				moving_left_boundaryY = pDetmY[detIndY];


			// Allocate memory for specif row of X map detector coords
			float* const pDetmXrow = (float*)mxMalloc(nDetXMap * sizeof(float));

			float overLapY = NULL;

			// Loop over Y intersections
			while ((detIndY < nDetY) && (pixIndY < nPixY)) {

				float alpha = (float)atan((pDetmY[detIndY] + (deltaDetmY / 2) - rtubeY) / rtubeZ);

				// Case A, when you jump to the next detector boundarie but stay
				// in the same pixel
				if (pDetmY[detIndY + 1] <= pPixmY[pixIndY + 1])
					overLapY = (pDetmY[detIndY + 1] - moving_left_boundaryY) / deltaDetmY; // Normalized overlap Calculation

				else
					// Case B, when you jump to the next pixel boundarie but stay
					// in the same detector
					overLapY = (pPixmY[pixIndY + 1] - moving_left_boundaryY) / deltaDetmY; // Normalized overlap Calculation

				//										***** X overlap *****
				int detIndX = detIstart;
				int pixIndX = pixIstart;


				// Get row / coll of X flipped, which correspond to that Y overlap det
				for (int x = 0, x_inv = nDetXMap - 1; x < nDetXMap; x++, x_inv--)
					pDetmXrow[x] = pDetmX[(x_inv * nDetYMap) + detIndY];

				// Mapped detecor length on X
				float deltaDetmX = pDetmXrow[detIstart + detIinc] - pDetmXrow[detIstart];

				// Case 1
				// Find first detector overlap maped with pixel maped on X
				if (pDetmXrow[detIndX] - pPixmX[pixIstart] < -deltaDetmX)
					while (pDetmXrow[detIndX] - pPixmX[pixIstart] < -deltaDetmX)
						detIndX = detIndX + detIinc;

				else
					// Case 2
					// Find first pixel overlap maped with detector maped on X
					if (pDetmXrow[detIstart] - pPixmX[pixIndX] > deltaPixmX)
						while (pDetmXrow[detIstart] - pPixmX[pixIndY] > deltaPixmX)
							pixIndX = pixIndX + pixIinc;

				float moving_left_boundaryX = NULL;

				// Get the left coordinate of the first overlap on X axis
				if (pDetmXrow[detIndX] < pPixmX[pixIndX])
					moving_left_boundaryX = pPixmX[pixIndX];
				else
					moving_left_boundaryX = pDetmXrow[detIndX];


				// Loop over X intersections
				while ((detIndX < nDetX) && (pixIndX < nPixX)) {

					float gamma = (float)atan((pDetmXrow[detIndX] + (deltaDetmX / 2) - tubeX) / rtubeZ);

					// Case A, when you jump to the next detector boundarie but stay
					// in the same pixel
					if (pDetmXrow[detIndX + 1] <= pPixmX[pixIndX + 1]) {

						float overLapX = (pDetmXrow[detIndX + 1] - moving_left_boundaryX) / deltaDetmX; // Normalized overlap Calculation

						pProjt[(p * nDetY * nDetX) + (detIndX * nDetY) + detIndY] += overLapX * overLapY * pSlice[pixIndX * nPixY + pixIndY] * dz / ((float)cos(alpha) * (float)cos(gamma));

						detIndX = detIndX + detIinc;
						moving_left_boundaryX = pDetmXrow[detIndX];
					}
					else {
						// Case B, when you jump to the next pixel boundarie but stay
						// in the same detector

						float overLapX = (pPixmX[pixIndX + 1] - moving_left_boundaryX) / deltaDetmX; // Normalized overlap Calculation

						pProjt[(p * nDetY * nDetX) + (detIndX * nDetY) + detIndY] += overLapX * overLapY * pSlice[pixIndX * nPixY + pixIndY] * dz / ((float)cos(alpha) * (float)cos(gamma));

						pixIndX = pixIndX + pixIinc;
						moving_left_boundaryX = pPixmX[pixIndX];

					}

				}
				//										***** Back to Y overlap *****

				// Case A, when you jump to the next detector boundarie but stay
				// in the same pixel
				if (pDetmY[detIndY + 1] <= pPixmY[pixIndY + 1]) {
					detIndY = detIndY + detIinc;
					moving_left_boundaryY = pDetmY[detIndY];
				}
				else {
					// Case B, when you jump to the next pixel boundarie but stay
					// in the same detector
					pixIndY = pixIndY + pixIinc;
					moving_left_boundaryY = pPixmY[pixIndY];
				}

			} // Y Overlap loop

			// Free memory
			mxFree(pDetmXrow);

		} // Loop end slices

	} // Loop end Projections


	// Free memory
	mxFree(pSlice);
	mxFree(pDetX);
	mxFree(pDetY);
	mxFree(pDetZ);
	mxFree(pObjX);
	mxFree(pObjY);
	mxFree(pObjZ);
	mxFree(pDetmY);
	mxFree(pDetmX);
	mxFree(pPixmX);
	mxFree(pPixmY);
	mxFree(pRdetY);
	mxFree(pRdetZ);

	// Output memory allocation
	mwSize dims[3] = { nDetY,nDetX,nProj };
	plhs[0] = mxCreateNumericArray(3, dims, mxSINGLE_CLASS, mxREAL);
	float* pProj = (float*)mxGetData(plhs[0]);

	// Flip projection back X(Img coord is reverse to Global)
	for (int p = 0; p < nProj; p++)
		for (int x = 0, x_inv = nDetX - 1; x < nDetX; x++, x_inv--)
			for (int y = 0; y < nDetY; y++)
				pProj[(p * nDetX * nDetY) + (x * nDetY) + y] = pProjt[(p * nDetX * nDetY) + (x_inv * nDetY) + y];

	// Free memory
	mxFree(pProjt);


	return;
}

// Map on XY plane
void mapp2xy(float* const pXmapp,
	float* const pYmapp,
	float tubeX,
	float tubeY,
	float tubeZ,
	float* const pXcoord,
	float* const pYcoord,
	float* const pZcoord,
	const int nXelem,
	const int nYelem) {


	for (int x = 0; x < nXelem; x++)
		for (int y = 0; y < nYelem; y++) {

			int ind = (x * nYelem) + y;
			pXmapp[ind] = pXcoord[x] + pZcoord[y] * (pXcoord[x] - tubeX) / (tubeZ - pZcoord[y]);

			if (x == 0)
				pYmapp[y] = pYcoord[y] + pZcoord[y] * (pYcoord[y] - tubeY) / (tubeZ - pZcoord[y]);
		}


	return;
}