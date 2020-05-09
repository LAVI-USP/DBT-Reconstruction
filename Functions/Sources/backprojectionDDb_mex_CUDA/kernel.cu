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
%	  - nProj = projection number to be backprojected
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
*							CUDA Kernels									*
****************************************************************************/
__global__ void pad_projections_kernel(double* d_img, const int nDetXMap, const int nDetYMap, const int nElem, unsigned int np) {

	const int threadGId = blockIdx.x * blockDim.x + threadIdx.x;

	// Make sure we don't try and access memory outside
	// by having any threads mapped there return early
	if (threadGId >= nElem)
		return;

	d_img[(np*nDetYMap *nDetXMap) + (threadGId*nDetYMap)] = 0;

	return;
}

__global__ void map_boudaries_kernel(double* d_pBound, const int nElem, const double valueLeftBound, const double sizeElem, const double offset) {

	const int threadGId = blockIdx.x * blockDim.x + threadIdx.x;

	// Make sure we don't try and access memory outside
	// by having any threads mapped there return early
	if (threadGId >= nElem)
		return;

	d_pBound[threadGId] = (threadGId - valueLeftBound) * sizeElem + offset;

	return;
}

__global__ void rot_detector_kernel(double* d_pRdetY, double* d_pRdetZ, double* d_pYcoord, double* d_pZcoord, const double yOffset, const double zOffset,
	const double phi, const int nElem) {

	const int threadGId = blockIdx.x * blockDim.x + threadIdx.x;

	// Make sure we don't try and access memory outside
	// by having any threads mapped there return early
	if (threadGId >= nElem)
		return;

	// cos and sin are in measured in radians.

	d_pRdetY[threadGId] = ((d_pYcoord[threadGId] - yOffset)* cos(phi) - (d_pZcoord[threadGId] - zOffset)* sin(phi)) + yOffset;
	d_pRdetZ[threadGId] = ((d_pYcoord[threadGId] - yOffset)* sin(phi) + (d_pZcoord[threadGId] - zOffset)* cos(phi)) + zOffset;

	return;

}

__global__ void mapDet2Slice_kernel(double* const pXmapp, double* const pYmapp, double tubeX, double tubeY, double tubeZ, double* const pXcoord,
	double* const pYcoord, double* const pZcoord, double* const pZSlicecoord, const int nDetXMap, const int nDetYMap, const unsigned int nz) {

	/*

	Note: Matlab has linear indexing as a column of elements, i.e, the elements are actually stored in memory as queued columns.
	So threads in X are launched in the column direction (DBT Y coord), and threads in Y are launched in the row direction (DBT X coord).

	*/

	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y);

	const int thread_1D_pos = thread_2D_pos.y * nDetYMap + thread_2D_pos.x;

	// Make sure we don't try and access memory outside the detector
	// by having any threads mapped there return early
	if (thread_2D_pos.x >= nDetYMap || thread_2D_pos.y >= nDetXMap)
		return;

	pXmapp[thread_1D_pos] = ((pXcoord[thread_2D_pos.y] - tubeX)*(pZSlicecoord[nz] - pZcoord[thread_2D_pos.x]) - (pXcoord[thread_2D_pos.y] * tubeZ) + (pXcoord[thread_2D_pos.y] * pZcoord[thread_2D_pos.x])) / (-tubeZ + pZcoord[thread_2D_pos.x]);

	if (thread_2D_pos.y == 0)
		pYmapp[thread_2D_pos.x] = ((pYcoord[thread_2D_pos.x] - tubeY)*(pZSlicecoord[nz] - pZcoord[thread_2D_pos.x]) - (pYcoord[thread_2D_pos.x] * tubeZ) + (pYcoord[thread_2D_pos.x] * pZcoord[thread_2D_pos.x])) / (-tubeZ + pZcoord[thread_2D_pos.x]);

	return;
}

__global__ void img_integration_kernel(double* d_img, const int nPixX, const int nPixY, bool direction, unsigned int offsetX, unsigned int offsetY, unsigned int nSlices) {

	/*

	Integration of 2D slices over the whole volume

	(S.1.Integration. - Liu et al(2017))

	** Perfom an inclusive scan **

	*/

	const int3 memory_2D_pos = make_int3(blockIdx.x * blockDim.x + threadIdx.x + offsetX,
		blockIdx.y * blockDim.y + threadIdx.y + offsetY,
		blockIdx.z * blockDim.z + threadIdx.z);

	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y);


	// Make sure we don't try and access memory outside the detector
	// by having any threads mapped there return early
	if (memory_2D_pos.x >= nPixY || memory_2D_pos.y >= nPixX || memory_2D_pos.z >= nSlices)
		return;


	if (direction == integrateXcoord) {

		for (int s = 1; s <= blockDim.y; s *= 2) {

			int spot = thread_2D_pos.y - s;

			double val = 0;

			if (spot >= 0) {
				val = d_img[(memory_2D_pos.z*nPixY*nPixX) + (offsetY + spot) * nPixY + memory_2D_pos.x];
			}
			__syncthreads();

			if (spot >= 0) {
				d_img[(memory_2D_pos.z*nPixY*nPixX) + (memory_2D_pos.y * nPixY) + memory_2D_pos.x] += val;
			}
			__syncthreads();
		}
	}
	else
	{

		for (int s = 1; s <= blockDim.x; s *= 2) {

			int spot = thread_2D_pos.x - s;

			double val = 0;

			if (spot >= 0) {
				val = d_img[(memory_2D_pos.z*nPixY*nPixX) + memory_2D_pos.y * nPixY + spot + offsetX];
			}
			__syncthreads();

			if (spot >= 0) {
				d_img[(memory_2D_pos.z*nPixY*nPixX) + (memory_2D_pos.y * nPixY) + memory_2D_pos.x] += val;
			}
			__syncthreads();
		}

	}
	return;
}

//__global__ void bilinear_interp_kernel(double* d_interpData, cudaTextureObject_t texObj, double* const d_pDetmX, double* const d_pDetmY, const int nDetXMap, const int nDetYMap, const double2 normFactor, const double2 offset)
//{
//
//	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
//										 blockIdx.y * blockDim.y + threadIdx.y);
//
//	// Make sure we don't try and access memory outside the detector
//	// by having any threads mapped there return early
//	if (thread_2D_pos.x >= nDetYMap || thread_2D_pos.y >= nDetXMap)
//		return;
//
//	// Normalize coordinates from 0-1 (coords on texture memory)
//	double xQuery = d_pDetmX[thread_2D_pos.y * nDetYMap + thread_2D_pos.x] / normFactor.x + offset.x + (-.5f * (1 / normFactor.x));
//	double yQuery = d_pDetmY[thread_2D_pos.x] / normFactor.y + offset.y + (.5f * (1 / normFactor.y));
//
//	// Read from texture and write to global memory
//	d_interpData[thread_2D_pos.y * nDetYMap + thread_2D_pos.x] = tex2D<double>(texObj, xQuery, yQuery);
//
//	return;
//}

__global__ void bilinear_interpolation_kernel_GPU(double* d_sliceI, double* d_pProj, double* d_pObjX, double* d_pObjY, double* d_pDetmX, double* d_pDetmY, const int nPixXMap, const int nPixYMap,
	const int nDetXMap, const int nDetYMap, const int nDetX, const int nDetY, const unsigned int np) {

	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
										 blockIdx.y * blockDim.y + threadIdx.y);

	// Make sure we don't try and access memory outside the detector
	// by having any threads mapped there return early
	if (thread_2D_pos.x >= nPixYMap || thread_2D_pos.y >= nPixXMap)
		return;

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

	----> Xcoord (thread_2D_pos.y)
	|							0___________
	v							|_d00_|_d10_|
	Ycoord (thread_2D_pos.x)	|_d01_|_d11_|

	Matlab code --
	objX = nDetX - (objX ./ detmX(1, end - 1));
	objY = (objY ./ detmX(1, end - 1)) - (detmY(1) / detmX(1, end - 1));

	*/

	// Adjust the mapped coordinates to cross the range of (0-nDetX).*duMap 
	//  Divide by pixelSize to get a unitary pixel size
	const double xNormData = nDetX - d_pObjX[thread_2D_pos.y] / d_pDetmX[0];
	const signed int    xData = floor(xNormData);
	const double  alpha = xNormData - xData;

	// Adjust the mapped coordinates to cross the range of (0-nDetY).*dyMap  
	//  Divide by pixelSize to get a unitary pixel size
	const double yNormData = (d_pObjY[thread_2D_pos.x] / d_pDetmX[0]) - (d_pDetmY[0] / d_pDetmX[0]);
	const signed int    yData = floor(yNormData);
	const double  beta = yNormData - yData;

	double d00, d01, d10, d11;
	if (((xNormData) >= 0) && ((xNormData) <= nDetX) && ((yNormData) >= 0) && ((yNormData) <= nDetY))	d00 = d_pProj[(np*nDetYMap*nDetXMap) + (xData*nDetYMap + yData)];						else    d00 = 0.0;
	if (((xData + 1) > 0) && ((xData + 1) <= nDetX) && ((yNormData) >= 0) && ((yNormData) <= nDetY))	d10 = d_pProj[(np*nDetYMap*nDetXMap) + ((xData + 1)*nDetYMap + yData)]; 				else    d10 = 0.0;
	if (((xNormData) >= 0) && ((xNormData) <= nDetX) && ((yData + 1) > 0) && ((yData + 1) <= nDetY))	d01 = d_pProj[(np*nDetYMap*nDetXMap) + (xData*nDetYMap + yData + 1)]; 				else    d01 = 0.0;
	if (((xData + 1) > 0) && ((xData + 1) <= nDetX) && ((yData + 1) > 0) && ((yData + 1) <= nDetY))		d11 = d_pProj[(np*nDetYMap*nDetXMap) + ((xData + 1)*nDetYMap + yData + 1)];			else    d11 = 0.0;

	double result_temp1 = alpha * d10 + (-d00 * alpha + d00);
	double result_temp2 = alpha * d11 + (-d01 * alpha + d01);

	d_sliceI[thread_2D_pos.y * nPixYMap + thread_2D_pos.x] = beta * result_temp2 + (-result_temp1 * beta + result_temp1);

}

__global__ void differentiation_kernel(double* d_pVolume, double* d_sliceI, double tubeX, double rtubeY, double rtubeZ,
	double* const d_pObjX, double* const d_pObjY, double* const d_pObjZ, const int nPixX, const int nPixY, const int nPixXMap,
	const int nPixYMap, const double du, const double dv, const double dx, const double dy, const double dz, const unsigned int nz) {

	const int2 thread_2D_pos = make_int2(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y);

	const int thread_1D_pos = (nPixX*nPixY*nz) + (thread_2D_pos.y * nPixY) + thread_2D_pos.x;

	// Make sure we don't try and access memory outside the detector
	// by having any threads mapped there return early
	if (thread_2D_pos.x >= nPixY || thread_2D_pos.y >= nPixX)
		return;

	/*

	S.3. Differentiation - Eq. 24 - Liu et al (2017)

	Detector integral projection
	___________
	|_A_|_B_|___|
	|_C_|_D_|___|
	|___|___|___|


	(thread_2D_pos.x,thread_2D_pos.y)
	________________
	|_A_|__B__|_____|
	|_C_|(0,0)|(0,1)|
	|___|(1,0)|(1,1)|

	Threads are lauched from D up to nPixX (thread_2D_pos.y) and nPixY (thread_2D_pos.x)
	i.e., they are running on the detector image. Thread (0,0) is on D.

	Coordinates on intergal projection:

	A = thread_2D_pos.y * nPixYMap + thread_2D_pos.x
	B = ((thread_2D_pos.y+1) * nPixYMap) + thread_2D_pos.x
	C = thread_2D_pos.y * nPixYMap + thread_2D_pos.x + 1
	D = ((thread_2D_pos.y+1) * nPixYMap) + thread_2D_pos.x + 1

	*/

	unsigned int coordA = thread_2D_pos.y * nPixYMap + thread_2D_pos.x;
	unsigned int coordB = ((thread_2D_pos.y + 1) * nPixYMap) + thread_2D_pos.x;
	unsigned int coordC = coordA + 1;
	unsigned int coordD = coordB + 1;

	// x - ray angle in X coord
	double gamma = atan((d_pObjX[thread_2D_pos.y] + (dx / 2.0) - tubeX) / (rtubeZ - d_pObjZ[nz]));

	// x - ray angle in Y coord
	double alpha = atan((d_pObjY[thread_2D_pos.x] + (dy / 2.0) - rtubeY) / (rtubeZ - d_pObjZ[nz]));


	double dA, dB, dC, dD;

	dA = d_sliceI[coordA];
	dB = d_sliceI[coordB];
	dC = d_sliceI[coordC];
	dD = d_sliceI[coordD];

	// Treat border of interpolated integral detector
	if (dC == 0 && dD == 0) {
		dC = dA;
		dD = dB;
	}


	// S.3.Differentiation - Eq. 24 - Liu et al(2017)
	d_pVolume[thread_1D_pos] += ((dD - dC - dB + dA)*(du*dv*dz / (cos(alpha)*cos(gamma)*dx*dy)));

	return;
}

__global__ void division_kernel(double* d_img, const int nPixX, const int nPixY, unsigned int nSlices, unsigned int nProj){

	const int3 thread_2D_pos = make_int3(blockIdx.x * blockDim.x + threadIdx.x,
		blockIdx.y * blockDim.y + threadIdx.y,
		blockIdx.z * blockDim.z + threadIdx.z);

	const int thread_1D_pos = (nPixX*nPixY*thread_2D_pos.z) + (thread_2D_pos.y * nPixY) + thread_2D_pos.x;

	// Make sure we don't try and access memory outside the detector
	// by having any threads mapped there return early
	if (thread_2D_pos.x >= nPixY || thread_2D_pos.y >= nPixX || thread_2D_pos.z >= nSlices)
		return;

	d_img[thread_1D_pos] /= (double) nProj;

}


/****************************************************************************
* function: backprojectionDDb() - CUDA backprojection Branchless Distance Driven.	*
****************************************************************************/

void backprojectionDDb(double* const h_pVolume,
	double* const h_pProj,
	double* const h_pTubeAngle,
	double* const h_pDetAngle,
	const double idXProj,
	const unsigned int nProj,
	const unsigned int nPixX,
	const unsigned int nPixY,
	const unsigned int nSlices,
	const unsigned int nDetX,
	const unsigned int nDetY,
	const double dx,
	const double dy,
	const double dz,
	const double du,
	const double dv,
	const double DSD,
	const double DDR,
	const double DAG)
{

	// Unique GPU
	int devID = 0;

	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, devID);
	cudaCheckErrors("cudaGetDeviceProperties");

	const unsigned int maxThreadsPerBlock = deviceProp.maxThreadsPerBlock;


	mexPrintf("GPU Device %d: \"%s\" with compute capability %d.%d has %d Multi-Processors and %zu bytes of global memory\n\n", devID,
		deviceProp.name, deviceProp.major, deviceProp.minor, deviceProp.multiProcessorCount, deviceProp.totalGlobalMem);

	
	// Create timmer variables and start it
	StopWatchInterface *timer = NULL;
	sdkCreateTimer(&timer);
	sdkStartTimer(&timer);
	
	//cudaStream_t stream1;
	//cudaStreamCreate(&stream1);

	dim3 threadsPerBlock(1, 1, 1);
	dim3 blockSize(1, 1, 1);


	// Number of mapped detector and pixels
	const int nDetXMap = nDetX + 1;
	const int nDetYMap = nDetY + 1;
	const int nPixXMap = nPixX + 1;
	const int nPixYMap = nPixY + 1;


	// Convention: h_ variables live on host
	// Convention: d_ variables live on device (GPU global mem)
	double* d_pProj;
	double* d_sliceI;
	double* d_pVolume;
	double* d_pTubeAngle;
	double* d_pDetAngle;
	
	
	// Allocate global memory on the device, place result in "d_----"
	cudaMalloc((void **)&d_pProj, nDetXMap*nDetYMap*nProj * sizeof(double)); 
	cudaMalloc((void **)&d_sliceI, nPixXMap*nPixYMap * sizeof(double));
	cudaMalloc((void **)&d_pVolume, nPixX*nPixY*nSlices * sizeof(double));
	cudaMalloc((void **)&d_pTubeAngle, nProj * sizeof(double));
	cudaMalloc((void **)&d_pDetAngle, nProj * sizeof(double));
	
	cudaCheckErrors("cudaMalloc Initial");
	
	// Copy data from host memory "h_----" to device memory "d_----"
	cudaMemcpy((void *)d_pTubeAngle, (void *)h_pTubeAngle, nProj * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy((void *)d_pDetAngle, (void *)h_pDetAngle, nProj * sizeof(double), cudaMemcpyHostToDevice);

	cudaCheckErrors("cudaMemcpy Angles");
	

	/*
	Copy projection data from host memory "h_----" to device memory "d_----" padding with zeros for image integation
	*/
	

	// Initialize first column and row with zeros
	double* h_pProj_tmp;
	double* d_pProj_tmp;

	threadsPerBlock.x = maxThreadsPerBlock;
	blockSize.x = (nDetXMap / maxThreadsPerBlock) + 1;

	for (unsigned int np = 0; np < nProj; np++) {

		// Pad on X coord direction
		pad_projections_kernel << <blockSize, threadsPerBlock >> > (d_pProj, nDetXMap, nDetYMap, nDetXMap, np);

		// Pad on Y coord direction
		d_pProj_tmp = d_pProj + (nDetXMap*nDetYMap*np) + 1;
		cudaMemset(d_pProj_tmp, 0, nPixY * sizeof(double));
		cudaCheckErrors("cudaMemset Padding Projections");
	}
	
	
	// Copy projections data from host memory
	for (unsigned int np = 0; np < nProj; np++)
		for (unsigned int c = 0; c < nDetX; c++) {

			h_pProj_tmp = h_pProj + (c *nDetY) + (nDetX*nDetY*np);
			d_pProj_tmp = d_pProj + (((c + 1) *nDetYMap) + 1) + (nDetXMap*nDetYMap*np);

			cudaMemcpy((void *)d_pProj_tmp, (void *)h_pProj_tmp, nDetY * sizeof(double), cudaMemcpyHostToDevice);
			cudaCheckErrors("cudaMemcpy Projections");
		}

	
	
	// Pointer for projections coordinates
	double* d_pDetX;
	double* d_pDetY;
	double* d_pDetZ;
	double* d_pObjX;
	double* d_pObjY;
	double* d_pObjZ;

	// Allocate global memory on the device for projections coordinates
	cudaMalloc((void **)&d_pDetX, nDetXMap * sizeof(double));
	cudaMalloc((void **)&d_pDetY, nDetYMap * sizeof(double));
	cudaMalloc((void **)&d_pDetZ, nDetYMap * sizeof(double));
	cudaMalloc((void **)&d_pObjX, nPixXMap * sizeof(double));
	cudaMalloc((void **)&d_pObjY, nPixYMap * sizeof(double));
	cudaMalloc((void **)&d_pObjZ, nSlices * sizeof(double));

	cudaCheckErrors("cudaMalloc Coordinates");



	// Pointer for mapped coordinates
	double* d_pDetmY;
	double* d_pDetmX;


	// Allocate global memory on the device for mapped coordinates
	cudaMalloc((void **)&d_pDetmY, nDetYMap * sizeof(double));
	cudaMalloc((void **)&d_pDetmX, nDetYMap * nDetXMap * sizeof(double));

	cudaCheckErrors("cudaMalloc Map-Coordinates");


	// Pointer for rotated detector coords
	double* d_pRdetY;
	double* d_pRdetZ;

	// Allocate global memory on the device for for rotated detector coords
	cudaMalloc((void **)&d_pRdetY, nDetYMap * sizeof(double));
	cudaMalloc((void **)&d_pRdetZ, nDetYMap * sizeof(double));

	cudaCheckErrors("cudaMalloc Rot-Coordinates");
	
	// Generate detector and object boudaries
	
	threadsPerBlock.x = maxThreadsPerBlock;

	blockSize.x = (nDetX / maxThreadsPerBlock) + 1;

	map_boudaries_kernel << <blockSize, threadsPerBlock >> > (d_pDetX, nDetXMap, (double)nDetX, -du, 0.0);

	blockSize.x = (nDetY / maxThreadsPerBlock) + 1;

	map_boudaries_kernel << <blockSize, threadsPerBlock >> > (d_pDetY, nDetYMap, nDetY / 2.0, dv, 0.0);

	blockSize.x = (nPixX / maxThreadsPerBlock) + 1;

	map_boudaries_kernel << <blockSize, threadsPerBlock >> > (d_pObjX, nPixXMap, (double)nPixX, -dx, 0.0);

	blockSize.x = (nPixY / maxThreadsPerBlock) + 1;

	map_boudaries_kernel << <blockSize, threadsPerBlock >> > (d_pObjY, nPixYMap, nPixY / 2.0, dy, 0.0);

	blockSize.x = (nSlices / maxThreadsPerBlock) + 1;

	map_boudaries_kernel << <blockSize, threadsPerBlock >> > (d_pObjZ, nSlices, 0.0, dz, DAG + (dz / 2.0));

	
	//mexPrintf("Map Boundaries -- Threads:%d Blocks:%d \n", threads, blocks);

	
	// Initiate variables value with 0
	cudaMemset(d_pDetZ, 0, nDetYMap * sizeof(double));
	cudaMemset(d_pVolume, 0, nPixX * nPixY * nSlices * sizeof(double));

	cudaCheckErrors("cudaMemset Zeros");
	
	//cudaDeviceSynchronize();

	// X - ray tube initial position
	double tubeX = 0;
	double tubeY = 0;
	double tubeZ = DSD;

	// Iso - center position
	double isoY = 0;
	double isoZ = DDR;
	

	// Integration of 2D projection over the whole projections
	// (S.1.Integration. - Liu et al(2017))

	// Naive integration o the X coord
	threadsPerBlock.x = 10;
	threadsPerBlock.y = 10;
	threadsPerBlock.z = 10;

	blockSize.x = (unsigned int)ceil(((double)nDetYMap / (threadsPerBlock.x - 1)));
	blockSize.y = 1;
	blockSize.z = (unsigned int)ceil((double)nProj / threadsPerBlock.z);

	//mexPrintf("Divisao: %f \n", (double)nDetY / (threadsPerBlock.x - 1));
	//mexPrintf("Divisao ceil: %d \n", (unsigned int)ceil((double)nDetY / (threadsPerBlock.x - 1)));

	for (int k = 0; k < ceil((double)nDetXMap / (threadsPerBlock.x - 1)); k++) {

		img_integration_kernel << <blockSize, threadsPerBlock >> > (d_pProj, nDetXMap, nDetYMap, integrateXcoord, 0, k * 9, nProj);

		//mexPrintf("integration kernel -- threadsX:%d blocksX:%d threadsY:%d blocksY:%d \n", threadsPerBlock.x, blockSize.x, threadsPerBlock.y, blockSize.y);

	}

	// Naive integration o the Y coord
	threadsPerBlock.x = 10;
	threadsPerBlock.y = 10;
	threadsPerBlock.z = 10;

	blockSize.x = 1;
	blockSize.y = (unsigned int)ceil((double)nDetXMap / (threadsPerBlock.y - 1));
	blockSize.z = (unsigned int)ceil((double)nProj / threadsPerBlock.z);


	for (int k = 0; k < ceil((double)nDetYMap / (threadsPerBlock.y - 1)); k++) {

		img_integration_kernel << <blockSize, threadsPerBlock >> > (d_pProj, nDetXMap, nDetYMap, integrateYcoord, k * 9, 0, nProj);

		//mexPrintf("integration kernel -- threadsX:%d blocksX:%d threadsY:%d blocksY:%d \n", threadsPerBlock.x, blockSize.x, threadsPerBlock.y, blockSize.y);

	}
	
	double* d_pDetmX_tmp = d_pDetmX + (nDetYMap * (nDetXMap-2));

	// Test if we will loop over all projs or not
	unsigned int projIni, projEnd, nProj2Run;
	if(idXProj == -1){
		projIni = 0;
		projEnd = nProj;
		nProj2Run = nProj;
	}
	else{
		nProj2Run = 1;
		projIni = (unsigned int) idXProj;
		projEnd = (unsigned int) idXProj + 1;
	}
	
	
	// For each projection
	for (unsigned int p = projIni; p < projEnd; p++) {
	
		// Get specif tube angle for the projection
		double theta = h_pTubeAngle[p] * M_PI / 180.0;

		// Get specif detector angle for the projection
		double phi = h_pDetAngle[p] * M_PI / 180.0;

		//mexPrintf("Tube angle:%f Det angle:%f\n", theta, phi);

		// Tube rotation
		double rtubeY = ((tubeY - isoY)*cos(theta) - (tubeZ - isoZ)*sin(theta)) + isoY;
		double rtubeZ = ((tubeY - isoY)*sin(theta) + (tubeZ - isoZ)*cos(theta)) + isoZ;

		//mexPrintf("R tube Y:%f R tube Z:%f\n", rtubeY, rtubeZ);

		// Detector rotation
		threadsPerBlock.x = maxThreadsPerBlock;
		threadsPerBlock.y = 1;
		threadsPerBlock.z = 1;

		blockSize.x = (nDetYMap / maxThreadsPerBlock) + 1;
		blockSize.y = 1;
		blockSize.z = 1;

		rot_detector_kernel << <blockSize, threadsPerBlock >> > (d_pRdetY, d_pRdetZ, d_pDetY, d_pDetZ, isoY, isoZ, phi, nDetYMap);

		//cudaDeviceSynchronize();

		//mexPrintf("Detector rotation -- Threads:%d Blocks:%d \n", threads, blocks);


		// For each slice
		for (unsigned int nz = 0; nz < nSlices; nz++) {
			
			/*

			Map detector onto XY plane(Inside proj loop in case detector rotates)

			*** Note: Matlab has linear indexing as a column of elements, i.e, the elements are actually stored in memory as queued columns.
			So threads in X are launched in the column direction (DBT Y coord), and threads in Y are launched in the row direction (DBT X coord).

			*/
		
			threadsPerBlock.x = 32;
			threadsPerBlock.y = 32;
			threadsPerBlock.z = 1;

			blockSize.x = (nDetYMap / threadsPerBlock.x) + 1;
			blockSize.y = (nDetXMap / threadsPerBlock.y) + 1;
			blockSize.z = 1;


			mapDet2Slice_kernel << <blockSize, threadsPerBlock >> > (d_pDetmX, d_pDetmY, tubeX, rtubeY, rtubeZ, d_pDetX, d_pRdetY, d_pRdetZ, d_pObjZ, nDetXMap, nDetYMap, nz);
			//cudaDeviceSynchronize();

			//mexPrintf("Map detector onto XY plane -- ThreadsX:%d ThreadsY:%d BlocksX:%d BlocksY:%d\n", threadsPerBlock.x, threadsPerBlock.y, blockSizeX, blockSizeY);

			/*
			Creates texture memory
			*/

			/*

			// Allocate CUDA array in device memory
			cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0,cudaChannelFormatKindFloat);

			cudaArray* d_imageArray;
			cudaMallocArray(&d_imageArray, &channelDesc, nPixXMap, nPixYMap);

			// Copy to texture memory
			checkCudaErrors(cudaMemcpyToArray(d_imageArray, 0, 0, d_pSlice, nPixXMap* nPixYMap * sizeof(double), cudaMemcpyDeviceToDevice));

			// Specify texture
			struct cudaResourceDesc resDesc;
			memset(&resDesc, 0, sizeof(resDesc));
			resDesc.resType = cudaResourceTypeArray;
			resDesc.res.array.array = d_imageArray;

			// Specify texture object parameters
			struct cudaTextureDesc texDesc;
			memset(&texDesc, 0, sizeof(texDesc));
			texDesc.addressMode[0] = cudaAddressModeBorder;
			texDesc.addressMode[1] = cudaAddressModeBorder;
			texDesc.filterMode = cudaFilterModeLinear;
			texDesc.readMode = cudaReadModeElementType;
			texDesc.normalizedCoords = true;

			// Create texture object
			cudaTextureObject_t texObj = 0;
			cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL);

			// Allocate result of transformation in device memory
			double* d_sliceI;
			cudaMalloc(&d_sliceI, nDetXMap*nDetYMap * sizeof(double));

			// Invoke kernel
			blockSize.x = (nDetYMap / threadsPerBlock.x) + 1;
			blockSize.y = (nDetXMap / threadsPerBlock.y) + 1;

			const double2 normFactor = make_double2(nPixX*-dx,
			nPixY*dy);

			const double2 offset = make_double2(1,
			0.5);

			bilinear_interp_kernel << <blockSize, threadsPerBlock >> > (d_sliceI, texObj, d_pDetmX, d_pDetmY, nDetXMap, nDetYMap, normFactor, offset);

			mexPrintf("Texture interpolation -- Threads:%d Blocks:%d \n", threads, blocks);

			Destroy texture object
			cudaDestroyTextureObject(texObj);
			*/


			/*
			S.2. Interpolation - Liu et al (2017)
			*/

			blockSize.x = (nPixYMap / threadsPerBlock.x) + 1;
			blockSize.y = (nPixXMap / threadsPerBlock.y) + 1;
			
			bilinear_interpolation_kernel_GPU << <blockSize, threadsPerBlock >> > (d_sliceI, d_pProj, d_pObjX, d_pObjY, d_pDetmX_tmp, d_pDetmY, nPixXMap, nPixYMap, nDetXMap, nDetYMap, nDetX, nDetY, p);


			/*
			S.3. Differentiation - Eq. 24 - Liu et al (2017)
			*/

			blockSize.x = (nPixY / threadsPerBlock.x) + 1;
			blockSize.y = (nPixX / threadsPerBlock.y) + 1;

			differentiation_kernel << <blockSize, threadsPerBlock >> > (d_pVolume, d_sliceI, tubeX, rtubeY, rtubeZ, d_pObjX, d_pObjY, d_pObjZ, nPixX, nPixY, nPixXMap, nPixYMap, du, dv, dx, dy, dz, nz);

		} // Loop end slices

	} // Loop end Projections


	// Normalize volume dividing by the number of projs
	threadsPerBlock.x = 10;
	threadsPerBlock.y = 10;
	threadsPerBlock.z = 10;

	blockSize.x = (nPixY / threadsPerBlock.x) + 1;
	blockSize.y = (nPixX / threadsPerBlock.y) + 1;
	blockSize.z = (nSlices / threadsPerBlock.z) + 1;

	division_kernel << <blockSize, threadsPerBlock >> > (d_pVolume, nPixX, nPixY, nSlices, nProj2Run);

	//double* const h_var1 = (double*)mxMalloc(1 * sizeof(double));
	//cudaMemcpy((void *)h_var1, (void *)d_pDetmX_tmp, 1 * sizeof(double), cudaMemcpyDeviceToHost);
	//mexPrintf("P: %d Nz: %d Var:%f \n",p,nz,h_var1[0]);

	//nDetYMap
	//double* const h_pPixmX = (double*)mxMalloc(nDetXMap*nDetYMap * sizeof(double));
	//double* const h_pPixmY = (double*)mxMalloc(nDetYMap * sizeof(double));
	//cudaMemcpy((void *)h_pPixmX, (void *)d_pDetmX, nDetXMap*nDetYMap * sizeof(double), cudaMemcpyDeviceToHost);
	//cudaMemcpy((void *)h_pPixmY, (void *)d_pDetmY, nDetYMap * sizeof(double), cudaMemcpyDeviceToHost);

	//for (int u = 0; u < nDetXMap; u++)
	//	for (int v = 0; v < nDetYMap; v++) {
	//		h_pProj[(0 * nDetYMap*nDetXMap) + (u*nDetYMap) + v] = h_pPixmX[(u*nDetYMap) + v];
	//		h_pProj[(1 * nDetYMap*nDetXMap) + (u*nDetYMap) + v] = h_pPixmY[v];
	//	}

	// d_sliceI
	//double* const h_var = (double*)mxMalloc(nPixXMap*nPixYMap * sizeof(double));
	//cudaMemcpy((void *)h_var, (void *)d_sliceI, nPixXMap*nPixYMap * sizeof(double), cudaMemcpyDeviceToHost);

	//for (int u = 0; u < nPixXMap; u++)
	//	for (int v = 0; v < nPixYMap; v++) {
	//		h_pVolume[(0 * nPixYMap*nPixXMap) + (u*nPixYMap) + v] = h_var[(u*nPixYMap) + v];
	//	}
	

	// d_pProj
	//cudaMemcpy((void *)h_pVolume, (void *)d_pProj, nProj*nDetXMap*nDetYMap * sizeof(double), cudaMemcpyDeviceToHost);

	// d_pVolume
	cudaMemcpy((void *)h_pVolume, (void *)d_pVolume, nSlices* nPixX * nPixY * sizeof(double), cudaMemcpyDeviceToHost);
	cudaCheckErrors("cudaMemcpy Final");
	

	cudaFree(d_pProj);
	cudaFree(d_sliceI);
	cudaFree(d_pVolume);
	cudaFree(d_pTubeAngle);
	cudaFree(d_pDetAngle);
	cudaFree(d_pDetX);
	cudaFree(d_pDetY);
	cudaFree(d_pDetZ);
	cudaFree(d_pObjX);
	cudaFree(d_pObjY);
	cudaFree(d_pObjZ);
	cudaFree(d_pDetmY);
	cudaFree(d_pDetmX);
	cudaFree(d_pRdetY);
	cudaFree(d_pRdetZ);
	
	cudaCheckErrors("cudaFree Final");
	
	sdkStopTimer(&timer);
	mexPrintf("Processing time: %f (ms)\n", sdkGetTimerValue(&timer));
	sdkDeleteTimer(&timer);

	cudaDeviceReset();

	return;

}