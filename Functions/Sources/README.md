# Compiling mex/mex-CUDA files

If you want to compile mex-CUDA files, make sure to follow the **(optional)** steps.

## What is necessary?

 1. Download [Visual Studio (VS) IDE](https://visualstudio.microsoft.com/vs/).

 2. A CUDA-capable GPU **(optional)**.
 
 2. Download [CUDA Toolkit](https://developer.nvidia.com/cuda-toolkit) **(optional)**.

## How to install?

 1. Install VS with C++ components.

 2. Install CUDA Toolkit with VS integration (default installation) **(optional)**.

## How to compile the codes?

 1. Open the project/solution on each folder.
 
 2. Change the paths according to your directories, for example:

     * Configuration properties -> VC++ Directories -> Include Directories:
        > C:\ProgramData\NVIDIA Corporation\CUDA Samples\v10.1\common\inc **(optional)**

        > C:\Program Files\MATLAB\R2017b\extern\include

     * Configuration properties -> Linker -> General:
        > C:\Program Files\MATLAB\R2017b\extern\lib\win64\microsoft

## How to test?

 1. Add the release path on MATLAB, for example: 

     * `addpath('C:\DBT-Reconstruction\Functions\Sources\projectionDDb_mex_CUDA\x64\Release')`
     * `addpath('C:\DBT-Reconstruction\Functions\Sources\backprojectionDDb_mex_CUDA\x64\Release')`

 2. Or place the **.mexw64** files on the _Functions_ folder.

 3. Run the main _Reconstruction.m_ file to get the parameters and to load projections.

 4. Run the following command to get a simple Back-projection:
     * `dataRecon3d_mex_CUDA = backprojectionDDb_mex_CUDA(double(dataProj),parameter);`

 5. You should get something like this:

     > GPU Device 0: "GeForce GTX 1050 Ti" with compute capability 6.1 has 6 Multi-Processors and 4294967296 bytes of global memory


**Note:** All these codes were developed and tested on:

| Software      | Version     
| ------------- |:-------------:
| MATLAB            | 2017b 
| Visual Studio IDE | 2019      
| CUDA Toolkit      | 10.1  
| GPU               | GeForce GTX 1050 Ti  
| Windows           | 10     

## Having problems?

Please report issues [here](https://github.com/LAVI-USP/DBT-Reconstruction/issues). I also can send you some already compiled mex files. 

---

Laboratory of Computer Vision ([Lavi](http://iris.sel.eesc.usp.br/lavi/))  
Department of Electrical and Computer Engineering  
São Carlos School of Engineering, University of São Paulo  
São Carlos - Brazil  