%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% 
%     Reference: 
%     https://www.codefull.net/2014/11/tutorial-use-cuda-code-in-matlab/
% 
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
%% Compilation Code

% Linux systems
if isunix
    
    fprintf('--- Starting complilation ---\n');
    
    matlabRootPath = matlabroot;
     
    % ------------------------------------------
    
    pathVar{1} = 'backprojectionDD_mex';
    pathVar{2} = 'projectionDD_mex';
    
    for k=1:2
        mex('-largeArrayDims',[pathVar{k} filesep pathVar{k} '.cpp'],['-L' matlabRootPath '/bin/glnxa64'],'-outdir','..');
    end
    
    % ------------------------------------------
    
    pathVar{1} = 'backprojectionDDb_mex';
    pathVar{2} = 'projectionDDb_mex';
    
    for k=1:2
        mex('-largeArrayDims',[pathVar{k} filesep pathVar{k} '.cpp'],['-L' matlabRootPath '/bin/glnxa64'],'CFLAGS="$CFLAGS -fopenmp"', 'LDFLAGS="$LDFLAGS -fopenmp"','-outdir','..');
    end
      
    % ------------------------------------------

    % https://www.mathworks.com/matlabcentral/answers/394830-nvcc-command-not-found
    cudaPath = '/usr/local/cuda';
    
    if(~system('cd /usr/local/cuda'))
        p = getenv('PATH');
        setenv('PATH', [p ':' cudaPath '/bin'])
        
        pathVar{1} = 'backprojectionDDb_mex_CUDA';  sourceFilePath{1} = 'backprojectionDDb_mex_cuda.cpp';
        pathVar{2} = 'projectionDDb_mex_CUDA';      sourceFilePath{2} = 'projectionDDb_mex_cuda.cpp';


        for k=1:2
            % https://stackoverflow.com/a/25197234/8682939
            system(['nvcc -O3 -c ' pathVar{k} '/kernel.cu -Xcompiler -fPIC -I' matlabRootPath '/extern/include -I' matlabRootPath '/toolbox/distcomp/gpu/extern/include -I' cudaPath '/samples/common/inc']);
            mex('-R2018a',[pathVar{k} filesep sourceFilePath{k}],'kernel.o',['-I' cudaPath '/include'],['-I' cudaPath '/samples/common/inc'],['-L' cudaPath '/lib64'],['-L' matlabRootPath '/bin/glnxa64'],'-lcudart', '-lcufft', '-lmwgpu','-outdir','..');    
        end

        system('rm kernel.o');
    
        flagCUDA = 1;
    else
        flagCUDA = 0;
    end
    
    % ------------------------------------------
    
    fprintf('\n\n\n\nC++ files compiled.\n');
    fprintf('C++ (OpenMP) files compiled.\n');
    if(flagCUDA)
        fprintf('CUDA files compiled.\n');
    else
        error('Either change CUDA toolkit path or install it.');
    end
    
    
    fprintf('\n--- Done!! --- \n');

end




