%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
---------------------------------------------------------------------------
                filterProj(proj,param)
---------------------------------------------------------------------------
    DESCRIPTION:
    This function filters the projection based on Feldkamp algorithms. It
    also applies a window filter to discourage noise and avoid ringing
    effect.

    The geometry is for DBT with half cone-beam. All parameters are set in 
    "ParameterSettings" code. 
 
    INPUT:

    - proj = 2D projection images 
    - param = Parameter of all geometry

    OUTPUT:

    - filteredProj = Filtered projections.

    Reference: Jiang Hsieh's book (second edition)
    Reference: Fessler Course -> (https://web.eecs.umich.edu/~fessler/course/516/l/c-tomo.pdf)

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
%% Filtering Code
function filteredProj = filterProj(proj,param)

filteredProj = zeros(size(proj),'single');

% Detector Coordinate sytem in (mm)
[uu,vv] = meshgrid(param.us,param.vs);

% Jiang Hsieh's book (second edition,page 97)
% w = (param.DSD)./sqrt((param.DSD)^2+uu.^2 + vv.^2);

% Compute weighted projections (Fessler Course Eq. (3.9.32))
weightFunction = (param.DSO.*sqrt(1+((uu./param.DSD).^2)))./...
                   (sqrt((param.DSD.^2)+(vv.^2)+(uu.^2)));

% Apply weighting function on each proj
for i=1:param.nProj
    filteredProj(:,:,i) = proj(:,:,i).* weightFunction;
end

% Increase filter length to two times nv 
h_Length = 2^nextpow2(2*param.nv);

% Builds ramp filter in space domain
ramp_kernel = ramp_builder(h_Length);

% Window filter in freq domain
H_filter = filter_window(ramp_kernel, h_Length);

% Replicate to all colluns to build a 2D filter kernel
H_filter = repmat(H_filter',1,param.nu);

% Filter each projection
for i=1:param.nProj
    
   H_proj = (zeros(h_Length,param.nu,'single'));
    
   H_proj(1:param.nv,:) = filteredProj(:,:,i);

   % Fourier transfor in projections
   H_proj = fftshift(fft(H_proj));   
    
   % Multiplication in frequency = convolution in space
   H_proj = H_proj.*H_filter;
    
   % Inverse Fourier transfor
   H_proj = (real(ifft(ifftshift(H_proj))));

   filteredProj(:,:,i) = H_proj(1:param.nv,:);

end


end

%% Function Ramp Filter
%{
The function builds Ramp Filter in space domain
Reference: Jiang Hsieh's book (second edition,page 73)
Reference: Fessler Course Eq.(3.4.14)
%}
function h = ramp_builder(h_Length)
n = (-(h_Length/2):(h_Length/2-1));
h = zeros(size(n),'single');
h(h_Length/2+1) = 1;    % Eq.(3.4.14)
odd = mod(n,2) == 1;    % Eq.(3.4.14)
h(odd) = -1 ./ (pi * n(odd)./2).^2; % Eq.(3.4.14)
end

%% Function Window Ramp Filter
%{
The function builds Ramp Filter apodizided in frequency domain
Reference: Fessler Course
%}
function H_filter = filter_window(ramp_kernel, h_Length)
H_ramp = abs(fftshift(fft(ramp_kernel)))./2;
H_window = hann(h_Length,'periodic')'; 
H_filter = H_ramp .* H_window;
end
