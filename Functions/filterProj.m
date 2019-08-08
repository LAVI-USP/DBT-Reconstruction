%% Author: Rodrigo de Barros Vimieiro
% Date: April, 2018
% rodrigo.vimieiro@gmail.com
% =========================================================================
%{
% -------------------------------------------------------------------------
%                 filterProj(proj,param,cutoff)
% -------------------------------------------------------------------------
%     DESCRIPTION:
%     This function filters the projection based on Feldkamp algorithms. It
%     also applies a window filter to discourage noise and avoid ringing
%     effect.
% 
%     The geometry is for DBT with half cone-beam. All parameters are set 
%      in "ParameterSettings" code. 
%  
%     INPUT:
% 
%     - proj = 2D projection images 
%     - param = Parameter of all geometry
%     - cutoff = Cut off frequency up to Nyquist
% 
%     OUTPUT:
% 
%     - filteredProj = Filtered projections.
% 
%     Reference: Jiang Hsieh's book (second edition)
%     Reference: Fessler, J. A.: Fundamentals of CT Reconstruction in 2D 
%     and 3D. In: Brahme, A. (eds.) Comprehensive Biomedical Physics, 
%     1st ed., vol. 2, pp 263-295. Elsevier, Netherlands (2014).
%     doi:10.1016/B978-0-444-53632-7.00212-4.
%     Reference: Fessler Book -> (http://web.eecs.umich.edu/~fessler/book/c-tomo.pdf)
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
%% Filtering Code
function filteredProj = filterProj(proj,param,cutoff)

filteredProj = zeros(size(proj),'single');

% Detector Coordinate sytem in (mm)
[uCoord,vCoord] = meshgrid(param.us,param.vs);

% Compute weighted projections (Fessler Book Eq. (3.10.6))
weightFunction = param.DSO ./ sqrt((param.DSD.^2)+(vCoord.^2)+(uCoord.^2));
                   
% Apply weighting function on each proj
for i=1:param.nProj
    filteredProj(:,:,i) = proj(:,:,i).* weightFunction;
end

% Increase filter length to two times nv to decrease freq step
h_Length = 2^nextpow2(2*param.nv);

% Builds ramp filter in space domain
ramp_kernel = ramp_builder(h_Length);

% Window filter in freq domain
H_filter = filter_window(ramp_kernel, h_Length, cutoff);

% Replicate to all colluns to build a 2D filter kernel
H_filter = repmat(H_filter',1,param.nu);

% Proj in freq domain
H_proj = (zeros(h_Length,param.nu,'single'));

% Filter each projection
for i=1:param.nProj
     
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
Reference: Jiang Hsieh's book (second edition,page 73) Eq. 3.29
Reference: Fessler Book Eq.(3.4.14)
%}
function h = ramp_builder(h_Length)

n = (-(h_Length/2):(h_Length/2-1));
h = zeros(size(n),'single');
h(h_Length/2+1) = 1/4;  % Eq. 3.29
odd = mod(n,2) == 1;    % Eq. 3.29
h(odd) = -1 ./ (pi * n(odd)).^2; % Eq. 3.29
end

%% Function Window Ramp Filter
%{
The function builds Ramp Filter apodizided in frequency domain
Reference: Fessler Book and MIRT
%}
function H_filter = filter_window(ramp_kernel, h_Length, cutoff)

H_ramp = abs(fftshift(fft(ramp_kernel))); % Bring filter to freq domain

w = round(h_Length*cutoff); % Cut off freq
n = (-(h_Length/2):(h_Length/2-1));
H_window = 0.5 * (1 + cos(2*pi*n/w)); % Hanning filter
H_window = H_window .* (abs(n) <w/2); % Apply cut off freq

H_filter = H_ramp .* H_window; % Apply window filter
end
