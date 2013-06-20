function f = chimaps_kernel(datasize, B0, B0dir)
% KERNEL = chimaps_kernel(DATASIZE, B0, B0DIR)
%
% Create spherical filter kernel
%
% DATASIZE      size of data array
% B0            strength of B0 field (in Tesla)
% B0DIR         direction of B0 field (see below, default is 3rd array 
%               dimension (z))
%
% Implements the kernel, F, in the equation 
%   
%       Delta(B) =  F x (chi_sphere - chi_external)
%
% where
%
%             B0    3(ZDIR.R)R - ZDIR   4 pi a^3
%       F =  ---- x ----------------- x --------
%            4 pi         |R|^3            3
%
% where chi_sphere is the susecptibility of the sphere
%       chi_external is the susceptibility outside the sphere
%       a is the radius of the sphere
%
% The second term is the dipole tensor and the third term is the volume of
% the sphere.
%
% FILTERSIZE & DUTY COVERAGE
% ==========================
%   Note that reduced duty coverage has a significant effect on the
%   accuracy of the kernel. A filter that covers 99% of the duty has a
%   very noticeable error in the resultant estimation of field shift. It is
%   therefore recommended to specify a duty coverage of 1. This will set
%   the filter size to double the datasize, thus ensuring the best duty
%   coverage possible.
%
% B0DIR
% =====
%   Indicates the direction of the B0 field. It can be set to
%       1: first array dimension
%       2: second array dimension
%       3: third array dimension (default)
%       [x y z]: vector indicating direction of B0 field

% Modified June 2013
% Simplified interface and removed VOXEL, FILTER_SIZE and SPHERE_VOLUME
% parameters

% Modified January 2011
% After playing around with including a DC offset to incorporate the
% lorentz sphere correction, I've decided that, for MRI brain imaging, the
% correction is not really necessary. Overall, the correction to the change
% in magnetic B field will be approximately B0/3*chi_water. This offset in
% the derived phase data is eliminated.
%   The parameters ACQMODE and VFACTOR have been removed. The code will now
% only produce a 3D kernel (given that 2D is not a realistic situation),
% and perturbations to the voxel size can be included when supplying the
% VOXEL parameter.

%% VALIDATE PARAMETERS
error(nargchk(2, 3, nargin, 'struct'));

filtersize = datasize*2+1;

if nargin < 3 || isempty(B0dir)
    B0dir = [0 0 1];
elseif length(B0dir) == 1
    switch B0dir
        case 1
            B0dir = [1 0 0];
        case 2
            B0dir = [0 1 0];
        case 3
            B0dir = [0 0 1];
        otherwise
            disp(B0dir)
            error 'Unrecognised B0DIR parameter'
    end
elseif length(B0dir) ~= 3
    error 'Unrecognised B0DIR parameter'
else
    B0dir = B0dir/norm(B0dir);
    B0dir = B0dir(:)';
end

%% CREATE FILTER KERNEL

[x y z] = ndgrid(1:filtersize(1),1:filtersize(2),1:filtersize(3));

x = (x-ceil(size(x,1)/2));
y = (y-ceil(size(y,2)/2));
z = (z-ceil(size(z,3)/2));

% dipole tensor
N = numel(x);
rvec = [x(:) y(:) z(:)];
r = repmat(sqrt(sum(rvec.^2,2)), 1, 3);
rhat = rvec ./ r;
rho = repmat(B0dir, N, 1);

clear x y z

f = B0 ./ (4*pi) .* (3 .* repmat(dot(rho, rhat, 2), 1, 3) .* rhat - rho) ./ r.^3; 
f(isnan(f)) = 0 ;

% return z direction only
f = dot(f,rho,2);
f = reshape(f, filtersize);
