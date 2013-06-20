function k = chimapc_kernel(datasize, B0, cylinderaxis, B0dir)
% KERNEL = chimapc_kernel(DATASIZE, B0, CYLINDERAXIS, B0DIR)
%
% Create cylindrical filter kernel
%
% DATASIZE      size of data array
% B0            strength of B0 field (in Tesla)
% CYLINDERAXIS  direction of cylinder axis (3 element vector)
% B0DIR         direction of B0 field (see below, default is 3rd array 
%               dimension (z))
%
% B0DIR
% =====
%   Indicates the direction of the B0 field. It can be set to
%       1: first array dimension
%       2: second array dimension
%       3: third array dimension (default)
%       [x y z]: vector indicating direction of B0 field

% Modified June 2013
% Simplified parameters

error(nargchk(3,4,nargin));

if numel(datasize) ~= 3
    error('DATASIZE must be a 3 element array');
end

if numel(cylinderaxis) ~= 3
    error('CYLINDERAXIS must be a 3 element array');
end

if nargin < 4 || isempty(B0dir)
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

% rename parameters
d = datasize;
b = B0dir;
c = cylinderaxis;

% constants
alpha = 0.47;

% ensure b and c are unit vectors
b = b / norm(b);
c = c / norm(c);

% calculate kernel values
[x y z] = ndgrid(-d(1):d(1), -d(2):d(2), -d(3):d(3));

bc = dot(b,c);
b_cbc = b - c*bc;
rr = x.*x + y.*y + z.*z;
rc = x.*c(1) + y.*c(2) + z.*c(3);
r_r_ = rr-rc.^2;

k = B0/2/pi ./ r_r_ .* (2 * (x.*b_cbc(1) + y.*b_cbc(2) + z.*b_cbc(3)).^2 ./...
    r_r_ - (1-bc^2)) ;
k(isinf(k)) = 0;
k(isnan(k)) = 0;
k(d(1)+1, d(2)+1, d(3)+1) = B0/6 * (3*bc^2-1);

a = alpha;
fun = @(rc) (2*a-abs(rc)).^2 .* (3*a-(2*a-abs(rc))) ./ (4*a^3);
R = fun(rc);
R(abs(rc) > 2*a) = 0;
k = k .* R;

% remove slices with only zero value elements
for n = 1:d(1)
    if nnz(k(1,:,:)) == 0 && nnz(k(end,:,:)) == 0
        k(1,:,:) = [];
        k(end,:,:) = [];
    end
end

for n = 1:d(2)
    if nnz(k(:,1,:)) == 0 && nnz(k(:,end,:)) == 0
        k(:,1,:) = [];
        k(:,end,:) = [];
    end
end

for n = 1:d(3)
    if nnz(k(:,:,1)) == 0 && nnz(k(:,:,end)) == 0
        k(:,:,1) = [];
        k(:,:,end) = [];
    end
end

