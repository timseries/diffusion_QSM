%% ========================================================================
%  DATA SETUP

%% Create kernels

n = 3;

[x y z] = ndgrid(-n:n,-n:n,-n:n);
r = sqrt(x.*x + y.*y + z.*z);
skernel = r .* (r<=2);
ckernel = r.*r .* (r<=2);

figure;imagesc3(skernel,ckernel);

%% Create masks

sz = [7 7 7];

cylinder_mask = zeros(sz);
cylinder_mask(2:3,2:3,2:3) = 1;

sphere_mask = ~cylinder_mask;

figure;imagesc3(cylinder_mask, sphere_mask)

%% Create susceptibility map (chi)

chi = cylinder_mask*1 + sphere_mask*2;
figure;imagesc3(chi)

%% Create convolution matrices

As = chimaps_fm(skernel, sz, true);
As(:,~sphere_mask(:)) = 0;
Ac = chimaps_fm(ckernel, sz, true);
Ac(:,~cylinder_mask(:)) = 0;

A = As+Ac;

%% Create delta b vector (b), including random noise
b = A*chi(:) + randn(prod(sz),1).*10^-6;

%% Calculate X using convolution matrix A

X = A'*(A*chi(:)-b);
figure;imagesc3(reshape(X,sz))

%% ========================================================================
%  PROOFS

%% Confirm convolution matrix As works the same as convolution function convn
p1 = As*chi(:)-b;
p2 = reshape(convn(chi.*sphere_mask, skernel, 'same'),[],1)-b;

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Confirm split convolution matrices (As and Ac) work the same as single conv matrix (A)
% ie A*x = As*x + Ac*x
p1 = A*chi(:);
p2 = As*chi(:) + Ac*chi(:);

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Confirm split conv matrices and masked b works the same as Ax - b
% ie A*x-b = As*x-b*ms + Ac*x-b*mc
p1 = A*chi(:) - b;
p2 = As*chi(:) - b.*sphere_mask(:) + Ac*chi(:) - b.*cylinder_mask(:);

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Confirm A'(Ax-b) = As'(As*x-b*ms) + Ac'(Ac*x-b*mc)
s = (As*chi(:)-b.*sphere_mask(:));
c = (Ac*chi(:)-b.*cylinder_mask(:));

sc = s + c;

Y = As'*(sc) + Ac'*(sc);

p1 = X;
p2 = Y;

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Confirm As'(y) works the same as masked convolution method: ms*conv(y,ks)
s = reshape(convn(chi.*sphere_mask,skernel,'same'),[],1)-b.*sphere_mask(:);
c = (Ac*chi(:)-b.*cylinder_mask(:));

sc = s + c;

p1 = As'*(sc);
p2 = sphere_mask(:).*reshape(convn(reshape(sc, sz), skernel, 'same'),[],1);

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Confirm ms*conv(y,ks) + Ac'(y) is the same as A'(Ax-b)
S = sphere_mask(:).*reshape(convn(reshape(sc, sz), skernel, 'same'),[],1);
C = Ac'*(sc);
Y = S + C;

p1 = X;
p2 = Y;

clf; hold on
plot([p1 p2]);
norm(p1-p2)

%% Introduce fourier transform convolution

% note that convn is a linear convolution, while ifftn(fftn) is a circular
% convolution. Therefore the ifftn(fftn) solution needs to involve zero
% padding to mimic linear convolution.

n = 11;
s = ifftn(fftn(chi,sz+n).*fftn(skernel,sz+n));

i1 = ceil(size(skernel)/2); 
i2 = n - i1 + 1; 

p1 = reshape(convn(chi,skernel,'same'),[],1);
p2 = reshape(s(i1(1):end-i2(1),i1(2):end-i2(2),i1(3):end-i2(3)),[],1);

clf; hold on
figure;plot([p1 p2]);
norm(p1-p2)

%% ========================================================================
%  FULL FINAL IMPLEMENTATION

%% Full implementation

n = 9;          % kernel size
sz = [19 19 19];   % data size

%kernels
[x y z] = ndgrid(-n:n,-n:n,-n:n);
r = sqrt(x.*x + y.*y + z.*z);
skernel = r .* (r<=2);
ckernel = r.*r .* (r<=2);

%masks
cylinder_mask = zeros(sz);
cylinder_mask(2:3,2:3,2:3) = 1;
sphere_mask = ~cylinder_mask;

%susceptibility map
chi = cylinder_mask*1 + sphere_mask*2;

%convolution matrices
As = chimaps_fm(skernel, sz, true);
As(:,~sphere_mask(:)) = 0;
Ac = chimaps_fm(ckernel, sz, true);
Ac(:,~cylinder_mask(:)) = 0;
A = As+Ac;

%delta b vector (b), including random noise
b = A*chi(:) + randn(prod(sz),1).*10^-6;

%old implementation
X = A'*(A*chi(:)-b);

%optimised implementation
n = size(skernel);
i1 = ceil(size(skernel)/2); 
i2 = n - i1 + 1; 

tmp = ifftn(fftn(chi.*sphere_mask,sz+n).*fftn(skernel,sz+n));
s = reshape(tmp(i1(1):end-i2(1),i1(2):end-i2(2),i1(3):end-i2(3)),[],1)-b.*sphere_mask(:);

%trying the mask afterward
%tmp = ifftn(fftn(chi,sz+n).*fftn(skernel,sz+n));
%s = reshape(tmp(i1(1):end-i2(1),i1(2):end-i2(2),i1(3):end-i2(3)),[],1).*sphere_mask-b.*sphere_mask(:);
c = (Ac*chi(:)-b.*cylinder_mask(:));
sc = s + c;

tmp = ifftn(fftn(reshape(sc,sz),sz+n).*fftn(skernel,sz+n));
S = sphere_mask(:).*reshape(tmp(i1(1):end-i2(1),i1(2):end-i2(2),i1(3):end-i2(3)),[],1);
C = Ac'*(sc);
Y = S + C;

%tim's optimised implementation using cheecky circshift trick
%n = size(skernel);
%
%nShift=floor(n);
%skernel=circshift(skernel,[-nShift(1) -nShift(2) -nShift(3)]);
%tmp = ifftn(fftn(chi.*sphere_mask).*fftn(skernel));
%s = reshape(tmp,[],1)-b.*sphere_mask(:);
%c = (Ac*chi(:)-b.*cylinder_mask(:));
%sc = s + c;
%
%tmp = ifftn(fftn(reshape(sc,sz)).*(fftn(skernel)));
%S = sphere_mask(:).*reshape(tmp,[],1);
%C = Ac'*(sc);
%Y = S + C;
%

%compare the implementation output
p1 = X;
p2 = Y;

clf; hold on
plot([p1 p2]);
norm(p1-p2)
