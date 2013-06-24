clear all; close all; clc

%% PARAMETERS - MODIFY THESE
dx = 12;                       % Size of data in first dimension

radius = 2;                    % Radius of cylinders/spheres
N = 3;                         % Number of cylinders/spheres
separation = 5;               % Separation between cylinder and sphere centres
chi_in = 4e-7:-1e-7:1e-7;      % Susceptibility of the cylinders/spheres

B0 = 4.7;                      % Magnetic B field (in Tesla)
B0dir = [0 0 1];               % Magnetic B field direction (currently does [0 0 1] only)
TE = 5e-3;                     % Echo time (in seconds)

% Check parameter settings
if radius * 6 >= dx
    error('Radius is too big for datasize dx. Suggest radius = %0.3f or dx = %d.', dx/6, ceil(radius*6))
end

if radius >= separation
    error('Separation is too small. Must be bigger than radius. Suggest separation = %d.', ceil(radius*2.5));
end

if numel(chi_in) ~= N
    error('N = %d cylinders/spheres have been specified, but only %d susceptibilities have been given in chi_in.', N, numel(chi_in))
end

B0dir = B0dir/norm(B0dir); % Ensure B0dir is a unit vector

%% VARIABLES - CHANGE WITH CAUTION. MAY AFFECT OTHER PARTS OF THE SCRIPT.

dy = dx + (N-1)*separation;    % Size of data in second dimension
dz = dx + separation;          % Size of data in third dimension
datasize = [dx dy dz];         % Size of data
cylaxis = [1 0 0];             % Cylinder axis

% main arrays
chi = zeros(dx,dy,dz);         % True susceptibility map
deltab = zeros(dx,dy,dz);      % Change in B field map
models = zeros(dx,dy,dz,3);    % Model map
bordermask = false(dx,dy,dz);
bordermask(2:end-1,2:end-1,2:end-1) = true;

% define (x,y,z) coordinates for each voxel
[x y z] = ndgrid(-dx/2+0.5:dx/2-0.5, -dx/2+0.5:(dy-dx/2)-0.5, -dx/2+0.5:(dz-dx/2)-0.5);

% create kernels
skernel = chimaps_kernel(datasize, B0, B0dir);
ckernel = chimapc_kernel(datasize, B0, cylaxis, B0dir);

%% CREATE SIMULATION

for n=1:N
    yoffset = (n-1)*separation;
    
    % SPHERE
    spheres(n).chi = zeros(dx,dy,dz);
    r = sqrt(x.^2 + (y - yoffset).^2 + z.^2);
    spheres(n).chi(r < radius) = chi_in(n);
    spheres(n).chi = spheres(n).chi;
    
    % identify boundary voxels
    boundaryvoxels = find(r+1 > radius & r-1 < radius);
    boundarymask = false(dx,dy,dz);
    boundarymask(boundaryvoxels) = true;    
    
    % correct boundary voxels
    fun = @(x,y,z) (sqrt(x.*x+y.*y+z.*z)<=radius).*chi_in(n);
    chib = zeros(numel(boundaryvoxels),1);
    for nb = 1:numel(boundaryvoxels)
        nc = boundaryvoxels(nb);
        chib(nb) = triplequad(fun, x(nc)-0.5, x(nc)+0.5, y(nc)-yoffset-0.5, y(nc)-yoffset+0.5, z(nc)-0.5, z(nc)+0.5);
    end
    spheres(n).chi(boundaryvoxels) = chib;
    
    % create a zero border around image
    spheres(n).chi = spheres(n).chi .* bordermask;
    
    % calculate deltab
    spheres(n).deltab = convn(spheres(n).chi, skernel, 'same') .* bordermask;
    
    % create mask
    spheres(n).mask = spheres(n).chi == max(abs(spheres(n).chi(:)));
    
    % Calculate cylinder
    cylinders(n).chi = zeros(dx,dy,dz);
    r = sqrt((y-yoffset).^2 + (z-separation).^2);
    cylinders(n).chi(r < radius) = chi_in(n);
    cylinders(n).chi = cylinders(n).chi;
    
    % identify boundary voxels
    boundaryvoxels = find(r+1 > radius & r-1 < radius);
    boundarymask = false(dx,dy,dz);
    boundarymask(boundaryvoxels) = true;    
    
    % correct boundary voxels
    %fun = @(x,y,z) (sqrt(y.*y+z.*z)<=radius).*chi_in(n);
    fun = @(y,z) (sqrt(y.*y+z.*z)<=radius).*chi_in(n);
    chib = zeros(numel(boundaryvoxels),1);
    for nb = 1:numel(boundaryvoxels)
        nc = boundaryvoxels(nb);
        %chib(nb) = triplequad(fun, x(nc)-0.5, x(nc)+0.5, y(nc)-yoffset-0.5, y(nc)-yoffset+0.5, z(nc)-separation-0.5, z(nc)-separation+0.5);
        chib(nb) = dblquad(fun, y(nc)-yoffset-0.5, y(nc)-yoffset+0.5, z(nc)-separation-0.5, z(nc)-separation+0.5);
    end
    cylinders(n).chi(boundaryvoxels) = chib;
    
    % create a zero border around image
    cylinders(n).chi = cylinders(n).chi .* bordermask;
    
    % calculate deltab
    cylinders(n).deltab = convn(cylinders(n).chi, ckernel, 'same') .* bordermask;

    % create mask
    cylinders(n).mask = cylinders(n).chi == max(abs(cylinders(n).chi(:)));
    
    % update main arrays
    chi = chi + spheres(n).chi + cylinders(n).chi;
    deltab = deltab + spheres(n).deltab + cylinders(n).deltab;
    
    models(:,:,:,1) = models(:,:,:,1) + cylinders(n).mask.*cylaxis(1);
    models(:,:,:,2) = models(:,:,:,2) + cylinders(n).mask.*cylaxis(2);
    models(:,:,:,3) = models(:,:,:,3) + cylinders(n).mask.*cylaxis(3);
    

end

%% Write data to file, save data
WriteChiMapDataToFile('data.bin', deltab, B0, B0dir, [1 1 1])
%WriteChiMapArrayToFile('mask.bin', bordervoxels)
WriteChiMapArrayToFile('mask.bin', bordermask)
WriteChiMapArrayToFile('models.bin', models)

save simulation.mat