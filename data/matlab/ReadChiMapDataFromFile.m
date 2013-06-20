function [data, B0, bdir] = ReadChiMapDataFromFile(filename)
% [data, B0, bdir] = ReadChiMapDataFromFile(filename)

    % open file
    fid = fopen(filename, 'r');
    if fid == -1
        error 'Could not open file'
    end
    
    % check endianness
    EndianCheckValue = fread(fid, 1, 'uint', 'l');
    if EndianCheckValue == hex2dec('12345678')
        machineformat = 'l';
    elseif EndianCheckValue == hex2dec(78563412)
        machineformat = 'b';
    else
        error 'File error: invalid endian check value'
    end
    
    % check version
    version = fread(fid, 1, 'uint', machineformat);
    if version ~= 1
        error 'File error: unrecognised format version'
    end
    
    % Number of elements
    N = fread(fid, 1, 'double', machineformat);
    
    % Data size
    sz = fread(fid, 3, 'double', machineformat)';
    if N ~= prod(sz)
        error 'File error: number of elements does not equal data size'
    end
    
    % B0
    B0 = fread(fid, 1, 'double', machineformat);
    
    % bdir    
    bdir = fread(fid, 3, 'double', machineformat);
    
    % voxel - this is a legacy field and is no longer used by any code
    voxel = fread(fid, 3, 'double', machineformat);
    
    % read data
    data = fread(fid, 'double', machineformat);
    
    if numel(data) ~= N
        error 'File error: not enough bytes in file'
    end
    
    data = reshape(data, sz);
        
    fclose(fid);
end
