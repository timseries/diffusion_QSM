function WriteChiMapDataToFile(filename, data, B0, bdir, voxel)
% WriteChiMapDataToFile(filename, data, B0, bdir, voxel)

% Change log
% 24 Oct 2011 - Added endian check value
% 25 Oct 2011 - Added bdir vector

    % write to file
    fid = fopen(filename, 'w');
    
    % write out endian check value
    fwrite(fid, hex2dec('12345678'), 'uint'); % endian check value
    
    % write out version
    fwrite(fid, 1, 'uint');
    
    % write out header
    fwrite(fid, numel(data), 'double');
    fwrite(fid, size(data), 'double');
    fwrite(fid, B0, 'double');
    fwrite(fid, bdir, 'double');
    fwrite(fid, voxel, 'double');

    % write out data
    fwrite(fid, data(:), 'double');
    
    fclose(fid);
end
