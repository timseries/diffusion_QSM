function WriteChiMapArrayToFile(filename, data)

    % write to file
    fid = fopen(filename, 'w');
    
    % write endian check value
    fwrite(fid, hex2dec('12345678'), 'uint'); % endian check value

    % write out version
    fwrite(fid, 1, 'uint');    
    
    % write out parameters
    fwrite(fid, numel(data), 'double');
    fwrite(fid, size(data), 'double');
    
    if islogical(data)
        fwrite(fid, data(:), 'uint8');
    else
        fwrite(fid, data(:), 'double');
    end
    
    fclose(fid);
end
