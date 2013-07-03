TmpChiMapOut.precision = 'single';
TmpChiMapOut.fid = fopen('out.bin', 'r', 'b');
if fread(TmpChiMapOut.fid, 1, 'int32') ~= 1
	fclose(TmpChiMapOut.fid);
	TmpChiMapOut.fid = fopen('out.bin', 'r', 'l');
	fread(TmpChiMapOut.fid, 1, 'int32');
end
deltab = reshape(fread(TmpChiMapOut.fid, 96000, TmpChiMapOut.precision), [40 40 60]);
mask = reshape(fread(TmpChiMapOut.fid, 96000, 'uint8'), [40 40 60]);
chi = reshape(fread(TmpChiMapOut.fid, 96000, TmpChiMapOut.precision), [40 40 60]);
fclose(TmpChiMapOut.fid);
clear TmpChiMapOut
