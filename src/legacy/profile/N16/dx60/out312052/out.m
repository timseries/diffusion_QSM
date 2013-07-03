TmpChiMapOut.precision = 'single';
TmpChiMapOut.fid = fopen('out.bin', 'r', 'b');
if fread(TmpChiMapOut.fid, 1, 'int32') ~= 1
	fclose(TmpChiMapOut.fid);
	TmpChiMapOut.fid = fopen('out.bin', 'r', 'l');
	fread(TmpChiMapOut.fid, 1, 'int32');
end
deltab = reshape(fread(TmpChiMapOut.fid, 288000, TmpChiMapOut.precision), [60 60 80]);
mask = reshape(fread(TmpChiMapOut.fid, 288000, 'uint8'), [60 60 80]);
chi = reshape(fread(TmpChiMapOut.fid, 288000, TmpChiMapOut.precision), [60 60 80]);
fclose(TmpChiMapOut.fid);
clear TmpChiMapOut
