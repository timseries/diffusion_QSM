TmpChiMapOut.precision = 'single';
TmpChiMapOut.fid = fopen('out.bin', 'r', 'b');
if fread(TmpChiMapOut.fid, 1, 'int32') ~= 1
	fclose(TmpChiMapOut.fid);
	TmpChiMapOut.fid = fopen('out.bin', 'r', 'l');
	fread(TmpChiMapOut.fid, 1, 'int32');
end
deltab = reshape(fread(TmpChiMapOut.fid, 16000, TmpChiMapOut.precision), [20 20 40]);
mask = reshape(fread(TmpChiMapOut.fid, 16000, 'uint8'), [20 20 40]);
chi = reshape(fread(TmpChiMapOut.fid, 16000, TmpChiMapOut.precision), [20 20 40]);
fclose(TmpChiMapOut.fid);
clear TmpChiMapOut
