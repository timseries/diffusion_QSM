function [fm algver] = chimaps_fm(f, datasize, showprogress)
% FM = CHIMAPS_FM(F,DATASIZE,SHOWPROGRESS)
% [FM ALGVER] = CHIMAPS_FM(F,DATASIZE,SHOWPROGRESS)
% ALGVER = CHIMAPS_FM([])
%
% Creates the filter matrix from the filter kernel F for a data array of
% size DATASIZE. 
%
% Set SHOWPROGRESS to 0 to suppress progress information (default is 1).

% ALGVER is the version of the algorithm. 
% Non-thresholded matrix is not sparse, so create full matrix.

    algver = 2;
    if isempty(f)
        fm = algver;
        return
    end
    
    %% validate parameters

    error(nargchk(2,3,nargin,'struct'));
    
    if nargin<3
        showprogress = 1;
    end

    %% initialise variables
    
    % zero thresholding not yet implemented. Need to look into the effect
    % on accuracy versus benefit in terms of sparse matrix computations.
    zero_threshold = 0;
    
    dataels = prod(datasize);
    filtersize = size(f);
    
    %% create indices
    fwidth = [floor((filtersize-1)/2)' floor(filtersize/2)'];
    [i j k] = ndgrid(1:datasize(1),1:datasize(2),1:datasize(3));
    indices = reshape(1:dataels, datasize);
    
    %apply filter zero threshold
    %f(f<zero_threshold) = 0;
    
    %% initialise sparse matrix

    if showprogress
        fprintf(1, 'Number elements: %d\n',  dataels^2);
        fprintf(1, 'Expected memory requirements for FM (approx): %0.1f GB\n',  min(dataels^2, dataels*nnz(f))*8/1024^3);
    end
    %fm = spalloc(dataels, dataels, min(dataels^2, dataels*nnz(f)));
    fm = zeros(dataels);
    
    nchar = 0;
    
    %% fill sparse matrix
    els = 0;
    myt = tic;
    for n = 1:dataels
        tmpd = (indices(max(i(n)-fwidth(1,1),1):min(i(n)+fwidth(1,2), datasize(1)), ...
                        max(j(n)-fwidth(2,1),1):min(j(n)+fwidth(2,2), datasize(2)), ...
                        max(k(n)-fwidth(3,1),1):min(k(n)+fwidth(3,2), datasize(3))));

        tmpf = f(max(2-i(n)+fwidth(1,1),1):min(fwidth(1,1)+1+datasize(1)-i(n),filtersize(1)), ...
                 max(2-j(n)+fwidth(2,1),1):min(fwidth(2,1)+1+datasize(2)-j(n),filtersize(2)), ...
                 max(2-k(n)+fwidth(3,1),1):min(fwidth(3,1)+1+datasize(3)-k(n),filtersize(3)));

        fm(tmpd(:),n) = tmpf(:);

        els = els + nnz(tmpf);
        
        if showprogress && mod(n,20)==0
            fprintf(1,repmat('\b',1,nchar));
            t = toc(myt);
            tl = t/n * (dataels-n);
            nchar = fprintf(1, 'n=%d of %d; els=%0.0f; time elapsed: %d:%02.0f:%02.0f; time left: %d:%02.0f:%02.0f; avg/calc: %0.3fs', ...
                n, dataels, els, ...
                floor(t/3600), floor(mod(t, 3600)/60), floor(mod(mod(t,3600),60)), ...
                floor(tl/3600), floor(mod(tl, 3600)/60), floor(mod(mod(tl,3600),60)), ...
                t/n);
        end
    end

    %% release excess memory
%     if showprogress
%          fprintf(1, '\nReleasing memory ... ',n,dataels);
%     end
%     
%     %fm = 1*fm;
%     fm = full(fm);

%     if showprogress
%          fprintf(1, 'done\n');
%          fprintf(1, 'Number non-zeros memory allocated for in FM: %d\n', nzmax(fm));
%     end

    %% display computation time
    if showprogress
        t = toc(myt);
        fprintf(1, 'Total time taken: %d:%02.0f:%02.0f\n\n', floor(t/3600), floor(mod(t, 3600)/60), floor(mod(mod(t,3600),60)));
    end

