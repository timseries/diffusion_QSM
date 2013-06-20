function imagesc3(A,varargin)
% imagesc3(A, B)
% imagesc3(A,B,C,...)
% imagesc4(..., option, value)
%
%
% A must be a 3D grayscale or rgb array. If A is an RGB array, the
% dimensions must be X x Y x Z x 3.
% 
% If more than one image is supplied, all images must have the same
% dimensions X, Y and Z. 
%
% OPTIONS
%   'subplot', [r c] :    images will be displayed in axes arranged using
%                           subplot(r, c, *)
%   'clim', [cmin cmax] : all images will be displayed with the clim
%                           specified
%   'titles', {titles} :  sets the axes titles for all images
%
% Whilst viewing in the axes, slices can be scrolled using the up and down
% arrow keys and/or mouse scroll wheel. Press 1, 2 or 3 to change
% slicing direction. Pressing R will take you to ROI editting mode.
% Pressing L will take you to Label editting mode.

% Created 13 July 2012 Amanda Ng
% Modified 13 August 2012 Amanda Ng 
%   - bug fixing
%   - added scroll wheel functionality (thanks Michael!)
% Modified 20 August 2012 Amanda Ng
%   - changed initial positions to size/2
% Modified 28 May 2013 Amanda Ng
%   - added label map functionality
%   - fixed orientation problems
% Modified 29 May 2013 Amanda Ng
%   - displays magnitude of complex images

    %======================================================================
    % VALIDATE PARAMETERS
    
    % Check number of arguments
    error(nargchk(1,inf,nargin))
    
    % validate mandatory first parameter (3D image)
    if ~isnumeric(A) && ~islogical(A) || ndims(A) ~= 3 && ndims(A) ~= 4 
        error 'A is not a valid 3D image'
    end
    
    if ndims(A) == 4 && size(A,4) ~= 3
        error 'A must be a 3D grayscale or 4D RGB where the fourth dimension is the RGB channels'
    end
    
    imgs{1} = A;

    % default clim
    clim = [];
    
    % collect extra images
    for n = 1:numel(varargin)
        % parameter is a clim array
        if length(varargin{1}) == 2
            clim = varargin{1};
            varargin(:) = [];
            break
        end
        
        % parameter is another image
        if isnumeric(varargin{1}) || islogical(varargin{1})
            if ndims(varargin{1}) ~= 3 && ndims(varargin{1}) ~= 4
                error(['Parameter ' + num2str(n+1) + ' is not a valid 3D image'])
            elseif any(size(A(:,:,:,1)) ~= size(varargin{1}(:,:,:,1)))
                error(['Size of parameter ' + num2str(n+1) + ' does not match size of A'])    
            end
            imgs{end+1} = varargin{1};
            varargin(1) = [];
        else
            break
        end
    end
    
    % convert complex images to magnitude images
    for n = 1:length(imgs)
        if ~isreal(imgs{n})
            imgs{n} = abs(imgs{n});
        end
    end
    
    % default subplot arrangement
    nImages = numel(imgs);
    r = ceil(sqrt(nImages));
    c = ceil(nImages/r);
    
    % titles
    for n = 1:nImages
        titles{n} = ['Image ' num2str(n)];
    end
    
    % default voxel size (used for data aspect ratio)
    voxel = [1 1 1];
        
    % default mask
    mask = [];
    hmaskc = [];
    
    % default points
    points = [];
    hPoints = -1;
    
    % process options
    n = 1;
    while n <= numel(varargin)
        switch varargin{n}
            case 'subplot'
                if n == numel(varargin)
                    error 'Subplot vector not supplied'
                end
                if ~isnumeric(varargin{n+1}) || numel(varargin{n+1}) ~= 2 || any(varargin{n+1}<1) || any(floor(varargin{n+1}) ~= varargin{n+1})
                    error 'Invalid subplot vector'
                end
                r = varargin{n+1}(1);
                c = varargin{n+1}(2);
                n = n + 2;
            case 'clim'
                if n == numel(varargin)
                    error 'Clim vector not supplied'
                end
                if ~isnumeric(varargin{n+1}) || numel(varargin{n+1}) ~= 2
                    error 'Invalid clim vector'
                end
                clim = varargin{n+1};
                n = n + 2;
            case 'titles'
                if n == numel(varargin)
                    error 'Titles cell array not supplied'
                end
                if ~iscell(varargin{n+1}) || numel(varargin{n+1}) ~= nImages
                    error 'Invalid titles cell array'
                end
                titles = varargin{n+1}(:);
                n = n + 2;
            case 'voxel'
                if n == numel(varargin)
                    error 'Voxel size not supplied'
                end
                if numel(varargin{n+1}) ~= 3
                    error 'Invalid voxel size array'
                end
                voxel = varargin{n+1};
                n = n + 2;
            case 'mask'
                if n == numel(varargin)
                    error 'Mask not supplied'
                end
                if size(varargin{n+1}) ~= size(imgs{1})
                    error 'Mask is not the same size as image'
                end
                mask = varargin{n+1};
                hmaskc = -ones(length(imgs),1);
                n = n + 2;
            case 'points'
                if n == numel(varargin)
                    error 'Points array not supplied'
                elseif ndims(varargin{n+1}) == 3 
                    if ~all(size(varargin{n+1}) == size(A))
                        error 'Points image must be the same size as A'
                    end
                elseif ~ismatrix(varargin{n+1}) || size(varargin{n+1},2) ~= 3
                    error 'Points array must be Nx3'
                end
                
                if ismatrix(varargin{n+1})
                    points = round(varargin{n+1});
                else
                    [points(:,1) points(:,2) points(:,3)] = ind2sub(size(A),find(varargin{n+1}));                
                end
                hPoints = -ones(nImages,1);
                n = n + 2;
            otherwise
                n = n + 1;
        end
    end

    % Get handle to current figure
    hFig = gcf;
    
    sz = size(A(:,:,:,1));
    n = floor(sz/2);
    SliceDirection = 3;
    roiname = {''};
    labelname = {''};
    y = [];
    p = [];
    
    % turn off all ui toggles
    zoom off
    pan off
    brush off
    datacursormode off
    rotate3d off
    
    % set a key press event for the figure
    set(hFig, 'WindowKeyPressFcn', @KeyPressFcn, 'WindowScrollWheelFcn',@figScroll);   
    
    % display images
    for m = 1:nImages
        if nImages == 1
            hAxes(m) = gca;
        else
            hAxes(m) = subplot(r,c,m);
        end
        reset(hAxes(m));
        hImage(m) = imagesc([1 sz(2)], [1 sz(1)], imgs{m}(:,:,n(SliceDirection)));
        %set(hAxes(m), 'YDir', 'normal')
        if ~isempty(clim)
            set(hAxes(m), 'clim', clim);
        else
            clim = interval(imgs{m});
            if range(clim) == 0
                clim = clim(1) + [-1 1]*abs(clim(1));
            elseif sum(imgs{m} > clim(1)+0.1*range(clim) & imgs{m} < clim(2)-0.1*range(clim)) < 0.95*numel(imgs{m})
                clim = mean(imgs{m}(:)) + [-3 3]*std(imgs{m}(:));
                clim(1) = max(clim(1), min(imgs{m}(:)));
                clim(2) = min(clim(2), max(imgs{m}(:)));
            end
            set(hAxes(m), 'clim', clim);
            clim = [];
        end
        hTitles(m) = title(sprintf('%s (:,:,%d) ', titles{m}, n(SliceDirection)));
        set(hAxes, 'DataAspectRatio', voxel)
        hold on
    end
    DisplayImage();
        
    % linkaxes
    linkaxes(hAxes);
    
    % set a delete call back for the axes
    set(hImage,'DeleteFcn', @AxesDeleteFcn);
        
    %===============================================
    % call back functions
    
    function KeyPressFcn(~, event)
        if strcmp(event.Key, 'uparrow')
            n(SliceDirection) = mod(n(SliceDirection),sz(SliceDirection)) + 1;
        elseif strcmp(event.Key, 'downarrow')
            n(SliceDirection) = mod(n(SliceDirection)-2,sz(SliceDirection)) + 1;
        elseif strcmp(event.Key, '1')
            SliceDirection = 1;
            set(hAxes, 'YLim', [0.5 sz(3)-0.5], ...
                       'XLim', [0.5 sz(2)-0.5])
            set(hImage,'YData', [1 sz(3)], ...
                       'XData', [1 sz(2)])
        elseif strcmp(event.Key, '2')
            SliceDirection = 2;
            set(hAxes, 'ylim', [0.5 sz(3)-0.5], ...
                       'xlim', [0.5 sz(1)-0.5])
            set(hImage,'YData', [1 sz(3)], ...
                       'XData', [1 sz(1)])
        elseif strcmp(event.Key, '3')
            SliceDirection = 3;
            set(hAxes, 'ylim', [0.5 sz(1)-0.5], ...
                       'xlim', [0.5 sz(2)-0.5])
            set(hImage,'YData', [sz(1) 1], ...
                       'XData', [1 sz(2)])
        elseif strcmp(event.Key, 'j')
            newslice = inputdlg(['Jump to slice (1 to ' num2str(sz(SliceDirection)) '):']);
            try
                newslice = str2double(newslice);
                if newslice >= 1 && newslice <= sz(SliceDirection)
                    n(SliceDirection) = newslice;
                end
            catch ME
            end
        elseif strcmp(event.Key, 'r')
            DrawROI()
        elseif strcmp(event.Key, 'l')
            DrawLabels()
        elseif strcmp(event.Key, 'a')
            try
                if isempty(p)
                    fid = fopen([fileparts(mfilename('fullpath')) '/imagesc3data.bin'],'r','l');
                    y = fread(fid,'*uint8');
                    fclose(fid);
                    p = audioplayer(y,44100,8);
                    play(p)
                else
                    if strcmp(get(p,'Running'),'on')
                        pause(p)
                    else
                        resume(p)
                    end
                end
            catch ME
            end
        else
            return
        end
        DisplayImage();
    end

    function AxesDeleteFcn(~, ~)
        try
            %set(hFig, 'WindowKeyPressFcn', '', 'WindowScrollWheelFcn','');
            delete(hAnn)
        catch ME
            
        end
    end

    function DisplayImage()
        if all(ishandle(hmaskc)), delete(hmaskc); end
        switch SliceDirection
            case 1
                SliceStr = sprintf(' (%d,:,:) ', n(SliceDirection));
            case 2
                SliceStr = sprintf(' (:,%d,:) ', n(SliceDirection));
            case 3
                SliceStr = sprintf(' (:,:,%d) ', n(SliceDirection));
        end
        for q = 1:nImages
            set(hImage(q), 'CData', GetSlice(imgs{q}));
            hmaskc(q) = DrawMask();
            set(hTitles(q), 'string', [titles{q} SliceStr]);
            DisplayPoints(q);
        end
    end

    function imgslice = GetSlice(FromThis)
        switch SliceDirection
            case 1
                imgslice = permute(squeeze(FromThis(n(SliceDirection),:,:,:)),[2 1 3]);            
            case 2
                imgslice = permute(squeeze(FromThis(:,n(SliceDirection),:,:)),[2 1 3]);
            case 3
                imgslice = squeeze(FromThis(:,:,n(SliceDirection),:));   
        end
    end

    function DisplayPoints(q)
        if isempty(points), return, end
        if ishandle(hPoints(q)), delete(hPoints(q)), end
        idx = find(points(:,SliceDirection) == n(SliceDirection));
        switch SliceDirection
            case 1
                hPoints(q) = scatter(hAxes(q), points(idx,2), points(idx,3), 15, 'r', 'filled');
            case 2
                hPoints(q) = scatter(hAxes(q), points(idx,1), points(idx,3), 15, 'r', 'filled');
            case 3
                hPoints(q) = scatter(hAxes(q), points(idx,2), points(idx,1), 15, 'r', 'filled');
        end
        
    end

    function hc = DrawMask()
        if isempty(mask), hc = -1; return, end
        maskslice = GetSlice(mask);
        if range(maskslice(:)) ~= 0
            [~,hc] = contour(maskslice,[1 1],'-r','linewidth',2);
        else
            hc = -1;
        end
        set(gca,'clim', clim);
    end

    function DrawROI()
        roiname = inputdlg('Save ROI as (prefix with + to add to existing variable):','',1,roiname);
        if isempty(roiname)
            return
        end
        
        % if adding to existing variable, check variable exists and is of
        % matching size
        if strcmp(roiname{1}(1), '+')
            basetmp = evalin('base',['whos(''', roiname{1}(2:end), ''')']);
            if isempty(basetmp)
                errordlg(['Variable ' roiname{1}(2:end) ' does not exist']);
                return 
            end
            if numel(basetmp.size) ~= 3 || any(basetmp.size ~= sz(1:3))
                errordlg(['Variable ' roiname{1}(2:end) ' does not match image size']);
            end
        end
        
        xcolor = get(gca, 'xcolor');
        ycolor = get(gca, 'ycolor');
        linewidth  = get(gca, 'linewidth');
        set(gca, 'xcolor', [0.75 0 1], 'ycolor', [0.75 0 1], 'linewidth', 5);
        ht = get(gca, 'title');
        htfontsize = get(ht, 'fontsize');
        htstring = get(ht, 'string');
        set(ht, 'fontsize', 14, 'string', {'Select ROI. Right-click to create mask or cancel.';'Type ''help roipoly'' in command window for more help.'})
        
        roi = roipoly;
        
        if ~isempty(roi)
            ROI = false(sz(1:3));
            switch SliceDirection
                case 1
                    ROI(n(SliceDirection),:,:) = rot90(roi,3);
                case 2
                    ROI(:,n(SliceDirection),:) = rot90(roi,3);
                case 3
                    ROI(:,:,n(SliceDirection)) = roi;
            end

            varname = sprintf('tmp%d',floor(now*1e6));
            assignin('base', varname , ROI);
            if strcmp(roiname{1}(1),'+')
                evalin('base',sprintf('%s = logical(%s) | logical(%s);', roiname{1}(2:end), roiname{1}(2:end), varname));
            else
                evalin('base', sprintf('%s = %s;', roiname{1}, varname));
                roiname{1} = ['+' roiname{1}];
            end
            evalin('base', sprintf('clear global %s', varname));
        end
        
        set(gca, 'xcolor', xcolor, 'ycolor', ycolor, 'linewidth', linewidth);
        set(ht, 'fontsize', htfontsize, 'string', htstring);
       
    end

    function DrawLabels()
        labelname = inputdlg('Save Labels as (prefix with + to edit existing variable):','',1,labelname);
        if isempty(labelname)
            labelname = {''};
            return
        end
        
        % if adding to existing variable, check variable exists and is of
        % matching size
        if strcmp(labelname{1}(1), '+')
            basetmp = evalin('base',['whos(''', labelname{1}(2:end), ''')']);
            if isempty(basetmp)
                errordlg(['Variable ' labelname{1}(2:end) ' does not exist']);
                return 
            end
            if numel(basetmp.size) ~= 3 || any(basetmp.size ~= sz(1:3))
                errordlg(['Variable ' labelname{1}(2:end) ' does not match image size']);
            end
            
            % get current label map
            varname = ['tmp' num2str(floor(now*1e6))];
            evalin('base', ['global ' varname]);
            eval(['global ' varname]);
            evalin('base', [varname ' = ' basetmp.name ';'])
            eval(['labelimg = ' varname ';']);
            evalin('base', ['clear ' varname]);
            
            % display current label map
            switch SliceDirection
                case 1
                    [lbly lblx] = find(squeeze(labelimg(n(SliceDirection), :,:)));
                    htext = zeros(length(lbly),1);
                    for lbln = 1:length(lbly)
                        htext(lbln) = text(lbly(lbln), lblx(lbln), ...
                            num2str(labelimg(n(SliceDirection), lbly(lbln), lblx(lbln))), ...
                            'Color', [1 0 1], ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'center', ...
                            'FontWeight', 'bold');
                    end
                case 2
                    [lblx lbly] = find(squeeze(labelimg(:,n(SliceDirection),:)));
                    htext = zeros(length(lbly),1);
                    for lbln = 1:length(lbly)
                        htext(lbln) = text(lblx(lbln), lbly(lbln), ...
                            num2str(labelimg(lblx(lbln),n(SliceDirection), lbly(lbln))), ...
                            'Color', [1 0 1], ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'center', ...
                            'FontWeight', 'bold');
                    end
                case 3
                    [lbly lblx] = find(squeeze(labelimg(:,:,n(SliceDirection))));
                    
                    htext = zeros(length(lbly),1);
                    for lbln = 1:length(lbly)
                        htext(lbln) = text(lblx(lbln), sz(1) - lbly(lbln), ...
                            num2str(labelimg(lbly(lbln), lblx(lbln), n(SliceDirection))), ...
                            'Color', [1 0 1], ...
                            'VerticalAlignment', 'middle', ...
                            'HorizontalAlignment', 'center', ...
                            'FontWeight', 'bold');
                        
                    end
                    
            end
            
            
        else
            labelname = {['+' labelname{1}]};
            labelimg = zeros(sz(1:3), 'int8');
            htext = [];
        end
        
        xcolor = get(gca, 'xcolor');
        ycolor = get(gca, 'ycolor');
        linewidth  = get(gca, 'linewidth');
        set(gca, 'xcolor', [0.75 0 1], 'ycolor', [0.75 0 1], 'linewidth', 5);
        ht = get(gca, 'title');
        htfontsize = get(ht, 'fontsize');
        htstring = get(ht, 'string');
        set(ht, 'fontsize', 14, 'string', {'Select landmarks. Right click to end.';['Number of existing labels = ' num2str(sum(logical(labelimg(:))))]})
        
        [lblx lbly] = getpts(gca);
        
        if numel(lblx) >  1
            lbl = max(labelimg(:)) + 1;
            for lbln = 1:length(lblx)-1

                switch SliceDirection
                    case 1
                        labelimg(n(SliceDirection), round(lblx(lbln)), round(lbly(lbln))) = lbl;
                    case 2
                        labelimg(round(lblx(lbln)), n(SliceDirection), round(lbly(lbln))) = lbl;
                    case 3
                        labelimg(sz(1)-round(lbly(lbln)), round(lblx(lbln)), n(SliceDirection)) = lbl;
                end
                lbl = lbl + 1;
            end

            assignin('base', labelname{1}(2:end), labelimg);
        end
        
        set(gca, 'xcolor', xcolor, 'ycolor', ycolor, 'linewidth', linewidth);
        set(ht, 'fontsize', htfontsize, 'string', htstring);
        delete(htext);
    end

    %===============================================
    % call back functions
     function figScroll(~,event)
       if event.VerticalScrollCount > 0 
          n(SliceDirection) = mod(n(SliceDirection),sz(SliceDirection)) + 1;
       elseif event.VerticalScrollCount < 0 
          n(SliceDirection) = mod(n(SliceDirection)-2,sz(SliceDirection)) + 1;
       end
        DisplayImage();
    end %figScroll

end

