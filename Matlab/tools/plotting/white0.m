function newmap = white0(color_cell,middle,m,lims)

%white0 adapted from BLUEWHITERED
%(http://www.mathworks.ch/matlabcentral/fileexchange/4058-bluewhitered)

% Output:

% Your very own, beautiful, custom colormap

% Inputs:

% color_cell - specify 2 or 3 colors on each side of the white 0, input
% cell array of 4 or 6 colors 
%(ex: colors = {[0.75 0 0] [1 0.5 0] [0 0.5 1] [0 0 1]}) (blue - white -
%orange - red)

% middle - default = [1 1 1] (white), color of zero of colormap

% m - length of colormap


%BLUEWHITERED   Blue, white, and red color map.
%   BLUEWHITERED(M) returns an M-by-3 matrix containing a blue to white
%   to red colormap, with white corresponding to the CAXIS value closest
%   to zero.  This colormap is most useful for images and surface plots
%   with positive and negative values.  BLUEWHITERED, by itself, is the
%   same length as the current colormap.
%
%   Examples:
%   ------------------------------
%   figure
%   imagesc(peaks(250));
%   colormap(bluewhitered(256)), colorbar
%
%   figure
%   imagesc(peaks(250), [0 8])
%   colormap(bluewhitered), colorbar
%
%   figure
%   imagesc(peaks(250), [-6 0])
%   colormap(bluewhitered), colorbar
%
%   figure
%   surf(peaks)
%   colormap(bluewhitered)
%   axis tight
%
%   See also HSV, HOT, COOL, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

if nargin < 2
    middle = [1 1 1];
end

if nargin < 3
    m = size(get(gcf,'colormap'),1);
end

if length(color_cell) == 6
    bottom = color_cell{1};
    bottom2 = color_cell{2};
    botmiddle = color_cell{3};
    topmiddle = color_cell{4};
    top2 = color_cell{5};
    top = color_cell{6};
elseif length(color_cell) == 8
    bottom = color_cell{1};
    bottom2 = color_cell{2};
    bottom3 = color_cell{3};
    botmiddle = color_cell{4};
    topmiddle = color_cell{5};
    top3 = color_cell{6};
    top2 = color_cell{7};
    top = color_cell{8};
elseif length(color_cell) == 10
    bottom = color_cell{1};
    bottom2 = color_cell{2};
    bottom3 = color_cell{3};
    bottom4 = color_cell{4};
    botmiddle = color_cell{5};
    topmiddle = color_cell{6};
    top4 = color_cell{7};
    top3 = color_cell{8};
    top2 = color_cell{9};
    top = color_cell{10};
elseif length(color_cell) == 12
    bottom = color_cell{1};
    bottom2 = color_cell{2};
    bottom3 = color_cell{3};
    bottom4 = color_cell{4};
    bottom5 = color_cell{5};
    botmiddle = color_cell{6};
    topmiddle = color_cell{7};
    top5 = color_cell{8};
    top4 = color_cell{9};
    top3 = color_cell{10};
    top2 = color_cell{11};
    top = color_cell{12};
else
    bottom = color_cell{1};
    botmiddle = color_cell{2};
    topmiddle = color_cell{3};
    top = color_cell{4};
end

% Find middle
if nargin < 4
    lims = get(gca, 'CLim');
end

% Find ratio of negative to positive
if (lims(1) < 0) & (lims(2) > 0)
    % It has both negative and positive
    % Find ratio of negative to positive
    ratio = abs(lims(1)) / (abs(lims(1)) + lims(2));
    neglen = round(m*ratio);
    poslen = m - neglen;
    
    % Just negative
    if length(color_cell) == 6
        new = [bottom; bottom2; botmiddle; middle];
    elseif length(color_cell) == 8
        new = [bottom; bottom2; bottom3; botmiddle; middle];
    elseif length(color_cell) == 10
        new = [bottom; bottom2; bottom3; bottom4; botmiddle; middle];
    elseif length(color_cell) == 12
        new = [bottom; bottom2; bottom3; bottom4; bottom5; botmiddle; middle];
    else
        new = [bottom; botmiddle; middle];
    end
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, neglen);
    newmap1 = zeros(neglen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap1(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % Just positive
    if length(color_cell) == 6
        new = [middle; topmiddle; top2; top];
    elseif length(color_cell) == 8
        new = [middle; topmiddle; top3; top2; top];
    elseif length(color_cell) == 10
        new = [middle; topmiddle; top4; top3; top2; top];
    elseif length(color_cell) == 12
        new = [middle; topmiddle; top5; top4; top3; top2; top];
    else
        new = [middle; topmiddle; top];
    end
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, poslen);
    newmap = zeros(poslen, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
    % And put 'em together
    newmap = [newmap1; newmap];
    
elseif lims(1) >= 0
    % Just positive
    if length(color_cell) == 6
        new = [middle; topmiddle; top2; top];
    elseif length(color_cell) == 8
        new = [middle; topmiddle; top3; top2; top];
    elseif length(color_cell) == 10
        new = [middle; topmiddle; top4; top3; top2; top];
    elseif length(color_cell) == 12
        new = [middle; topmiddle; top5; top4; top3; top2; top];
    else
        new = [middle; topmiddle; top];
    end
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
else
    % Just negative
    if length(color_cell) == 6
        new = [bottom; bottom2; botmiddle; middle];
    elseif length(color_cell) == 8
        new = [middle; topmiddle; top3; top2; top];
    else
        new = [bottom; botmiddle; middle];
    end
    len = length(new);
    oldsteps = linspace(0, 1, len);
    newsteps = linspace(0, 1, m);
    newmap = zeros(m, 3);
    
    for i=1:3
        % Interpolate over RGB spaces of colormap
        newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
    end
    
end

% 
% m = 64;
% new = [bottom; botmiddle; middle; topmiddle; top];
% % x = 1:m;
% 
% oldsteps = linspace(0, 1, 5);
% newsteps = linspace(0, 1, m);
% newmap = zeros(m, 3);
% 
% for i=1:3
%     % Interpolate over RGB spaces of colormap
%     newmap(:,i) = min(max(interp1(oldsteps, new(:,i), newsteps)', 0), 1);
% end
% 
% % set(gcf, 'colormap', newmap), colorbar