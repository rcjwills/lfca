function [h] = plot_field_div_replace(lon,lat,field,ctrs,filter,cmap,nan_mask)

edge_darkening = 0.7;

field = squeeze(field);

%colors = {[0.6 0 0] [1 0.5 0] [0 0.6 1] [0 0 0.6]}; %red-blue
%colors = {[0 0 0.6] [0 0.6 1] [1 0.5 0] [0.6 0 0]}; %blue-red
colors = {[0 0 0.55] [0 0.6 1] [1 0.55 0] [0.55 0 0]}; %blue-red
%colors = {[0 0 1] [0 0.5 1] [1 0.5 0] [0.75 0 0]}; %blue-red_old

%colors = {[55 116 183]/256 [107 203 236]/256 [246 148 50]/256 [205 33 41]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [49 123 184]/256 [72 153 199]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [223 57 41]/256 [189 34 39]/256 [146 27 30]/256}; % alternate blue-red
%colors = {[13 87 150]/256 [49 123 184]/256 [72 153 199]/256 [125 185 215]/256 [168 210 228]/256 [210 228 239]/256 [251 230 145]/256 [248 162 52]/256 [242 102 46]/256 [223 57 41]/256 [189 34 39]/256 [146 27 30]/256}; % alternate blue-red

if nargin > 4
    try
        field = conv2(field,filter,'same');
    catch
        if ~strcmp(filter,'none')
            disp('Warning: unable to inerpret input filter, no filtering applied')
        end
    end
end

if nargin > 6
    field(nan_mask) = nan;
end

if nargin < 4
    ctrs = linspace(nanmin(nanmin(field)),nanmax(nanmax(field)),21);
end

field(field<ctrs(1)) = ctrs(1);

lcmap = length(ctrs)-1;
ch = [min(ctrs) max(ctrs)];

[~,h] = contourf(lon,lat',field',ctrs,'linestyle','none'); colorbar; 
%cmap = colormap(jet(lcmap)); 
caxis(ch); 
if nargin > 5
    colormap(cmap);
else
    cmap = colormap(white0(colors,[1 1 1],lcmap));
    %cmap = colormap(diverging0(colors,[1 1 1],lcmap));
    %cmap = colormap(coldwarm_white(lcmap));
end
hold on;
for i = 1:length(cmap)
    hold on; contour(lon,lat',field',[ctrs(i) ctrs(i)],'LineColor',cmap(i,:)*edge_darkening);
end
%if round(length(lon)/length(lat)) == 2
if max(lat) > 60 && max(lon) > 200 && max(lat) < 100 && min(lon) < 5
    pretty_figure(600,250,'none','none',30:30:330,-90:30:90,16,{'','','90°E','','','180','','','90°W','',''},{'','60°S','30°S','EQ','30°N','60°N',''});
elseif max(lon) > 30 && max(lon) < 40 
    pretty_figure(600,250,'none','none',-30:10:50,30:5:70,16,{'30°W','','10°W','','10°E','','30°E'},{'30°N','','40°N','','50°N','','60°N','','70°N'});
else
    pretty_figure(600,250,'none','none','none','none',16);
end