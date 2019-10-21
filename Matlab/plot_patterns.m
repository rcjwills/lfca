function [] = plot_patterns(Xdisc,patternsf,ndisc,time,lon,lat,ctrs,option,title_text)

if nargin < 8
    option = 'none';
end

s(1) = length(lon);
s(2) = length(lat);
    
pattern = reshape(patternsf(ndisc,:),s(1),s(2));

pattern(abs(pattern)>1e5) = nan;

if nargin > 5
    if ~strcmp(option,'split')
        figure; subplot(3,1,[1 2]);
    end
    plot_field_div_replace(lon,lat,pattern,ctrs); pcontinents; set(gca,'xlim',[0 357.5]);
    switch option
        case 'Pacific'
            set(gca,'xlim',[100 295]); set(gca,'ylim',[-44 70])
            set(gca,'xtick',120:30:270); set(gca,'xticklabel',{'120°E','','180°','','120°W',''});
        case 'Atlantic'
            set(gca,'ylim',[-45 86]); set(gca,'xlim',[-80 80])
    end
    colorbar off; originalSize = get(gca, 'Position'); hc = colorbar('eastoutside'); 
    if max(ctrs)==0.2
        set(hc,'ytick',-0.2:0.1:0.2)
    end
    %set(gca,'xticklabel',''); set(gca,'yticklabel',''); 
    set(gca,'fontsize',14)
    if strcmp(option,'split')
        set(gcf,'position',[0 1000 600 250]);
    else
        set(gcf,'position',[0 1000 600 450]);
    end
    originalSize(1) = originalSize(1)-0.04;
    set(gca, 'Position', originalSize);
    if nargin > 9
        title(title_text,'fontsize',14)
    end
    if strcmp(option,'split')
        figure
    else
        subplot(3,1,3)
    end
else
    figure;
end    

%year1 = 1900;
%l = size(Xdisc,1);
%time = (1:l)/12+year1;
year1 = floor(time(1));
year2 = ceil(time(end));

plot(time,Xdisc(:,ndisc),'k');

originalSize = get(gca, 'Position');
originalSize(1) = originalSize(1)-0.04;
set(gca, 'Position', originalSize);

if year1 == 0
    set(gca,'xlim',[0 year2])
elseif year1 == 1900
    set(gca,'xlim',[1900 2017])
elseif year1 == 1984
    set(gca,'xlim',[1984 2009])
elseif year1 < 1900 && year1 > 1849
    set(gca,'xlim',[year1 2016])
else
    set(gca,'xlim',[year1 year2])
end
set(gca,'ylim',[floor(min(Xdisc(:,ndisc))) ceil(max(Xdisc(:,ndisc)))])
set(gca,'fontsize',14)
ylabel('Standard Deviations','fontsize',14);
xlabel('Year','fontsize',14);
if max(ctrs) == 0.6
    set(hc,'ytick',-0.6:0.2:0.6)
end
