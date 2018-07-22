function [] = pretty_figure(xsize,ysize,xlab,ylab,xtick,ytick,fontsize,xticklab,yticklab,ctick)

if nargin > 9
    if ~strcmp(ctick,'none')
        hc = colorbar;
        set(hc,'ytick',ctick);
        set(hc,'yticklabel',ctick);
    end
end

set(gcf,'position',[0 1000 xsize ysize])

if ~strcmp(xlab,'none')
    xlabel(gca,xlab,'fontsize',fontsize)
end
if ~strcmp(ylab,'none')
    ylabel(gca,ylab,'fontsize',fontsize)
end
try
    set(gca,'xtick',xtick)
end
try
    set(gca,'ytick',ytick)
end
set(gca,'fontsize',fontsize)
%set(gca,'xgrid','on')
%set(gca,'ygrid','on')
set(gca,'tickdir','out')
%set(gca,'box','off')

if nargin > 7
    if ~strcmp(xticklab,'none')
        set(gca,'xticklabel',xticklab)
    end
    if ~strcmp(yticklab,'none')
        set(gca,'yticklabel',yticklab)
    end
end


