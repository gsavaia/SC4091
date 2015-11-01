function customBodePlot(sys, title, color)    
    opts = bodeoptions('cstprefs');
    opts.Title.String = title;
    
    bodeplot(sys, 'blue', opts);
    
    lineHandle = findobj(gcf,'Type','line','-and','Color','blue');
    set(lineHandle,'Color',color);
end