function plot_group_effect(dv, grouping_factor, handle, xlabell, ylabell, xticklabelss,logg, model) 
% dv, dependent variable data
% handle: handle of an existing figure plot or subplot
% xlabell, label for x axis
% ylabell, label for y axis
% optional - xticklabelss, labels for grouping variable on x axis (if empty or not provided, read them from grouping_factor levels directly)
% logg, if 1, y is in log scale, 0 by default
% model is used to remove flagged outliers from the data
% ex of usage: 
% h=subplot(1,4,4); 
% plot_group_effect(data.final_orient, data.meditation, h, 'Meditation group', 'final orientation threshold', '',0, model)
if ~exist('logg','var'); logg=0; end
if ~exist('model','var'); model.exclude = []; end
factor_levels = unique(grouping_factor);
if ~exist('xticklabelss','var')||isempty(xticklabelss); xticklabelss = {char(factor_levels(1)),char(factor_levels(2))}; end

    % exclude outliers
    if ~isempty(model.exclude) 
       dv(model.exclude) = []; 
       grouping_factor(model.exclude) = []; 
    end
    
    nbLevels = numel(factor_levels);
    colors = {'b','r','g','m','c','y'};
    medians=nan(1,nbLevels);
    for i=1:nbLevels
        lev = factor_levels(i);
        plot(handle,i.*ones(numel(dv(grouping_factor==lev))),dv(grouping_factor==lev),'Color',colors{i},'Marker','.','LineStyle','none'); hold on; 
        medians(i) = nanmedian(dv(grouping_factor==lev));
    end
    plot(handle,1:nbLevels,medians,'-k');
    xticks(1:nbLevels);xticklabels(xticklabelss); 
    xlabel(xlabell); ylabel(ylabell); 
    xlim([0.5,nbLevels+0.5]);
    if logg==1 
        set(gca, 'YScale', 'log'); 
    else
        % if nothing negative, set the minimum y to 0 and keep the maximum, otherwise, likely an effect with positive and negative values
        if ~any(dv<0)
            c = ylim(); % retrieve current axis limits
            ylim([0 c(2)]);
        end
    end
end