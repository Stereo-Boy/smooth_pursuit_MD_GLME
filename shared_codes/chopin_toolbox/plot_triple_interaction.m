function plot_triple_interaction(dv, x_grouping_factor, legend_grouping_factor, plot_grouping_factor, xlabell, ylabell, legendLabels, titleLabels, model) 
% Plot triple interaction effect for a dependent variable between a covariate (in x axis) and a legend grouping factor (put in legend), and two separate plots according to a plot grouping factor
% The covariate can also be categorical but the grouping factor has to have two levels only.
% dv, dependent variable data
% x_grouping_factor data (2 levels only) to be shown on x axis
% legend_grouping_factor data (2 levels only): this will shown in legend
% plot_grouping_factor (2 levels only): another grouping factor, to be plot on separate plots
% handle: handle of an existing figure plot or subplot
% xlabell, label for x axis
% ylabell, label for y axis
% optional legendLabels, label for grouping variable in legend, if not provided, read if from legend_grouping_factor values
% optional titleLabels, title labels for each plot, if not provided, read if from plot_grouping_factor values
% the model structure with the field exclude
% ex of usage: 
% plot_triple_interaction(data_Hands.Speed_Metric_T_Hands, data_Hands.ageGroup, data_Hands.stereo, data_Hands.hand, 'Age group', 'Smoothness (speed metrics)', [], [], model) 

if ~exist('model','var'); model.exclude = []; end
legend_levels = unique(legend_grouping_factor);
plot_levels = unique(plot_grouping_factor);
if ~exist('legendLabels','var')||isempty(legendLabels); legendLabels = {char(legend_levels(1)),char(legend_levels(2))}; end
if ~exist('titleLabels','var')||isempty(titleLabels); titleLabels = {char(plot_levels(1)),char(plot_levels(2))}; end

try
    % exclude outliers
    if ~isempty(model.exclude) 
       dv(model.exclude) = []; 
       x_grouping_factor(model.exclude) = []; 
       legend_grouping_factor(model.exclude) = []; 
       plot_grouping_factor(model.exclude) = []; 
    end
    figure('Color', 'w', 'units','normalized','outerposition',[0 0.1 1 0.5]);
    % first plot
    subplot(1,2,1);
    plotSubplot(1);
    
    % second plot
    subplot(1,2,2);
    plotSubplot(2);
catch err
    keyboard
end
    function plotSubplot(plotlevel)
        % plot for legend_grouping_factor = first value
        p1=plotLegend(plotlevel,1);

        % plot for legend_grouping_factor = second value
        p2=plotLegend(plotlevel,2);
        title(titleLabels(plotlevel))
        xlabel(xlabell); ylabel(ylabell);
        % plot legend
        legend([p1,p2],legendLabels);
    end
    function p=plotLegend(plotlevel,legendlevel)
        data2plot = (legend_grouping_factor==legend_levels(legendlevel))&(plot_grouping_factor==plot_levels(plotlevel));
        x = x_grouping_factor(data2plot); 
        y =  dv(data2plot);
        xlevels = unique(x);
        if legendlevel==1
            p=plot(x,y,'b.'); hold on; 
        else
            p=plot(x,y,'r.'); hold on; 
        end
        medians(1) = nanmedian(y(x==xlevels(1))); medians(2) = nanmedian(y(x==xlevels(2)));
        if legendlevel==1
            plot(xlevels,medians,'b-'); 
        else
            plot(xlevels,medians,'r-'); 
        end
        % if nothing negative, set the minimum y to 0 and keep the maximum, otherwise, likely an effect with positive and negative values
        if ~any(dv<0)
            c = ylim(); % retrieve current axis limits
            ylim([0 c(2)]);
        end
    end
    
end