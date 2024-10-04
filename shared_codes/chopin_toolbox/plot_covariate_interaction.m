function plot_covariate_interaction(dv, covariate, legend_grouping_factor, handle, xlabell, ylabell, legendLabels, loggX, loggY, mdl, plotModel,model)  
% dv, dependent variable data
% covariate data (continuous)
% legend_grouping_factor - interaction factor (levels)
% handle: handle of an existing figure plot or subplot
% xlabell, label for x axis
% ylabell, label for y axis
% legendLabels, label for grouping variable in legend, if left empty, read it from legend_grouping_factor values
% loggX, if 1, x is in log scale, 0 by default, optional
% loggY, if 1, y is in log scale, 0 by default, optional
% mdl, the mdl structure (useful only if plotting model predictions with plotModel)
% plotModel, 0 or 1, if 1, will use data predictions in mdl structure and plot the model data, optional
% model is used to remove flagged outliers from the data
% ex of usage: 
% h=subplot(1,4,1); 
% plot_covariate_effect(data.initial_orient, data.music, data.ageGroup, h, 'Music practice (hours)', 'initial orientation threshold', [], 0, 0, mdls{1}, 1, model)
if ~exist('model','var'); model.exclude = []; plotModel=0; end
if ~exist('loggX','var'); loggX=0; end
if ~exist('loggY','var'); loggY=0; end
if ~exist('plotModel','var'); plotModel=0; end
if ~exist('mdl','var')||isempty(mdl); plotModel=0; end
legend_levels = unique(legend_grouping_factor);
if ~exist('legendLabels','var')||isempty(legendLabels)
    legendLabels = {};
    for i=1:numel(legend_levels)
        legendLabels(i) = {char(legend_levels(i))}; 
    end
end

try
    x = covariate; y =  dv; colors = 'brgcmk';
    % exclude outliers
    if ~isempty(model.exclude)
       x(model.exclude) = []; 
       y(model.exclude) = []; 
    end
    % select the correct model predictions in case of plotModel
    if plotModel==1 
        if model.glme == 0 %glm
            y2 =  mdl.Fitted.Response;
        else % glme
            y2 =  mdl.fitted;
        end
    end
    % plot data
    counter = 1; legendLabelsList = {};
    for i=1:numel(legend_levels)
        h(counter)=plot(handle,x(legend_grouping_factor==legend_levels(i)),y(legend_grouping_factor==legend_levels(i)),'MarkerEdgeColor',colors(i),'Marker','.','LineStyle','none'); hold on; 
        legendLabelsList(counter) = legendLabels(i);  counter = counter+1;
        if plotModel==1
            covariate_levels = unique(covariate);
            for j=1:numel(covariate_levels)
                medianX(j) = median(x(covariate==covariate_levels(j)&legend_grouping_factor==legend_levels(i)));
                medianY(j) = median(y2(covariate==covariate_levels(j)&legend_grouping_factor==legend_levels(i)));
            end
            h(counter)=plot(handle,medianX,medianY,'Color',colors(i),'Marker','none','LineStyle','-'); hold on;
            legendLabelsList(counter) = {[legendLabels{i},' model']};  counter = counter+1;
        end
        
    end
    xlabel(xlabell); ylabel(ylabell);
    legend(h,legendLabelsList,'Location','South');
    if loggX==1; set(gca, 'XScale', 'log'); end
    if loggY==1; set(gca, 'YScale', 'log'); end

catch err
    keyboard
end