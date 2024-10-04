function plot_interaction(dv, grouping_factor, covariate, handle, xlabell, ylabell, legendLabels, mdl, plotModel, model) 
% Plot interaction effect for a dependent variable between a covariate (in x axis) and a grouping factor (put in legend)
% The covariate can also be categorical but the grouping factor has to have two levels only.
% dv, dependent variable data
% grouping_factor data (2 levels only) - will be shown in legend
% covariate data (continuous or categorical) - will be shown in x axis [Warning: this code is untested for a continuous covariate]
% handle: handle of an existing figure plot or subplot
% xlabell, label for x axis
% ylabell, label for y axis
% optional - legendLabels, labels for grouping variable in legend (if empty or not provided, read them from grouping_factor levels directly)
% mdl, the mdl structure (useful only if plotting model predictions with plotModel)
% plotModel, 0 or 1, if 1, will use data predictions in mdl structure and plot the model data
% ex of usage: 
% h=subplot(1,4,4); 
% plot_interaction(data.Time, data.stereo,data.ageGroup, h, 'Age group (younger / older)','Task completion time (sec)', '',mdls{1}, 1, model)  
if ~exist('model','var'); model.exclude = []; end
if ~exist('plotModel','var'); plotModel=0; end
if ~exist('mdl','var')||isempty(mdl); plotModel=0; end
levels = unique(grouping_factor);
if ~exist('legendLabels','var')||isempty(legendLabels); legendLabels = {char(levels(1)),char(levels(2))}; end

try
    % exclude outliers
    if ~isempty(model.exclude) 
       dv(model.exclude) = []; 
       grouping_factor(model.exclude) = []; 
       covariate(model.exclude) = []; 
    end 
    
    % plot for grouping_factor = first value
    x = covariate(grouping_factor==levels(1)); y =  dv(grouping_factor==levels(1));
    xlevels = unique(x);
    p1=plot(handle,x,y,'b.'); hold on; 
    if plotModel
        if model.glme == 0 %glm
            y2 =  mdl.Fitted.Response(grouping_factor==levels(1)); 
        else %glme
            y2 =  mdl.fitted; y2 = y2(grouping_factor==levels(1)); % in two steps to avoid a bug with model object methods
        end
        m1 = plot(x,y2,'bo'); 

        % plot line trends
        ab=robustfit(x,y2); 
        plot(handle,sort(x),ab(2).*sort(x)+ab(1),'b-');
    end

   if numel(xlevels)<=2
        medians(1) = nanmedian(y(x==xlevels(1))); medians(2) = nanmedian(y(x==xlevels(2)));
        plot(handle,xlevels,medians,'b-');
   end
   
    % plot for grouping_factor = second value
    x = covariate(grouping_factor==levels(2)); y =  dv(grouping_factor==levels(2));
    xlevels = unique(x);
    p2=plot(handle,x,y,'r.'); hold on;
    if plotModel
        if model.glme == 0 %glm
            y2 =  mdl.Fitted.Response(grouping_factor==levels(2)); 
        else %glme
            y2 =  mdl.fitted; y2 = y2(grouping_factor==levels(2)); % in two steps to avoid a bug with model object methods
        end
        m2 = plot(x,y2,'ro'); 
        % plot line trends
        ab=robustfit(x,y2); 
        plot(handle,sort(x),ab(2).*sort(x)+ab(1),'r-');
    end
    xlabel(xlabell); ylabel(ylabell);  
    if numel(xlevels)<=2
        xticklabels(xlevels);
        medians(1) = nanmedian(y(x==xlevels(1))); medians(2) = nanmedian(y(x==xlevels(2)));
        plot(handle,xlevels,medians,'r-');
    end
    
    % plot legend
    if plotModel
        legend([p1,p2,m1,m2],{legendLabels{1},legendLabels{2},'Model estimates','Model estimates'},'Location','southeast');  
    else
        legend([p1,p2],legendLabels,'Location','southeast');
    end

    % if nothing negative, set the minimum y to 0 and keep the maximum, otherwise, likely an effect with positive and negative values
    if ~any(dv<0)
        % Set the minimum y to 0 and keep the maximum
        c = ylim(); % retrieve current axis limits
        ylim([0 c(2)]);
    end
catch err 
    % DEBUGGING
    % write rethrow(err) in the command window to know the error
    keyboard
end