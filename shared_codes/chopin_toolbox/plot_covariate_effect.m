function plot_covariate_effect(dv, covariate, handle, xlabell, ylabell, loggX, loggY, mdl, plotModel,model)  
% dv, dependent variable data
% covariate data (continuous)
% handle: handle of an existing figure plot or subplot
% xlabell, label for x axis
% ylabell, label for y axis
% loggX, if 1, x is in log scale, 0 by default, optional
% loggY, if 1, y is in log scale, 0 by default, optional
% mdl, the mdl structure (useful only if plotting model predictions with plotModel)
% plotModel, 0 or 1, if 1, will use data predictions in mdl structure and plot the model data, optional
% model is used to remove flagged outliers from the data
% ex of usage: 
% h=subplot(1,4,1); 
% plot_covariate_effect(data.initial_orient, data.music, h, 'Music practice (hours)', 'initial orientation threshold', 0, 0, mdls{1},1, model)
if ~exist('model','var'); model.exclude = []; plotModel=0; end
if ~exist('loggX','var'); loggX=0; end
if ~exist('loggY','var'); loggY=0; end
if ~exist('plotModel','var'); plotModel=0; end
if ~exist('mdl','var')||isempty(mdl); plotModel=0; end
try
    x = covariate; y =  dv;
    % exclude outliers
    if ~isempty(model.exclude)
       x(model.exclude) = []; 
       y(model.exclude) = []; 
    end
    plot(handle,x,y,'k.'); hold on; 
    if plotModel==1
        if model.glme == 0 %glm
            y2 =  mdl.Fitted.Response;
        else % glme
            y2 =  mdl.fitted;
        end
    else
        y2 = zeros(size(x)); 
    end
    if plotModel 
        plot(x,y2,'ro'); 
        ab=robustfit(x,y2); 
        plot(handle,sort(x),ab(2).*sort(x)+ab(1),'r-');
        line([x,x]',[y,y2]','Color','r');
    end
    xlabel(xlabell); ylabel(ylabell);
    xlim([min(x).*0.95,max(x).*1.05]);
    if loggX==1; set(gca, 'XScale', 'log'); end
    if loggY==1
        set(gca, 'YScale', 'log')
    else
        % if nothing negative, set the minimum y to 0 and keep the maximum, otherwise, likely an effect with positive and negative values
        if ~any(dv<0)
            c = ylim(); % retrieve current axis limits
            ylim([0 c(2)]);
        end
    end
catch err
    keyboard
end